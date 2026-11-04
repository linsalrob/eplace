#!/usr/bin/env python
"""
NBDL-aware ePLACE command line wrapper.

This module adds a small command that can search BLAST databases built from
CSIRO/NBDL-style FASTA headers such as::

    >NBDL-...|MT-RNR2 [organism=Dascyllus reticulatus] ... |taxid=80951

The normal ePLACE BLAST workflow asks BLAST for staxid/staxids. Custom BLAST
DBs may not return those fields unless the database was built with a taxid map.
Here we additionally request the subject title (stitle) and, when BLAST does not
provide a usable taxid, parse ``|taxid=...`` from the FASTA header.
"""

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .alignment import (
    IQTreeBuilder,
    check_alignment_consistency,
    create_grouped_fasta_with_queries,
    group_hits_by_group_rank,
    process_grouped_alignment_and_tree_parallel,
    process_query_alignment_and_tree_parallel,
)
from .blast_analysis import BlastHit, BlastRunner, FastaReader, _parse_nbdl_custom_header
from .taxonomy import (
    VALID_RANKS,
    generate_classification_summary,
    process_blast_results_for_taxonomy,
    rewrite_blast_hits,
)

logger = logging.getLogger(__name__)

RANK_CHOICES = [rank for rank in VALID_RANKS if rank != "domain"]
SEARCH_OUTFMT_WITH_TITLE = (
    "6 qseqid sseqid pident length qlen slen qstart qend "
    "sstart send evalue bitscore staxid staxids stitle"
)


def _safe_id(value: str) -> str:
    return value.replace("|", "_").replace("/", "_").replace(" ", "_")


def _is_missing_taxid(value: Optional[str]) -> bool:
    return value is None or value.strip() in {"", "0", "N/A", "NA", "-"}


def _parse_blast_results_with_nbdl_titles(blast_output: Path) -> List[BlastHit]:
    """Parse BLAST tabular output and recover taxids from NBDL subject titles.

    The expected output format is ``SEARCH_OUTFMT_WITH_TITLE``. Older BLAST
    result files without ``stitle`` are still accepted, but no NBDL fallback is
    possible for those lines.
    """
    if not blast_output.exists():
        raise FileNotFoundError(f"BLAST output file not found: {blast_output}")

    hits: List[BlastHit] = []
    with open(blast_output, "r") as handle:
        for line_num, line in enumerate(handle, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 12:
                raise ValueError(
                    f"Invalid BLAST output format at line {line_num}: "
                    f"expected at least 12 fields, got {len(fields)}"
                )

            query_id = fields[0]
            subject_id = fields[1]
            percent_identity = float(fields[2])
            alignment_length = int(fields[3])
            query_length = int(fields[4])
            subject_length = int(fields[5])
            query_start = int(fields[6])
            query_end = int(fields[7])
            subject_start = int(fields[8])
            subject_end = int(fields[9])
            evalue = float(fields[10])
            bit_score = float(fields[11])
            staxid = fields[12] if len(fields) > 12 else ""
            staxids = fields[13] if len(fields) > 13 else ""
            subject_title = fields[14] if len(fields) > 14 else ""

            if _is_missing_taxid(staxid) and subject_title:
                parsed_subject_id, _description, parsed_taxid = _parse_nbdl_custom_header(subject_title)
                if parsed_taxid:
                    staxid = parsed_taxid
                    staxids = parsed_taxid if _is_missing_taxid(staxids) else staxids
                    if subject_id in {"", "N/A", "NA"}:
                        subject_id = parsed_subject_id
                    logger.debug(
                        "Recovered taxid %s for subject %s from NBDL header",
                        parsed_taxid,
                        subject_id,
                    )

            query_coverage = (abs(query_end - query_start) + 1) / query_length * 100
            hits.append(
                BlastHit(
                    query_id=query_id,
                    subject_id=subject_id,
                    percent_identity=percent_identity,
                    alignment_length=alignment_length,
                    query_length=query_length,
                    subject_length=subject_length,
                    query_start=query_start,
                    query_end=query_end,
                    subject_start=subject_start,
                    subject_end=subject_end,
                    evalue=evalue,
                    bit_score=bit_score,
                    query_coverage=query_coverage,
                    subject_taxid=staxid,
                    subject_taxids=staxids,
                )
            )

    return hits


def run_nbdl_blast_search(
    query_fasta: Path,
    output_file: Path,
    min_identity: float,
    min_coverage: float,
    database: str,
    blastdb_path: Optional[Path],
    num_threads: int,
    skip_existing: bool,
) -> Tuple[bool, List[BlastHit]]:
    """Run BLAST using ``stitle`` so NBDL taxids can be parsed from headers."""
    runner = BlastRunner(blastdb_path)

    if output_file.exists() and skip_existing:
        logger.info("Using existing BLAST output: %s", output_file)
    else:
        success = runner.run_blastn(
            query_fasta=query_fasta,
            output_file=output_file,
            database=database,
            num_threads=num_threads,
            outfmt=SEARCH_OUTFMT_WITH_TITLE,
        )
        if not success:
            return False, []

    try:
        hits = _parse_blast_results_with_nbdl_titles(output_file)
        filtered_hits = runner.filter_blast_hits(
            hits,
            min_identity=min_identity,
            min_coverage=min_coverage,
        )
        return True, filtered_hits
    except Exception as exc:
        logger.error("Error processing NBDL-aware BLAST results: %s", exc)
        return False, []


def _process_taxonomy_and_representatives(args, filtered_hits: List[BlastHit]) -> Dict[str, Optional[Path]]:
    return process_blast_results_for_taxonomy(
        blast_hits=filtered_hits,
        output_dir=args.output_dir,
        rank=args.rank,
        database=args.database,
        blastdb_path=args.blastdb_path,
    )


def search_command(args) -> int:
    if not args.query_fasta.exists():
        logger.error("Query FASTA file not found: %s", args.query_fasta)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    sequences = FastaReader.read_fasta(args.query_fasta)
    search_output = args.output_dir / "blast_results.txt"

    success, filtered_hits = run_nbdl_blast_search(
        query_fasta=args.query_fasta,
        output_file=search_output,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        database=args.database,
        blastdb_path=args.blastdb_path,
        num_threads=args.num_threads,
        skip_existing=not args.overwrite_existing_blast,
    )
    if not success:
        return 1

    _process_taxonomy_and_representatives(args, filtered_hits)
    rewrite_blast_hits(filtered_hits, args.output_dir / "blast_results_annotated.txt", header=True)

    tree_files_map = {}
    if not args.skip_alignment:
        hits_by_query = defaultdict(list)
        for hit in filtered_hits:
            hits_by_query[hit.query_id].append(hit)

        tree_jobs = []
        query_job_info = {}
        for query_id, query_hits in hits_by_query.items():
            query_dir = args.output_dir / _safe_id(query_id)
            result = process_query_alignment_and_tree_parallel(
                query_id=query_id,
                query_dir=query_dir,
                blast_hits=query_hits,
                taxonomic_rank=args.tree_label_rank,
                query_fasta=args.query_fasta,
                num_threads=args.num_threads,
                background_tree=True,
            )
            if result["tree_job"]:
                tree_jobs.append(result["tree_job"])
                query_job_info[str(result["tree_file"])] = {
                    "query_id": query_id,
                    "tree_file": result["tree_file"],
                    "labeled_tree_path": result["labeled_tree_path"],
                    "blast_hits": result["blast_hits"],
                    "taxonomic_rank": result["taxonomic_rank"],
                }

        if tree_jobs:
            tree_results = IQTreeBuilder.wait_for_tree_jobs(tree_jobs)
            for tree_path, tree_success in tree_results.items():
                if not tree_success or tree_path not in query_job_info:
                    continue
                job_info = query_job_info[tree_path]
                tree_files_map[job_info["query_id"]] = job_info["tree_file"]
                if args.tree_label_rank == "no_rank":
                    logger.info("Skipping tree relabeling in no_rank mode for %s", job_info["query_id"])
                    continue
                IQTreeBuilder.relabel_tree_with_taxonomy(
                    tree_file=job_info["tree_file"],
                    blast_hits=job_info["blast_hits"],
                    output_tree=job_info["labeled_tree_path"],
                    taxonomic_rank=job_info["taxonomic_rank"],
                )

    output_classification = args.output_classification or args.output_dir / f"{args.query_fasta.stem}_classification.tsv"
    generate_classification_summary(
        sequences=sequences,
        blast_hits=filtered_hits,
        output_file=output_classification,
        rank=args.rank,
        group_rank=args.rank,
        tree_label_rank=args.tree_label_rank,
        tree_files=tree_files_map or None,
    )
    logger.info("NBDL-aware search workflow completed")
    return 0


def _group_hits_for_nbdl(filtered_hits: List[BlastHit], group_rank: str):
    if group_rank == "no_rank":
        grouped = defaultdict(lambda: defaultdict(list))
        for hit in filtered_hits:
            grouped["no_rank"][hit.query_id].append(hit)
        return dict(grouped)
    return group_hits_by_group_rank(filtered_hits, group_rank)


def grouped_command(args) -> int:
    if not args.query_fasta.exists():
        logger.error("Query FASTA file not found: %s", args.query_fasta)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    sequences = FastaReader.read_fasta(args.query_fasta)
    search_output = args.output_dir / "blast_results.txt"

    success, filtered_hits = run_nbdl_blast_search(
        query_fasta=args.query_fasta,
        output_file=search_output,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        database=args.database,
        blastdb_path=args.blastdb_path,
        num_threads=args.num_threads,
        skip_existing=not args.overwrite_existing_blast,
    )
    if not success:
        return 1

    _process_taxonomy_and_representatives(args, filtered_hits)
    rewrite_blast_hits(filtered_hits, args.output_dir / "blast_results_annotated.txt", header=True)

    consistency = check_alignment_consistency(filtered_hits, tolerance=args.alignment_tolerance)
    inconsistent = sum(1 for status in consistency.values() if not status)
    if inconsistent:
        logger.warning("Found %s reference sequences with inconsistent alignments", inconsistent)

    grouped_hits = _group_hits_for_nbdl(filtered_hits, args.group_rank)
    if not grouped_hits:
        logger.error("No groups found after grouping by %s", args.group_rank)
        return 1

    group_results = {}
    for group_id, query_hits_map in grouped_hits.items():
        group_name = group_id
        group_dir = args.output_dir / _safe_id(group_name)
        group_dir.mkdir(parents=True, exist_ok=True)
        combined_fasta = group_dir / f"{_safe_id(group_name)}_combined.fasta"

        success = create_grouped_fasta_with_queries(
            group_tid=group_id,
            group_name=group_name,
            query_hits_map=query_hits_map,
            labeling_rank=args.rank,
            query_fasta=args.query_fasta,
            output_fasta=combined_fasta,
            database=args.database,
            blastdb_path=args.blastdb_path,
        )
        if success:
            group_results[group_id] = {
                "name": group_name,
                "dir": group_dir,
                "fasta": combined_fasta,
                "query_ids": list(query_hits_map.keys()),
                "hits": [hit for hits in query_hits_map.values() for hit in hits],
            }

    tree_files_map = {}
    if not args.skip_alignment:
        tree_jobs = []
        group_job_info = {}
        for group_id, group_info in group_results.items():
            result = process_grouped_alignment_and_tree_parallel(
                group_name=group_info["name"],
                group_dir=group_info["dir"],
                taxonomic_rank=args.tree_label_rank,
                blast_hits=group_info["hits"],
                query_ids=group_info["query_ids"],
                num_threads=args.num_threads,
                background_tree=True,
            )
            if result["tree_job"]:
                tree_jobs.append(result["tree_job"])
                group_job_info[str(result["tree_file"])] = {
                    "group_id": group_id,
                    "group_name": group_info["name"],
                    "tree_file": result["tree_file"],
                    "labeled_tree_path": result["labeled_tree_path"],
                    "blast_hits": result["blast_hits"],
                    "taxonomic_rank": result["taxonomic_rank"],
                }

        if tree_jobs:
            tree_results = IQTreeBuilder.wait_for_tree_jobs(tree_jobs)
            for tree_path, tree_success in tree_results.items():
                if not tree_success or tree_path not in group_job_info:
                    continue
                job_info = group_job_info[tree_path]
                for query_id in group_results[job_info["group_id"]]["query_ids"]:
                    tree_files_map[query_id] = job_info["tree_file"]
                if args.tree_label_rank == "no_rank":
                    logger.info("Skipping tree relabeling in no_rank mode for %s", job_info["group_name"])
                    continue
                IQTreeBuilder.relabel_tree_with_taxonomy(
                    tree_file=job_info["tree_file"],
                    blast_hits=job_info["blast_hits"],
                    output_tree=job_info["labeled_tree_path"],
                    taxonomic_rank=job_info["taxonomic_rank"],
                )

    output_classification = args.output_classification or args.output_dir / f"{args.query_fasta.stem}_classification.tsv"
    generate_classification_summary(
        sequences=sequences,
        blast_hits=filtered_hits,
        output_file=output_classification,
        rank=args.rank,
        group_rank=args.group_rank,
        tree_label_rank=args.tree_label_rank,
        tree_files=tree_files_map or None,
    )
    logger.info("NBDL-aware grouped workflow completed")
    return 0


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("query_fasta", type=Path, help="Path to query FASTA file")
    parser.add_argument("output_dir", type=Path, help="Output directory for results")
    parser.add_argument("--rank", default="genus", choices=RANK_CHOICES, help="Representative rank, or no_rank for flat databases")
    parser.add_argument("--tree-label-rank", default="genus", choices=RANK_CHOICES, help="Tree labeling rank, or no_rank")
    parser.add_argument("--min-identity", type=float, default=90.0, help="Minimum percent identity")
    parser.add_argument("--min-coverage", type=float, default=80.0, help="Minimum query coverage percentage")
    parser.add_argument("--database", default="core_nt", help="BLAST database name")
    parser.add_argument("--blastdb-path", type=Path, default=None, help="Path to BLAST database directory")
    parser.add_argument("--num-threads", type=int, default=1, help="Threads for BLAST/alignment")
    parser.add_argument("--overwrite-existing-blast", action="store_true", help="Overwrite existing BLAST results")
    parser.add_argument("--skip-alignment", action="store_true", help="Skip alignment and tree building")
    parser.add_argument("--output-classification", type=Path, default=None, help="Path to output classification TSV")


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="eplace-nbdl",
        description="Run ePLACE with NBDL/CSIRO FASTA header taxid parsing",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging verbosity",
    )
    subparsers = parser.add_subparsers(dest="command")

    search_parser = subparsers.add_parser("search", help="Run NBDL-aware individual search workflow")
    _add_common_args(search_parser)

    grouped_parser = subparsers.add_parser("grouped", help="Run NBDL-aware grouped workflow")
    _add_common_args(grouped_parser)
    grouped_parser.add_argument("--group-rank", default="class", choices=RANK_CHOICES, help="Grouping rank, or no_rank")
    grouped_parser.add_argument("--alignment-tolerance", type=int, default=50, help="Maximum coordinate difference for consistency checks")

    args = parser.parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    if args.command == "search":
        return search_command(args)
    if args.command == "grouped":
        return grouped_command(args)

    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
