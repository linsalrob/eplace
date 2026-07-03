#!/usr/bin/env python
"""
ePLACE: environmental Phylogenetic Localisation and Clade Estimation

Main command-line interface for ePLACE toolkit.
Provides unified access to database download, BLAST analysis, and grouped workflows.
"""

import sys
import json
import argparse
import logging
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from collections import defaultdict

from .ncbi_download import (
    setup_ncbi_database,
    setup_mmseqs_database,
    setup_mmseqs_taxonomy,
    check_available_memory_gb
)
from .blast_analysis import run_blast_search, run_mmseqs_search, validate_mmseqs_memory_limit, FastaReader
from .taxonomy import (
    process_blast_results_for_taxonomy,
    rewrite_blast_hits,
    generate_classification_summary
)
from .alignment import (
    process_query_alignment_and_tree_parallel,
    process_grouped_alignment_and_tree_parallel,
    check_alignment_consistency,
    group_hits_by_group_rank,
    create_grouped_fasta_with_queries,
    IQTreeBuilder,
    concatenate_all_groups_and_build_tree
)

# Configure logging (level is overridden at runtime via --log-level)
RANK_CHOICES = ['phylum', 'class', 'order', 'family', 'genus', 'species', 'no_rank']
logging.basicConfig(
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

MMSEQS_DOWNLOAD_MIN_MEMORY_GB = 64.0
MMSEQS_TAXONOMY_MIN_MEMORY_GB = 128.0
MMSEQS_DB_PATH_FALLBACK_MSG = "$MMSEQS_DB_DIR, then $MMSEQS2DB, or ~/mmseqs2db"


def _looks_like_nt_database_name(name: str) -> bool:
    """Return True when a database name/path component denotes an NT database."""
    normalized = name.strip().lower()
    return (
        normalized == "nt"
        or normalized.startswith("nt.")
        or normalized.startswith("nt_")
    )


def _mmseqs_memory_limit_type(value: str) -> str:
    """Argparse type for ``--mmseqs-memory-limit``.

    Delegates to :func:`~eplace_lib.blast_analysis.validate_mmseqs_memory_limit`
    and converts ``ValueError`` into ``argparse.ArgumentTypeError`` so that the
    full error message is shown to the user.

    Args:
        value: Raw string value from the command line.

    Returns:
        The validated memory limit string unchanged.

    Raises:
        argparse.ArgumentTypeError: If ``value`` is not a valid MMseqs2
            memory limit string.
    """
    try:
        return validate_mmseqs_memory_limit(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def _get_cli_version() -> str:
    """Get package version from installed metadata."""
    try:
        return version("eplace")
    except PackageNotFoundError:
        return "unknown"


def _write_search_metadata(
    output_dir: Path,
    search_backend: str,
    database_name: str,
    database_path: str,
    database_source: str
) -> None:
    """Write backend and database provenance to a JSON metadata file.

    The file is written to ``<output_dir>/search_metadata.json`` and contains
    the following keys:

    * ``search_backend`` – the sequence search tool used (``blast`` or ``mmseqs2``)
    * ``database_name``  – the database name passed to the search tool
    * ``database_path``  – the directory that holds the database files
    * ``database_source`` – a human-readable label describing where the
      database originates (e.g. ``ncbi_core_nt`` or a user-supplied string
      such as ``ncbi_core_nt_derived_mmseqs2``)

    Args:
        output_dir: Output directory for the current run.
        search_backend: Name of the search backend (``blast`` or ``mmseqs2``).
        database_name: Name of the database used for the search.
        database_path: Path to the database directory.
        database_source: Provenance label for the database.
    """
    metadata = {
        "search_backend": search_backend,
        "database_name": database_name,
        "database_path": database_path,
        "database_source": database_source,
    }
    metadata_file = output_dir / "search_metadata.json"
    try:
        with open(metadata_file, "w") as f:
            json.dump(metadata, f, indent=2)
        logger.info(f"Search metadata written to: {metadata_file}")
    except OSError as e:
        logger.warning(f"Could not write search metadata file: {e}")


def _effective_mmseqs_database(args) -> str:
    """Return the effective MMseqs2 database name from parsed arguments.

    If ``--mmseqs-database`` was given explicitly, that value is used.
    Otherwise the value of ``--database`` (which defaults to ``core_nt``) is
    returned as a fallback so that the same reference corpus is targeted as
    in the BLAST workflow.

    Args:
        args: Parsed argument namespace from argparse.

    Returns:
        The database name to pass to MMseqs2Runner.
    """
    return args.mmseqs_database if args.mmseqs_database else args.database


def _log_mmseqs_database_warnings(args, mmseqs_database: str) -> None:
    """Emit advisory warnings when MMseqs2 database configuration is incomplete.

    Warns the user if any of the following are not set:
    - ``--mmseqs-database`` (defaulting to ``--database``)
    - ``--mmseqs-db-path`` (falling back to environment / home directory)
    - ``--mmseqs-db-source`` (no provenance label provided)

    Also warns when the configured database path or name looks like the full
    NCBI NT database, which requires hundreds of GB of RAM.

    Args:
        args: Parsed argument namespace from argparse.
        mmseqs_database: The effective MMseqs2 database name (after fallback).
    """
    if not args.mmseqs_database:
        logger.warning(
            "No explicit --mmseqs-database provided; defaulting to '%s'. "
            "Ensure that an MMseqs2-formatted database named '%s' exists at "
            "the configured database location, built from the same sequence "
            "content as BLAST core_nt.",
            mmseqs_database, mmseqs_database
        )
    if not args.mmseqs_db_path:
        logger.warning(
            "No --mmseqs-db-path provided. The MMseqs2 database directory "
            "will be resolved from %s. For "
            "reproducibility, explicitly specify --mmseqs-db-path pointing "
            "to a database built from the same sequence universe as BLAST "
            "core_nt.",
            MMSEQS_DB_PATH_FALLBACK_MSG
        )
    if not args.mmseqs_db_source:
        logger.warning(
            "No --mmseqs-db-source provided. Use this flag to document the "
            "provenance of your MMseqs2 database (e.g. "
            "'ncbi_core_nt_derived_mmseqs2') so that results remain "
            "reproducible and comparable with the BLAST workflow."
        )

    # Warn when the database looks like the full NCBI NT collection.
    nt_name_candidates = [mmseqs_database]
    if args.mmseqs_db_path:
        nt_name_candidates.append(Path(args.mmseqs_db_path).name)

    _looks_like_nt = any(
        _looks_like_nt_database_name(candidate)
        for candidate in nt_name_candidates
        if candidate
    )
    if _looks_like_nt:
        logger.warning(
            "You appear to be searching against the full NCBI NT database. "
            "This can require hundreds of GB of RAM. "
            "Consider using --mmseqs-memory-limit 400G or a smaller curated "
            "16S database (e.g. SILVA, GTDB rRNA, RDP, or NCBI 16S)."
        )


def _write_backend_search_metadata(args, mmseqs_database: str) -> None:
    """Resolve backend paths and write provenance metadata to the output directory.

    For BLAST the ``database_source`` is taken from ``--blast-db-source`` when
    provided; if not provided it defaults to the value of ``--database`` (e.g.
    ``core_nt``).  For MMseqs2 the ``database_source`` is the value of
    ``--mmseqs-db-source`` (empty string if not provided).

    Args:
        args: Parsed argument namespace from argparse.
        mmseqs_database: The effective MMseqs2 database name (after fallback).
    """
    if args.search_tool == 'mmseqs2':
        from .blast_analysis import MMseqs2Runner
        resolved_path = MMseqs2Runner(args.mmseqs_db_path).db_path
        _write_search_metadata(
            output_dir=args.output_dir,
            search_backend="mmseqs2",
            database_name=mmseqs_database,
            database_path=str(resolved_path),
            database_source=args.mmseqs_db_source or ""
        )
    else:
        from .blast_analysis import BlastRunner
        resolved_path = BlastRunner(args.blastdb_path).blastdb_path
        _write_search_metadata(
            output_dir=args.output_dir,
            search_backend="blast",
            database_name=args.database,
            database_path=str(resolved_path),
            database_source=args.blast_db_source or args.database
        )


def download_command(args):
    """Handle the download subcommand."""
    logger.info("=" * 60)
    logger.info("ePLACE Database Download")
    logger.info("=" * 60)

    if args.add_taxonomy and args.target == "blast":
        logger.error("--add-taxonomy is only valid when --target is mmseqs2 or both")
        return 1

    if args.target in ("blast", "both"):
        success, message = setup_ncbi_database(force_download=args.force)
        if not success:
            logger.error(f"✗ {message}")
            return 1
        logger.info(f"✓ {message}")

    mmseqs_db = None
    if args.target in ("mmseqs2", "both"):
        if not args.skip_memory_check:
            has_sufficient_memory, available_memory_gb = check_available_memory_gb(MMSEQS_DOWNLOAD_MIN_MEMORY_GB)
            if not has_sufficient_memory:
                logger.error(
                    "MMseqs2 NT database download requires at least %.0f GiB RAM "
                    "(detected %.1f GiB). Use --skip-memory-check to override.",
                    MMSEQS_DOWNLOAD_MIN_MEMORY_GB,
                    available_memory_gb
                )
                return 1

        success, message, mmseqs_db = setup_mmseqs_database(
            force_download=args.force,
            threads=args.mmseqs_threads,
            db_dir=args.mmseqs_db_dir
        )
        if not success:
            logger.error(f"✗ {message}")
            return 1
        logger.info(f"✓ {message}")

        if args.add_taxonomy:
            if not args.ncbi_taxonomy:
                logger.error("--ncbi-taxonomy is required when --add-taxonomy is used")
                return 1

            if mmseqs_db is None:
                logger.error("Unable to resolve MMseqs2 NT database path for taxonomy integration")
                return 1

            if not args.skip_memory_check:
                has_sufficient_memory, available_memory_gb = check_available_memory_gb(MMSEQS_TAXONOMY_MIN_MEMORY_GB)
                if not has_sufficient_memory:
                    logger.error(
                        "MMseqs2 taxonomy integration requires at least %.0f GiB RAM "
                        "(detected %.1f GiB). Use --skip-memory-check to override.",
                        MMSEQS_TAXONOMY_MIN_MEMORY_GB,
                        available_memory_gb
                    )
                    return 1

            success, message = setup_mmseqs_taxonomy(
                mmseqs_db=mmseqs_db,
                ncbi_taxonomy=args.ncbi_taxonomy,
                threads=args.mmseqs_threads,
                acc2taxid_dir=args.acc2taxid_dir,
                taxonomy_workdir=args.taxonomy_workdir,
                db_dir=args.mmseqs_db_dir
            )
            if not success:
                logger.error(f"✗ {message}")
                return 1
            logger.info(f"✓ {message}")

    return 0


def blast_command(args):
    """Handle the blast subcommand - individual workflow."""
    # Validate input file
    if not args.query_fasta.exists():
        logger.error(f"Query FASTA file not found: {args.query_fasta}")
        return 1
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set default output classification file if not provided
    if args.output_classification is None:
        base_name = args.query_fasta.stem
        for ext in ['.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn']:
            if base_name.endswith(ext):
                base_name = base_name[:-len(ext)]
                break
        args.output_classification = args.output_dir / f"{base_name}_classification.tsv"
    if not args.output_classification.is_absolute():
        args.output_classification = args.output_dir / args.output_classification
    
    skip_existing = not args.overwrite_existing_blast

    # Determine effective MMseqs2 database name early so it can be logged
    mmseqs_database = _effective_mmseqs_database(args)

    logger.info("=" * 60)
    logger.info("ePLACE Search Workflow")
    logger.info("=" * 60)
    logger.info(f"Query FASTA: {args.query_fasta}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Search tool: {args.search_tool}")
    logger.info(f"Taxonomic rank: {args.rank}")
    logger.info(f"Taxonomic rank for tree labeling: {args.tree_label_rank}")
    logger.info(f"Classification output file: {args.output_classification}")
    logger.info(f"Min identity: {args.min_identity}%")
    logger.info(f"Min coverage: {args.min_coverage}%")
    logger.info(f"Threads: {args.num_threads}")
    if args.search_tool == 'mmseqs2':
        logger.info(f"Search backend: mmseqs2")
        logger.info(f"MMseqs2 database name: {mmseqs_database}")
        logger.info(
            "MMseqs2 database path: %s",
            args.mmseqs_db_path or f"(default: {MMSEQS_DB_PATH_FALLBACK_MSG})"
        )
        logger.info(f"MMseqs2 database source: {args.mmseqs_db_source or '(not specified)'}")
        _log_mmseqs_database_warnings(args, mmseqs_database)
    else:
        logger.info(f"Search backend: blast")
        logger.info(f"Database: {args.database}")
        logger.info(f"BLAST database path: {args.blastdb_path or '(default: $BLASTDB or ~/blastdb)'}")
        logger.info(f"BLAST database source: {args.blast_db_source or args.database}")
    logger.info("=" * 60)

    # Step 1: Read query sequences
    logger.info("\n[Step 1/5] Reading query sequences...")
    try:
        sequences = FastaReader.read_fasta(args.query_fasta)
        logger.info(f"Found {len(sequences)} query sequences")
        for seq_id in sequences.keys():
            logger.info(f"  - {seq_id} ({len(sequences[seq_id])} bp)")
    except Exception:
        logger.exception(f"Error reading FASTA file")
        return 1
    
    # Step 2: Run sequence search
    # Determine output path and whether the search will actually execute so
    # that provenance metadata is written only when a fresh search runs.
    if args.search_tool == 'mmseqs2':
        logger.info("\n[Step 2/5] Running MMseqs2 search...")
        search_output = args.output_dir / "mmseqs_results.txt"
        search_ran = not (search_output.exists() and skip_existing)
        try:
            success, filtered_hits = run_mmseqs_search(
                query_fasta=args.query_fasta,
                output_file=search_output,
                min_identity=args.min_identity,
                min_coverage=args.min_coverage,
                database=mmseqs_database,
                db_path=args.mmseqs_db_path,
                num_threads=args.num_threads,
                sensitivity=args.mmseqs_sensitivity,
                skip_existing=skip_existing,
                search_type=args.mmseqs_search_type,
                memory_limit=args.mmseqs_memory_limit,
                timeout=args.timeout
            )

            if not success:
                logger.error("MMseqs2 search failed")
                return 1

            logger.info(f"MMseqs2 search completed successfully")
            logger.info(f"Found {len(filtered_hits)} hits after filtering")

        except Exception:
            logger.exception(f"Error during MMseqs2 search")
            return 1
    else:
        logger.info("\n[Step 2/5] Running BLAST search...")
        search_output = args.output_dir / "blast_results.txt"
        search_ran = not (search_output.exists() and skip_existing)
        try:
            success, filtered_hits = run_blast_search(
                query_fasta=args.query_fasta,
                output_file=search_output,
                min_identity=args.min_identity,
                min_coverage=args.min_coverage,
                database=args.database,
                blastdb_path=args.blastdb_path,
                num_threads=args.num_threads,
                skip_existing=skip_existing
            )

            if not success:
                logger.error("BLAST search failed")
                return 1

            logger.info(f"BLAST search completed successfully")
            logger.info(f"Found {len(filtered_hits)} hits after filtering")

        except Exception:
            logger.exception(f"Error during BLAST search")
            return 1

    # Write search metadata for provenance tracking only when a fresh search
    # was executed; skip when reusing existing results to avoid overwriting
    # accurate provenance from the original run with potentially different
    # CLI args.
    if search_ran:
        _write_backend_search_metadata(args, mmseqs_database)

    # Step 3: Group hits by query and display summary
    logger.info("\n[Step 3/5] Analyzing search results...")
    hits_by_query = defaultdict(int)
    for hit in filtered_hits:
        hits_by_query[hit.query_id] += 1
    
    for query_id, count in hits_by_query.items():
        logger.info(f"  {query_id}: {count} hits")
    
    # Step 4: Extract representative sequences
    logger.info(f"\n[Step 4/5] Extracting representative sequences (rank: {args.rank})...")
    
    try:
        results = process_blast_results_for_taxonomy(
            blast_hits=filtered_hits,
            output_dir=args.output_dir,
            rank=args.rank,
            database=args.database,
            blastdb_path=args.blastdb_path
        )
        
        logger.info("\nResults:")
        for query_id, output_fasta in results.items():
            if output_fasta:
                logger.info(f"  {query_id}: {output_fasta}")
            else:
                logger.warning(f"  {query_id}: Failed to extract sequences")
        
    except Exception:
        logger.exception(f"Error extracting sequences")
        return 1

    logger.info("Rewriting the blast output file with the new annotations")
    try:
        rewrite_blast_hits(
            blast_hits=filtered_hits,
            output_file=args.output_dir / "blast_results_annotated.txt",
            header=True
        )
    except Exception:
        logger.exception(f"Error rewriting the blast hits")
        return 1

    # Step 5: Align sequences and build trees (if not skipped)
    if not args.skip_alignment:
        logger.info("\n[Step 5/5] Aligning sequences and building phylogenetic trees...")
        
        # Group hits by query for processing
        hits_by_query_map = defaultdict(list)
        for hit in filtered_hits:
            hits_by_query_map[hit.query_id].append(hit)
        
        # Process all queries to do trimming and alignment
        # and start tree building in background
        tree_jobs = []
        query_job_info = {}
        
        for query_id, query_hits in hits_by_query_map.items():
            logger.info(f"\nProcessing alignment for query: {query_id}")
            
            # Get the query directory
            safe_query_id = query_id.replace('|', '_').replace('/', '_')
            query_dir = args.output_dir / safe_query_id
            
            if not query_dir.exists():
                logger.warning(f"Query directory not found: {query_dir}")
                continue
            
            # Process alignment and tree
            try:
                result = process_query_alignment_and_tree_parallel(
                    query_id=query_id,
                    query_dir=query_dir,
                    blast_hits=query_hits,
                    taxonomic_rank=args.rank,
                    query_fasta=args.query_fasta,
                    num_threads=args.num_threads,
                    background_tree=True
                )
                
                # Log progress
                if result['alignment']:
                    logger.info(f"  ✓ Alignment: {result['alignment']}")
                if result['trimmed_fasta']:
                    logger.info(f"  ✓ Trimmed sequences: {result['trimmed_fasta']}")
                
                # Collect tree jobs for later waiting
                if result['tree_job']:
                    tree_jobs.append(result['tree_job'])
                    query_job_info[str(result['tree_file'])] = {
                        'query_id': query_id,
                        'tree_file': result['tree_file'],
                        'labeled_tree_path': result['labeled_tree_path'],
                        'blast_hits': result['blast_hits'],
                        'taxonomic_rank': result['taxonomic_rank']
                    }
                    logger.info(f"  ✓ Tree building started in background")
                    
            except Exception as e:
                logger.error(f"Error processing {query_id}: {e}")
                continue
        
        # Wait for all tree building jobs to complete
        if tree_jobs:
            logger.info(f"\n{'='*60}")
            logger.info(f"Waiting for {len(tree_jobs)} tree building jobs to complete...")
            logger.info(f"This may take a while depending on the size of the alignments.")
            logger.info(f"{'='*60}\n")
            
            tree_results = IQTreeBuilder.wait_for_tree_jobs(tree_jobs)
            
            # Relabel trees that completed successfully
            logger.info("\nRelabeling trees with taxonomic names...")
            for tree_path, success in tree_results.items():
                if success and tree_path in query_job_info:
                    job_info = query_job_info[tree_path]
                    query_id = job_info['query_id']
                    tree_file = job_info['tree_file']
                    labeled_tree = job_info['labeled_tree_path']
                    blast_hits = job_info['blast_hits']
                    taxonomic_rank = job_info['taxonomic_rank']
                    
                    logger.info(f"  Processing {query_id}...")
                    if IQTreeBuilder.relabel_tree_with_taxonomy(
                        tree_file=tree_file,
                        blast_hits=blast_hits,
                        output_tree=labeled_tree,
                        taxonomic_rank=taxonomic_rank
                    ):
                        logger.info(f"    ✓ Labeled tree: {labeled_tree}")
                    else:
                        logger.warning(f"    ! Failed to relabel tree, but unlabeled tree is available: {tree_file}")
        
        logger.info(f"\nAlignment and tree building completed for {len(hits_by_query_map)} queries")
    else:
        logger.info("\n[Step 5/5] Skipping alignment and tree building (--skip-alignment)")

    # Step 6: Generate classification summary TSV
    logger.info("\nGenerating classification summary TSV file...")
    
    # Collect tree file paths for each query (if trees were built)
    tree_files_map = {}
    if not args.skip_alignment and 'query_job_info' in locals():
        for tree_path, job_info in query_job_info.items():
            query_id = job_info['query_id']
            tree_file = job_info['tree_file']
            # Use the tree file (unlabeled version) as it has the original subject IDs
            if tree_file and tree_file.exists():
                tree_files_map[query_id] = tree_file
    
    try:
        
        print("DEBUG classification output:", args.output_classification)
        print("DEBUG number sequences:", len(sequences))
        print("DEBUG number filtered_hits:", len(filtered_hits))
        print("DEBUG tree_files_map:", tree_files_map)
        
        success = generate_classification_summary(
            sequences=sequences,
            blast_hits=filtered_hits,
            output_file=args.output_classification,
            rank=args.rank,
            group_rank=args.rank,  # For individual workflow, group_rank same as rank
            tree_label_rank=args.tree_label_rank,
            tree_files=tree_files_map if tree_files_map else None
        )
        
        print("DEBUG classification success:", success)
        
        if success:
            logger.info(f"✓ Classification summary: {args.output_classification}")
        else:
            logger.warning("Failed to generate classification summary")
    except Exception as e:
        logger.error(f"Error generating classification summary: {e}")

    logger.info("\n" + "=" * 60)
    logger.info("Workflow completed successfully!")
    logger.info("=" * 60)
    
    return 0


def relabel_command(args):
    """Handle the relabel subcommand - relabel tree with taxonomy."""
    # Validate input files
    if not args.blast_output.exists():
        logger.error(f"BLAST output file not found: {args.blast_output}")
        return 1
    
    if not args.tree_file.exists():
        logger.error(f"Tree file not found: {args.tree_file}")
        return 1
    
    logger.info("=" * 60)
    logger.info("ePLACE Tree Relabeling")
    logger.info("=" * 60)
    logger.info(f"BLAST output: {args.blast_output}")
    logger.info(f"Tree file: {args.tree_file}")
    logger.info(f"Taxonomic rank: {args.rank}")
    logger.info(f"Output tree: {args.output_tree}")
    logger.info("=" * 60)
    
    # Step 1: Parse BLAST results
    logger.info("\n[Step 1/3] Parsing BLAST results...")
    try:
        from .blast_analysis import BlastRunner
        runner = BlastRunner(blastdb_path=args.blastdb_path)
        blast_hits = runner.parse_blast_results(args.blast_output)
        logger.info(f"Found {len(blast_hits)} BLAST hits")
    except Exception:
        logger.exception("Error parsing BLAST results")
        return 1
    
    # Step 2: Extract taxonomy information
    logger.info("\n[Step 2/3] Extracting taxonomy information...")
    try:
        from .taxonomy import TaxonomyExtractor
        tax_extractor = TaxonomyExtractor()
        
        # Get all unique subject tax IDs
        subject_taxids = {hit.subject_taxid for hit in blast_hits}
        logger.info(f"Found {len(subject_taxids)} unique taxonomy IDs")
        
        # Parse taxonomic information
        tax_dict = tax_extractor.parse_taxids(list(subject_taxids))
        
        # Add taxonomy to all hits
        for hit in blast_hits:
            hit.subject_taxonomy = tax_dict.get(hit.subject_taxid)
        
        logger.info(f"Taxonomy information added to all hits")
    except Exception:
        logger.exception("Error extracting taxonomy information")
        return 1
    
    # Step 3: Relabel tree
    logger.info("\n[Step 3/3] Relabeling tree...")
    try:
        # Read the tree file
        with open(args.tree_file, 'r') as f:
            tree_string = f.read()
        
        # Create mapping of sequence accession to taxonomic name
        label_map = {}
        for hit in blast_hits:
            accession = hit.get_accession()
            
            if not hit.subject_taxonomy:
                logger.warning(f"No taxonomy found for {hit.subject_id}")
                continue
            
            # Handle special case for species rank - use "genus species"
            if args.rank == 'species':
                genus_info = hit.subject_taxonomy.get('genus')
                species_info = hit.subject_taxonomy.get('species')
                
                if genus_info and species_info:
                    genus_name = genus_info[1]
                    species_name = species_info[1]
                    label = f"{genus_name} {species_name}"
                elif species_info:
                    label = species_info[1]
                elif genus_info:
                    label = genus_info[1]
                else:
                    logger.warning(f"No genus or species taxonomy for {hit.subject_id}")
                    continue
            else:
                # For other ranks, just use the rank name
                if args.rank not in hit.subject_taxonomy:
                    logger.warning(f"Rank {args.rank} not found for {hit.subject_id}")
                    continue
                label = hit.subject_taxonomy[args.rank][1]
            
            # Clean up the label for tree format (Newick format constraints)
            clean_label = (label.replace(' ', '_')
                          .replace(':', '_')
                          .replace('(', '_')
                          .replace(')', '_')
                          .replace(',', '_')
                          .replace(';', '_'))
            
            label_map[accession] = clean_label
        
        logger.info(f"Created label mapping for {len(label_map)} sequences")
        
        # Replace sequence IDs with taxonomic names
        for seq_id, tax_name in label_map.items():
            # Handle normal sequences (not reversed)
            tree_string = tree_string.replace(f"({seq_id}:", f"({tax_name}:")
            tree_string = tree_string.replace(f",{seq_id}:", f",{tax_name}:")
            tree_string = tree_string.replace(f" {seq_id}:", f" {tax_name}:")
            
            # Handle sequences with _R_ prefix (reversed by MAFFT)
            reversed_seq_id = f"_R_{seq_id}"
            reversed_label = f"{tax_name}_R"
            tree_string = tree_string.replace(f"({reversed_seq_id}:", f"({reversed_label}:")
            tree_string = tree_string.replace(f",{reversed_seq_id}:", f",{reversed_label}:")
            tree_string = tree_string.replace(f" {reversed_seq_id}:", f" {reversed_label}:")
        
        # Write the relabeled tree
        with open(args.output_tree, 'w') as f:
            f.write(tree_string)
        
        logger.info(f"✓ Relabeled tree saved to: {args.output_tree}")
        
    except Exception:
        logger.exception("Error relabeling tree")
        return 1
    
    logger.info("\n" + "=" * 60)
    logger.info("Tree relabeling completed successfully!")
    logger.info("=" * 60)
    
    return 0


def grouped_command(args):
    """Handle the grouped subcommand - grouped workflow."""
    # Validate input file
    if not args.query_fasta.exists():
        logger.error(f"Query FASTA file not found: {args.query_fasta}")
        return 1
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set default output classification file if not provided
    if args.output_classification is None:
        base_name = args.query_fasta.stem
        for ext in ['.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn']:
            if base_name.endswith(ext):
                base_name = base_name[:-len(ext)]
                break
        args.output_classification = args.output_dir / f"{base_name}_classification.tsv"
    if not args.output_classification.is_absolute():
        args.output_classification = args.output_dir / args.output_classification
    
    skip_existing = not args.overwrite_existing_blast

    # Determine effective MMseqs2 database name early so it can be logged
    mmseqs_database = _effective_mmseqs_database(args)

    logger.info("=" * 60)
    logger.info("ePLACE Grouped Search Workflow")
    logger.info("=" * 60)
    logger.info(f"Query FASTA: {args.query_fasta}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Search tool: {args.search_tool}")
    logger.info(f"Representative rank: {args.rank}")
    logger.info(f"Grouping rank: {args.group_rank}")
    logger.info(f"Tree labeling rank: {args.tree_label_rank}")
    logger.info(f"Classification output file: {args.output_classification}")
    logger.info(f"Min identity: {args.min_identity}%")
    logger.info(f"Min coverage: {args.min_coverage}%")
    logger.info(f"Threads: {args.num_threads}")
    if args.search_tool == 'mmseqs2':
        logger.info(f"Search backend: mmseqs2")
        logger.info(f"MMseqs2 database name: {mmseqs_database}")
        logger.info(
            "MMseqs2 database path: %s",
            args.mmseqs_db_path or f"(default: {MMSEQS_DB_PATH_FALLBACK_MSG})"
        )
        logger.info(f"MMseqs2 database source: {args.mmseqs_db_source or '(not specified)'}")
        _log_mmseqs_database_warnings(args, mmseqs_database)
    else:
        logger.info(f"Search backend: blast")
        logger.info(f"Database: {args.database}")
        logger.info(f"BLAST database path: {args.blastdb_path or '(default: $BLASTDB or ~/blastdb)'}")
        logger.info(f"BLAST database source: {args.blast_db_source or args.database}")
    logger.info("=" * 60)

    # Step 1: Read query sequences
    logger.info("\n[Step 1/9] Reading query sequences...")
    try:
        sequences = FastaReader.read_fasta(args.query_fasta)
        logger.info(f"Found {len(sequences)} query sequences")
        for seq_id in sequences.keys():
            logger.info(f"  - {seq_id} ({len(sequences[seq_id])} bp)")
    except Exception:
        logger.exception(f"Error reading FASTA file")
        return 1
    
    # Step 2: Run sequence search
    # Determine output path and whether the search will actually execute so
    # that provenance metadata is written only when a fresh search runs.
    if args.search_tool == 'mmseqs2':
        logger.info("\n[Step 2/9] Running MMseqs2 search...")
        search_output = args.output_dir / "mmseqs_results.txt"
        search_ran = not (search_output.exists() and skip_existing)
        try:
            success, filtered_hits = run_mmseqs_search(
                query_fasta=args.query_fasta,
                output_file=search_output,
                min_identity=args.min_identity,
                min_coverage=args.min_coverage,
                database=mmseqs_database,
                db_path=args.mmseqs_db_path,
                num_threads=args.num_threads,
                sensitivity=args.mmseqs_sensitivity,
                skip_existing=skip_existing,
                search_type=args.mmseqs_search_type,
                memory_limit=args.mmseqs_memory_limit,
                timeout=args.timeout
            )

            if not success:
                logger.error("MMseqs2 search failed")
                return 1

            logger.info(f"MMseqs2 search completed successfully")
            logger.info(f"Found {len(filtered_hits)} hits after filtering")

        except Exception:
            logger.exception(f"Error during MMseqs2 search")
            return 1
    else:
        logger.info("\n[Step 2/9] Running BLAST search...")
        search_output = args.output_dir / "blast_results.txt"
        search_ran = not (search_output.exists() and skip_existing)
        try:
            success, filtered_hits = run_blast_search(
                query_fasta=args.query_fasta,
                output_file=search_output,
                min_identity=args.min_identity,
                min_coverage=args.min_coverage,
                database=args.database,
                blastdb_path=args.blastdb_path,
                num_threads=args.num_threads,
                skip_existing=skip_existing
            )

            if not success:
                logger.error("BLAST search failed")
                return 1

            logger.info(f"BLAST search completed successfully")
            logger.info(f"Found {len(filtered_hits)} hits after filtering")

        except Exception:
            logger.exception(f"Error during BLAST search")
            return 1

    # Write search metadata for provenance tracking only when a fresh search
    # was executed; skip when reusing existing results to avoid overwriting
    # accurate provenance from the original run with potentially different
    # CLI args.
    if search_ran:
        _write_backend_search_metadata(args, mmseqs_database)

    # Step 3: Process taxonomy information
    logger.info(f"\n[Step 3/9] Processing taxonomy information (rank: {args.rank})...")
    
    try:
        results = process_blast_results_for_taxonomy(
            blast_hits=filtered_hits,
            output_dir=args.output_dir,
            rank=args.rank,
            database=args.database,
            blastdb_path=args.blastdb_path
        )
        
        logger.info("\nExtracted representative sequences per query:")
        for query_id, output_fasta in results.items():
            if output_fasta:
                logger.info(f"  {query_id}: {output_fasta}")
            else:
                logger.warning(f"  {query_id}: Failed to extract sequences")
        
    except Exception:
        logger.exception(f"Error extracting sequences")
        return 1

    logger.info("Rewriting the BLAST output file with the new annotations")
    try:
        rewrite_blast_hits(
            blast_hits=filtered_hits,
            output_file=args.output_dir / "blast_results_annotated.txt",
            header=True
        )
    except Exception as e:
        logger.exception(f"Error rewriting the blast hits: {e}")
        return 1

    # Step 4: Check alignment consistency
    logger.info("\n[Step 4/9] Checking alignment consistency...")
    consistency = check_alignment_consistency(
        blast_hits=filtered_hits,
        tolerance=args.alignment_tolerance
    )
    
    inconsistent_count = sum(1 for is_consistent in consistency.values() if not is_consistent)
    if inconsistent_count > 0:
        logger.warning(
            f"Found {inconsistent_count} reference sequences with inconsistent alignments "
            f"(tolerance: {args.alignment_tolerance} bp)"
        )
    else:
        logger.info("All reference sequences have consistent alignments across queries")

    # Step 5: Group hits by group_rank
    logger.info(f"\n[Step 5/9] Grouping hits by {args.group_rank}...")
    grouped_hits = group_hits_by_group_rank(filtered_hits, args.group_rank)
    
    if not grouped_hits:
        logger.error("No groups found after grouping by rank")
        return 1

    # Step 6: Create grouped FASTA files
    logger.info(f"\n[Step 6/9] Creating grouped FASTA files...")
    
    group_results = {}
    for group_tid, query_hits_map in grouped_hits.items():
        group_name = group_tid
        
        if not group_name:
            logger.warning(f"No group name found for group {group_tid}, skipping")
            continue
        
        logger.info(f"\nProcessing group: {group_name} ({group_tid})")
        logger.info(f"  Queries in this group: {', '.join(query_hits_map.keys())}")
        
        # Create group directory
        safe_group_name = group_name.replace(' ', '_').replace('/', '_').replace('|', '_')
        group_dir = args.output_dir / safe_group_name
        group_dir.mkdir(parents=True, exist_ok=True)
        
        # Create combined FASTA file
        combined_fasta = group_dir / f"{safe_group_name}_combined.fasta"
        success = create_grouped_fasta_with_queries(
            group_tid=group_tid,
            group_name=group_name,
            query_hits_map=query_hits_map,
            labeling_rank=args.rank,
            query_fasta=args.query_fasta,
            output_fasta=combined_fasta,
            database=args.database,
            blastdb_path=args.blastdb_path
        )
        
        if success:
            logger.info(f"  ✓ Created: {combined_fasta}")
            group_results[group_tid] = {
                'name': group_name,
                'dir': group_dir,
                'fasta': combined_fasta,
                'query_ids': list(query_hits_map.keys()),
                'hits': [hit for hits in query_hits_map.values() for hit in hits]
            }
        else:
            logger.error(f"  ✗ Failed to create grouped FASTA for {group_name}")

    # Step 7: Build alignments and trees for each group
    if not args.skip_alignment:
        logger.info("\n[Step 7/9] Building alignments and phylogenetic trees for each group...")
        
        # Process all groups to do trimming and alignment
        # and start tree building in background
        tree_jobs = []
        group_job_info = {}
        
        for group_tid, group_info in group_results.items():
            group_name = group_info['name']
            group_dir = group_info['dir']
            query_ids = group_info['query_ids']
            blast_hits = group_info['hits']
            
            logger.info(f"\nProcessing group: {group_name} ({group_tid})")
            
            try:
                result = process_grouped_alignment_and_tree_parallel(
                    group_name=group_name,
                    group_dir=group_dir,
                    taxonomic_rank=args.tree_label_rank,
                    blast_hits=blast_hits,
                    query_ids=query_ids,
                    num_threads=args.num_threads,
                    background_tree=True
                )
                
                # Log progress
                if result['alignment']:
                    logger.info(f"  ✓ Alignment: {result['alignment']}")
                if result['trimmed_fasta']:
                    logger.info(f"  ✓ Trimmed sequences: {result['trimmed_fasta']}")
                
                # Collect tree jobs for later waiting
                if result['tree_job']:
                    tree_jobs.append(result['tree_job'])
                    group_job_info[str(result['tree_file'])] = {
                        'group_tid': group_tid,
                        'group_name': group_name,
                        'tree_file': result['tree_file'],
                        'labeled_tree_path': result['labeled_tree_path'],
                        'blast_hits': result['blast_hits'],
                        'taxonomic_rank': result['taxonomic_rank']
                    }
                    logger.info(f"  ✓ Tree building started in background")
                    
            except Exception as e:
                logger.error(f"Error processing {group_name}: {e}")
                continue
        
        # Wait for all tree building jobs to complete
        if tree_jobs:
            logger.info(f"\n{'='*60}")
            logger.info(f"Waiting for {len(tree_jobs)} tree building jobs to complete...")
            logger.info(f"This may take a while depending on the size of the alignments.")
            logger.info(f"{'='*60}\n")
            
            tree_results = IQTreeBuilder.wait_for_tree_jobs(tree_jobs)
            
            # Relabel trees that completed successfully
            logger.info("\nRelabeling trees with taxonomic names...")
            for tree_path, success in tree_results.items():
                if success and tree_path in group_job_info:
                    job_info = group_job_info[tree_path]
                    group_name = job_info['group_name']
                    tree_file = job_info['tree_file']
                    labeled_tree = job_info['labeled_tree_path']
                    blast_hits = job_info['blast_hits']
                    taxonomic_rank = job_info['taxonomic_rank']
                    
                    logger.info(f"  Processing {group_name}...")
                    if IQTreeBuilder.relabel_tree_with_taxonomy(
                        tree_file=tree_file,
                        blast_hits=blast_hits,
                        output_tree=labeled_tree,
                        taxonomic_rank=taxonomic_rank
                    ):
                        logger.info(f"    ✓ Labeled tree: {labeled_tree}")
                    else:
                        logger.warning(f"    ! Failed to relabel tree, but unlabeled tree is available: {tree_file}")
        
        logger.info(f"\nAlignment and tree building completed for {len(group_results)} groups")
    else:
        logger.info("\n[Step 7/9] Skipping alignment and tree building (--skip-alignment)")

    # Step 8: Generate classification summary TSV
    logger.info("\n[Step 8/9] Generating classification summary TSV file...")
    
    # Collect tree file paths for each query (if trees were built)
    # In grouped workflow, multiple queries may share the same group tree
    tree_files_map = {}
    if not args.skip_alignment and 'group_job_info' in locals() and 'group_results' in locals():
        for tree_path, job_info in group_job_info.items():
            group_tid = job_info.get('group_tid')
            tree_file = job_info['tree_file']
            
            # Find which queries belong to this group
            if group_tid in group_results:
                query_ids = group_results[group_tid]['query_ids']
                # Map each query to this tree file
                if tree_file and tree_file.exists():
                    for query_id in query_ids:
                        tree_files_map[query_id] = tree_file
    
    try:
        success = generate_classification_summary(
            sequences=sequences,
            blast_hits=filtered_hits,
            output_file=args.output_classification,
            rank=args.rank,
            group_rank=args.group_rank,
            tree_label_rank=args.tree_label_rank,
            tree_files=tree_files_map if tree_files_map else None
        )
        
        if success:
            logger.info(f"✓ Classification summary: {args.output_classification}")
        else:
            logger.warning("Failed to generate classification summary")
    except Exception as e:
        logger.error(f"Error generating classification summary: {e}")

    # Step 9: Build combined tree from all groups (optional)
    if not args.skip_alignment and group_results and args.combined_tree_label_rank is not None:
        logger.info("\n[Step 9/9] Building combined tree from all groups...")
        try:
            combined_results = concatenate_all_groups_and_build_tree(
                output_dir=args.output_dir,
                query_fasta=args.query_fasta,
                classification_file=args.output_classification,
                blast_hits=filtered_hits,
                combined_tree_label_rank=args.combined_tree_label_rank,
                num_threads=args.num_threads
            )
            
            if combined_results['alignment']:
                logger.info(f"✓ Combined alignment: {combined_results['alignment']}")
            if combined_results['tree']:
                logger.info(f"✓ Combined tree: {combined_results['tree']}")
            if combined_results['labeled_tree']:
                logger.info(f"✓ Combined labeled tree: {combined_results['labeled_tree']}")
        except Exception as e:
            logger.error(f"Error building combined tree: {e}")
    elif not args.skip_alignment and group_results and args.combined_tree_label_rank is None:
        logger.info("\n[Step 9/9] Skipping combined tree building (not requested)")
        logger.info("To build a combined tree from all groups, use --combined-tree-label-rank option")

    logger.info("\n" + "=" * 60)
    logger.info("Grouped workflow completed successfully!")
    logger.info(f"Processed {len(group_results)} taxonomic groups")
    logger.info("=" * 60)
    
    return 0


_LOG_LEVEL_CHOICES = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']


def _add_log_level_argument(p, *, is_top_level=False):
    """Add a --log-level argument to *p*.

    For the top-level parser the default is 'INFO'.  For subparsers the
    default is suppressed so that an explicit value on the top-level parser
    is not overwritten when the subparser is invoked without the flag.
    """
    p.add_argument(
        '--log-level',
        default='INFO' if is_top_level else argparse.SUPPRESS,
        choices=_LOG_LEVEL_CHOICES,
        help='Set logging verbosity level (default: INFO)'
    )


def main():
    """Main entry point for the ePLACE CLI."""
    parser = argparse.ArgumentParser(
        prog='eplace',
        description='ePLACE: environmental Phylogenetic Localisation and Clade Estimation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download NCBI database
  eplace download
  
  # Run individual search workflow
  eplace search query.fasta output_dir
  
  # Run grouped BLAST workflow
  eplace grouped query.fasta output_dir --group-rank order
  
  # Relabel a tree with taxonomic names (use blast_results.txt for BLAST, mmseqs_results.txt for MMseqs2)
  eplace relabel blast_results.txt input.treefile output.treefile --rank genus
  
For detailed help on each subcommand:
  eplace download --help
  eplace search --help
  eplace grouped --help
  eplace relabel --help

Documentation: https://github.com/linsalrob/eplace
        """
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {_get_cli_version()}'
    )

    _add_log_level_argument(parser, is_top_level=True)
    
    subparsers = parser.add_subparsers(
        dest='command',
        title='Available commands',
        description='Use these commands to access different ePLACE functionalities',
        help='Command to run'
    )
    
    # Download subcommand
    download_parser = subparsers.add_parser(
        'download',
        help='Download BLAST and/or MMseqs2 databases',
        description='Download and setup NCBI core_nt BLAST and MMseqs2 NT databases',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download BLAST core_nt database to default location ($BLASTDB or ~/blastdb)
  eplace download
  
  # Download MMseqs2 NT database to $MMSEQS_DB_DIR (or fallback path)
  eplace download --target mmseqs2 --mmseqs-db-dir /path/to/mmseqs_db

  # Download MMseqs2 NT and add taxonomy sidecar files
  eplace download --target mmseqs2 --add-taxonomy --ncbi-taxonomy /path/to/ncbi/taxonomy/current

  # Force redownload even if database exists (applies to selected target)
  eplace download --force

Notes:
  - BLAST DB location uses $BLASTDB or ~/blastdb
  - MMseqs2 DB location uses $MMSEQS_DB_DIR, then $MMSEQS2DB, or ~/mmseqs2db
  - MMseqs2 NT download typically requires at least 64 GiB RAM
  - MMseqs2 taxonomy integration typically requires at least 128 GiB RAM
        """
    )
    download_parser.add_argument(
        '--target',
        type=str,
        default='blast',
        choices=['blast', 'mmseqs2', 'both'],
        help='Database backend to download (default: blast)'
    )
    download_parser.add_argument(
        '--force',
        action='store_true',
        help='Force redownload even if database exists'
    )
    download_parser.add_argument(
        '--mmseqs-db-dir',
        type=Path,
        default=None,
        help='Path to MMseqs2 database root directory for NT downloads'
    )
    download_parser.add_argument(
        '--mmseqs-threads',
        type=int,
        default=1,
        help='Number of threads for MMseqs2 download/taxonomy commands (default: 1)'
    )
    download_parser.add_argument(
        '--add-taxonomy',
        action='store_true',
        help='Add NCBI taxonomy sidecar files to MMseqs2 NT after download'
    )
    download_parser.add_argument(
        '--ncbi-taxonomy',
        type=Path,
        default=None,
        help='Path to NCBI taxonomy dump directory containing nodes.dmp, names.dmp, merged.dmp'
    )
    download_parser.add_argument(
        '--acc2taxid-dir',
        type=Path,
        default=None,
        help='Path to accession2taxid files (defaults to $ACC2TAXID_DIR or <ncbi-taxonomy>/accession2taxid)'
    )
    download_parser.add_argument(
        '--taxonomy-workdir',
        type=Path,
        default=None,
        help='Working directory for MMseqs2 taxonomy mapping files'
    )
    download_parser.add_argument(
        '--skip-memory-check',
        action='store_true',
        help='Skip RAM preflight checks for MMseqs2 download and taxonomy integration'
    )
    _add_log_level_argument(download_parser)
    
    # Search subcommand (individual workflow)
    search_parser = subparsers.add_parser(
        'search',
        aliases=['blast'],
        help='Run sequence search with individual taxonomy analysis',
        description='Run sequence search and extract representative sequences per query',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with default parameters
  eplace search query.fasta output_dir
  
  # Specify taxonomic rank and custom thresholds
  eplace search query.fasta output_dir --rank genus --min-identity 95
  
  # Search + taxonomy/extraction only (skip alignment and tree building)
  eplace search query.fasta output_dir --skip-alignment
  
  # Use custom BLAST database location
  eplace search query.fasta output_dir --blastdb-path /path/to/blastdb

Notes:
  - Creates one phylogenetic tree per query sequence
  - Default filtering: 90% identity, 80% query coverage
  - Requires BLAST+ or MMseqs2 tools and optionally MAFFT and IQTree
        """
    )
    search_parser.add_argument(
        'query_fasta',
        type=Path,
        help='Path to query FASTA file'
    )
    search_parser.add_argument(
        'output_dir',
        type=Path,
        help='Output directory for results'
    )
    search_parser.add_argument(
        '--rank',
        type=str,
        default='genus',
        choices=RANK_CHOICES,
        help='Taxonomic rank for representative selection (default: genus)'
    )
    search_parser.add_argument(
        '--tree-label-rank',
        type=str,
        default='genus',
        choices=RANK_CHOICES,
        help='Taxonomic rank for tree labeling (default: genus)'
    )
    search_parser.add_argument(
        '--min-identity',
        type=float,
        default=90.0,
        help='Minimum percent identity for BLAST hits (default: 90.0)'
    )
    search_parser.add_argument(
        '--min-coverage',
        type=float,
        default=80.0,
        help='Minimum query coverage percentage (default: 80.0)'
    )
    search_parser.add_argument(
        '--database',
        type=str,
        default='core_nt',
        help='BLAST database name (default: core_nt)'
    )
    search_parser.add_argument(
        '--blastdb-path',
        type=Path,
        default=None,
        help='Path to BLAST database directory'
    )
    search_parser.add_argument(
        '--num-threads',
        type=int,
        default=1,
        help='Number of threads for search and alignment (default: 1)'
    )
    search_parser.add_argument(
        '--timeout',
        type=int,
        default=3600,
        help='Maximum runtime for MMseqs2 search in seconds (default: 3600). '
             'Only used when --search-tool mmseqs2 is specified.'
    )
    search_parser.add_argument(
        '--search-tool',
        type=str,
        default='blast',
        choices=['blast', 'mmseqs2'],
        help='Sequence search tool to use (default: blast)'
    )
    search_parser.add_argument(
        '--blast-db-source',
        type=str,
        default=None,
        help='Provenance label for the BLAST database, recorded in '
             'search_metadata.json for reproducibility '
             "(e.g. 'ncbi_core_nt' for the default NCBI core_nt, or a "
             "custom label when using a non-standard BLAST database). "
             'If not provided, defaults to the value of --database. '
             'Use this flag when --database points to a database other than '
             'the standard NCBI core_nt to document the database origin.'
    )
    search_parser.add_argument(
        '--mmseqs-database',
        type=str,
        default=None,
        help='MMseqs2 database name (default: same as --database). '
             'Only used when --search-tool mmseqs2 is specified. '
             'The recommended database is one built from the same sequence '
             'collection as BLAST core_nt (there is no official pre-built '
             'MMseqs2 core_nt database; users must create their own).'
    )
    search_parser.add_argument(
        '--mmseqs-db-path',
        type=Path,
        default=None,
        help='Path to the MMseqs2 database directory. '
             'Only used when --search-tool mmseqs2 is specified. '
             f'If not provided, falls back to {MMSEQS_DB_PATH_FALLBACK_MSG}. '
             'For results comparable with BLAST, this directory should contain '
             'a database built from the same sequence universe as core_nt.'
    )
    search_parser.add_argument(
        '--mmseqs-db-source',
        type=str,
        default=None,
        help='Provenance label for the MMseqs2 database, recorded in '
             'search_metadata.json for reproducibility '
             "(e.g. 'ncbi_core_nt_derived_mmseqs2'). "
             'Only used when --search-tool mmseqs2 is specified. '
             'This label documents where the database originates so that '
             'downstream interpretation can account for any differences '
             'from the BLAST core_nt reference corpus.'
    )
    search_parser.add_argument(
        '--mmseqs-sensitivity',
        type=float,
        default=5.7,
        help='MMseqs2 sensitivity setting, 1–7.5 (default: 5.7). '
             'Only used when --search-tool mmseqs2 is specified.'
    )
    search_parser.add_argument(
        '--mmseqs-search-type',
        type=int,
        default=3,
        help='MMseqs2 search type passed as --search-type to easy-search. '
             'Commonly used values: 2 (translated), 3 (nucleotide), '
             '4 (translated nucleotide backtrace). Default is 3 (nucleotide). '
             'See MMseqs2 documentation for all valid values. '
             'Only used when --search-tool mmseqs2 is specified.'
    )
    search_parser.add_argument(
        '--mmseqs-memory-limit',
        type=_mmseqs_memory_limit_type,
        default='400G',
        metavar='LIMIT',
        help='Maximum RAM for the MMseqs2 prefilter step, passed as '
             '--split-memory-limit to mmseqs easy-search (default: 400G). '
             'On smaller hosts, lower this explicitly (e.g. 16G or 32G). '
             'Accepts MMseqs2-style memory strings such as 64G, 400G, or 1T. '
             'Only used when --search-tool mmseqs2 is specified.'
    )
    search_parser.add_argument(
        '--overwrite-existing-blast',
        action='store_true',
        help='Overwrite existing search results'
    )
    search_parser.add_argument(
        '--skip-alignment',
        action='store_true',
        help='Skip alignment and tree building steps'
    )
    search_parser.add_argument(
        '--output-classification',
        type=Path,
        default=None,
        help='Path to output classification TSV file'
    )
    _add_log_level_argument(search_parser)
    
    # Grouped subcommand
    grouped_parser = subparsers.add_parser(
        'grouped',
        help='Run BLAST search with grouped taxonomy analysis',
        description='Run BLAST search and group sequences by taxonomic rank',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (group by class, default)
  eplace grouped query.fasta output_dir
  
  # Group by different taxonomic rank
  eplace grouped query.fasta output_dir --group-rank order
  
  # Specify both representative and grouping ranks
  eplace grouped query.fasta output_dir --rank genus --group-rank family
  
  # Skip alignment and tree building
  eplace grouped query.fasta output_dir --skip-alignment

Notes:
  - Creates one phylogenetic tree per taxonomic group
  - Groups queries that match the same taxonomic rank together
  - Useful for analyzing multiple related queries in one phylogenetic context
  - Default filtering: 90% identity, 80% query coverage
        """
    )
    grouped_parser.add_argument(
        'query_fasta',
        type=Path,
        help='Path to query FASTA file'
    )
    grouped_parser.add_argument(
        'output_dir',
        type=Path,
        help='Output directory for results'
    )
    grouped_parser.add_argument(
        '--rank',
        type=str,
        default='genus',
        choices=RANK_CHOICES,
        help='Taxonomic rank for representative selection (default: genus)'
    )
    grouped_parser.add_argument(
        '--group-rank',
        type=str,
        default='class',
        choices=RANK_CHOICES,
        help='Taxonomic rank for grouping sequences (default: class)'
    )
    grouped_parser.add_argument(
        '--tree-label-rank',
        type=str,
        default='genus',
        choices=RANK_CHOICES,
        help='Taxonomic rank for tree labeling (default: genus)'
    )
    grouped_parser.add_argument(
        '--combined-tree-label-rank',
        type=str,
        default=None,
        choices=RANK_CHOICES,
        help='Taxonomic rank for tree labeling for the combined tree. If not provided, the combined tree will not be built (to save time with large datasets).'
    )
    grouped_parser.add_argument(
        '--min-identity',
        type=float,
        default=90.0,
        help='Minimum percent identity for search hits (default: 90.0)'
    )
    grouped_parser.add_argument(
        '--min-coverage',
        type=float,
        default=80.0,
        help='Minimum query coverage percentage (default: 80.0)'
    )
    grouped_parser.add_argument(
        '--database',
        type=str,
        default='core_nt',
        help='BLAST database name (default: core_nt)'
    )
    grouped_parser.add_argument(
        '--blastdb-path',
        type=Path,
        default=None,
        help='Path to BLAST database directory'
    )
    grouped_parser.add_argument(
        '--num-threads',
        type=int,
        default=1,
        help='Number of threads for search and alignment (default: 1)'
    )
    grouped_parser.add_argument(
        '--timeout',
        type=int,
        default=3600,
        help='Maximum runtime for MMseqs2 search in seconds (default: 3600). '
             'Only used when --search-tool mmseqs2 is specified.'
    )
    grouped_parser.add_argument(
        '--search-tool',
        type=str,
        default='blast',
        choices=['blast', 'mmseqs2'],
        help='Sequence search tool to use (default: blast)'
    )
    grouped_parser.add_argument(
        '--blast-db-source',
        type=str,
        default=None,
        help='Provenance label for the BLAST database, recorded in '
             'search_metadata.json for reproducibility '
             "(e.g. 'ncbi_core_nt' for the default NCBI core_nt, or a "
             "custom label when using a non-standard BLAST database). "
             'If not provided, defaults to the value of --database. '
             'Use this flag when --database points to a database other than '
             'the standard NCBI core_nt to document the database origin.'
    )
    grouped_parser.add_argument(
        '--mmseqs-database',
        type=str,
        default=None,
        help='MMseqs2 database name (default: same as --database). '
             'Only used when --search-tool mmseqs2 is specified. '
             'The recommended database is one built from the same sequence '
             'collection as BLAST core_nt (there is no official pre-built '
             'MMseqs2 core_nt database; users must create their own).'
    )
    grouped_parser.add_argument(
        '--mmseqs-db-path',
        type=Path,
        default=None,
        help='Path to the MMseqs2 database directory. '
             'Only used when --search-tool mmseqs2 is specified. '
             f'If not provided, falls back to {MMSEQS_DB_PATH_FALLBACK_MSG}. '
             'For results comparable with BLAST, this directory should contain '
             'a database built from the same sequence universe as core_nt.'
    )
    grouped_parser.add_argument(
        '--mmseqs-db-source',
        type=str,
        default=None,
        help='Provenance label for the MMseqs2 database, recorded in '
             'search_metadata.json for reproducibility '
             "(e.g. 'ncbi_core_nt_derived_mmseqs2'). "
             'Only used when --search-tool mmseqs2 is specified. '
             'This label documents where the database originates so that '
             'downstream interpretation can account for any differences '
             'from the BLAST core_nt reference corpus.'
    )
    grouped_parser.add_argument(
        '--mmseqs-sensitivity',
        type=float,
        default=5.7,
        help='MMseqs2 sensitivity setting, 1–7.5 (default: 5.7). '
             'Only used when --search-tool mmseqs2 is specified.'
    )
    grouped_parser.add_argument(
        '--mmseqs-search-type',
        type=int,
        default=3,
        help='MMseqs2 search type passed as --search-type to easy-search. '
             'Commonly used values: 2 (translated), 3 (nucleotide), '
             '4 (translated nucleotide backtrace). Default is 3 (nucleotide). '
             'See MMseqs2 documentation for all valid values. '
             'Only used when --search-tool mmseqs2 is specified.'
    )
    grouped_parser.add_argument(
        '--mmseqs-memory-limit',
        type=_mmseqs_memory_limit_type,
        default='400G',
        metavar='LIMIT',
        help='Maximum RAM for the MMseqs2 prefilter step, passed as '
             '--split-memory-limit to mmseqs easy-search (default: 400G). '
             'Accepts MMseqs2-style memory strings such as 64G, 400G, or 1T. '
             'Only used when --search-tool mmseqs2 is specified.'
    )
    grouped_parser.add_argument(
        '--overwrite-existing-blast',
        action='store_true',
        help='Overwrite existing search results'
    )
    grouped_parser.add_argument(
        '--skip-alignment',
        action='store_true',
        help='Skip alignment and tree building steps'
    )
    grouped_parser.add_argument(
        '--alignment-tolerance',
        type=int,
        default=50,
        help='Maximum coordinate difference for alignment consistency (default: 50)'
    )
    grouped_parser.add_argument(
        '--output-classification',
        type=Path,
        default=None,
        help='Path to output classification TSV file'
    )
    _add_log_level_argument(grouped_parser)
    
    # Relabel subcommand
    relabel_parser = subparsers.add_parser(
        'relabel',
        help='Relabel tree with taxonomic names from BLAST results',
        description='Relabel phylogenetic tree nodes with taxonomic names',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Relabel tree with genus names
  eplace relabel blast_results.txt input.treefile output.treefile --rank genus
  
  # Relabel tree with species names (genus + species)
  eplace relabel blast_results.txt input.treefile output.treefile --rank species
  
  # Relabel tree with family names
  eplace relabel blast_results.txt input.treefile output.treefile --rank family

Notes:
  - BLAST results must be in tabular format with taxonomy information
  - When using --rank species, the output will use "genus species" format
  - Tree file must be in Newick format
        """
    )
    relabel_parser.add_argument(
        'blast_output',
        type=Path,
        help='Path to BLAST output file (tabular format with taxonomy)'
    )
    relabel_parser.add_argument(
        'tree_file',
        type=Path,
        help='Path to input tree file (Newick format)'
    )
    relabel_parser.add_argument(
        'output_tree',
        type=Path,
        help='Path to output relabeled tree file'
    )
    relabel_parser.add_argument(
        '--rank',
        type=str,
        default='genus',
        choices=['phylum', 'class', 'order', 'family', 'genus', 'species'],
        help='Taxonomic rank for tree labeling (default: genus)'
    )
    relabel_parser.add_argument(
        '--blastdb-path',
        type=Path,
        default=None,
        help='Path to BLAST database directory (optional, not required for relabeling)'
    )
    _add_log_level_argument(relabel_parser)
    
    # Parse arguments
    args = parser.parse_args()
    
    # Set logging level: subcommand --log-level takes precedence over top-level
    log_level = getattr(args, 'log_level', None)
    log_level = log_level if log_level is not None else 'INFO'
    logging.getLogger().setLevel(getattr(logging, log_level))
    
    # If no command provided, show help
    if not args.command:
        parser.print_help()
        return 1
    
    # Route to appropriate command handler
    if args.command == 'download':
        return download_command(args)
    elif args.command in ('search', 'blast'):
        return blast_command(args)
    elif args.command == 'grouped':
        return grouped_command(args)
    elif args.command == 'relabel':
        return relabel_command(args)
    else:
        parser.print_help()
        return 1


if __name__ == '__main__':
    sys.exit(main())
