"""
Taxonomy extraction and sequence retrieval module.

This module provides functionality for extracting taxonomic information from BLAST results,
selecting representative sequences per taxonomic rank, and extracting sequences from databases.
"""

import os
import re
import subprocess
import logging
import sys
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Callable, Mapping
from collections import defaultdict

from .blast_analysis import BlastHit, normalize_sequence_id, _parse_nbdl_custom_header

from .placement import QueryPlacementPlan

import pytaxonkit


def _subject_id_matches(subject_id: str, target_id: str) -> bool:
    """Return True if *subject_id* refers to the same sequence as *target_id*.

    Exact equality is checked first so that non-NCBI pipe-delimited labels
    (e.g. ``sampleA|42``) are never conflated with an unrelated sequence that
    happens to share the same trailing segment.  Normalized comparison is used
    only as a fallback to handle cases where the same accession appears in
    different formats (e.g. ``gi|...|gb|HQ641676.1|`` vs ``HQ641676.1``, or a
    MAFFT ``_R_`` reverse-complement marker).
    """
    if subject_id == target_id:
        return True
    return normalize_sequence_id(subject_id) == normalize_sequence_id(target_id)

# Configure module logger
logger = logging.getLogger(__name__)

# Valid taxonomic ranks supported by the library
# 'no_rank' bypasses taxonomy lookup for custom databases without taxids
VALID_RANKS = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'no_rank']


def extract_custom_subject_id_and_taxid(
    subject_id: str,
    header: Optional[str] = None,
    custom_header_parser: Optional[Callable[[str], Tuple[str, str, str]]] = None
) -> Tuple[str, str]:
    """
    Extract the subject_id and taxid from a BLAST subject ID, optionally parsing
    a custom header format.
    
    This function handles both standard NCBI-formatted IDs and custom database headers
    (e.g., NBDL format). If a custom_header_parser is provided and a header is given,
    it uses the parser to extract the sequence ID and taxid. Otherwise, it returns
    the original subject_id and an empty taxid.
    
    Args:
        subject_id: The original BLAST subject ID (typically from BLAST output).
        header: The full FASTA header (optional), used with custom_header_parser.
        custom_header_parser: Optional parser function that takes a header string and
                             returns a tuple of (seq_id, description, taxid).
    
    Returns:
        Tuple of (canonical_subject_id, taxid).
        If no custom parser or header is provided, returns (subject_id, '').
    """
    if custom_header_parser and header:
        try:
            seq_id, description, taxid = custom_header_parser(header)
            return (seq_id, taxid)
        except Exception as e:
            logger.warning(f"Error parsing custom header: {e}. Using original subject_id.")
            return (subject_id, '')
    
    return (subject_id, '')


class TaxonomyExtractor:
    """
    Class for extracting taxonomic information from sequence IDs.
    """

    def parse_taxids(self, tax_ids: list[str]) -> dict[str, dict[str, tuple[str, str]]]:
        """
        Parse taxonomic information from the taxonomy IDs from the BLAST hits

        Args:
            tax_ids: the taxonomy IDs reported by BLAST
            
        Returns:
           dictionary containing the rank and a tuple of the taxonomy ID and the name
        """
        # make sure that duplicate taxids are removed before we look them up
        tax_ids = list(set(tax_ids))
        taxonomy_dict = {}

        # Handle empty taxid list (e.g., when using custom databases without taxonomy)
        if not tax_ids or all(tid == '' for tid in tax_ids):
            logger.warning("No valid taxonomy IDs provided. Skipping taxonomy lookup.")
            return taxonomy_dict

        # we need to get the whole lineage, and then convert it to a dict
        try:
            df = pytaxonkit.lineage(tax_ids)
        except Exception:
            logger.exception("Error retrieving taxonomic lineages")
            sys.exit(1)

        df['names'] = df['FullLineage'].str.split(';')
        df['taxids'] = df['FullLineageTaxIDs'].str.split(';')
        df['ranks'] = df['FullLineageRanks'].str.split(';')
        long_df = df.explode(['names', 'taxids', 'ranks'])
        filtered = long_df[long_df['ranks'].isin(VALID_RANKS)]

        for tid, rank, taxid, name in (
                filtered[['TaxID', 'ranks', 'taxids', 'names']]
                        .drop_duplicates()
                        .itertuples(index=False, name=None)
        ):
            tid = str(tid)
            taxid = str(taxid)

            taxonomy_dict.setdefault(tid, {})[rank] = (taxid, name)
        return taxonomy_dict
    
    def group_hits_by_query(
        self,
        hits: list[BlastHit]
    ) -> dict[str, list[BlastHit]]:
        """
        Group BLAST hits by query sequence.
        
        Args:
            hits: list of BlastHit objects
            
        Returns:
            dictionary mapping query IDs to lists of hits
        """
        grouped = defaultdict(list)
        for hit in hits:
            grouped[hit.query_id].append(hit)
        return dict(grouped)
    
    def select_representatives_by_rank(
        self,
        hits: list[BlastHit],
        rank: str,
        max_per_rank: int = 1,
        preferred_representatives: Optional[Dict[str, str]] = None
    ) -> list[BlastHit]:
        """
        Select representative sequences per taxonomic rank.
        
        When rank='no_rank', bypasses taxonomy lookup and groups by subject_id instead.
        This mode is useful for custom databases that don't have taxid information.
        
        Args:
            hits: list of BlastHit objects for a single query
            rank: Taxonomic rank for representative selection. Use 'no_rank' to group
                  by sequence ID instead of taxonomic rank (useful for custom databases).
            max_per_rank: Maximum number of representatives per rank (default: 1)
            preferred_representatives: Optional dictionary mapping rank_tid to preferred subject_id
                                       to ensure consistent representatives across queries
            
        Returns:
            list of representative BlastHit objects
        """
        
        if preferred_representatives is None:
            preferred_representatives = {}
        
        # Group hits by taxonomic rank (or by subject_id if rank='no_rank')
        rank_groups = defaultdict(list)
        
        reported_hits = set()
        for hit in hits:
            # Special handling for 'no_rank': group by subject_id instead of taxonomy
            if rank == 'no_rank':
                logger.info(
                    f"Using no_rank mode: grouping by subject_id (custom database mode)"
                )
                rank_groups[hit.subject_id].append(hit)
                continue
            
            if not hit.subject_taxonomy:
                # No taxonomy available (e.g. MMseqs2 database without taxonomy or custom database).
                # Fall back to grouping by subject_id so the hit still contributes
                # a representative rather than being silently dropped.
                logger.info(
                    f"No taxonomy for hit {hit.subject_id} (query {hit.query_id}); "
                    f"using subject_id as fallback group key"
                )
                rank_groups[hit.subject_id].append(hit)
                continue
            if rank not in hit.subject_taxonomy:
                logger.info(
                    f"We did not find {rank} in the taxonomy of {hit.query_id} which has subject taxid of {hit.subject_taxid}")
                continue
            if not hit.subject_taxonomy[rank]:
                logger.warning(
                    f"Hit {hit.subject_id} for query {hit.query_id} has no taxonomic information at rank {rank}")
                continue

            if isinstance(hit.subject_taxonomy[rank], tuple):
                # Log the first time we see each rank name
                if hit.subject_taxonomy[rank][1] not in reported_hits:
                    logger.info(f"Found a hit for {hit.query_id} at rank {rank}: {hit.subject_taxonomy[rank][1]} ({hit.subject_taxonomy[rank][0]})")
                    reported_hits.add(hit.subject_taxonomy[rank][1])
                # Add all hits with taxonomic information to rank_groups
                rank_groups[hit.subject_taxonomy[rank][1]].append(hit)
            else:
                logger.warning(f"Not really sure what {hit.subject_taxonomy[rank]} of type {type(hit.subject_taxonomy[rank])} is supposed to be")
        
        # Select best representative from each rank
        representatives = []
        for rank_key, rank_hits in rank_groups.items():
            # Check if we have a preferred representative for this rank
            preferred_subject_id = preferred_representatives.get(rank_key)
            
            if preferred_subject_id:
                # Look for the preferred representative in the current hits.
                # Try exact match first; fall back to normalized comparison to handle
                # NCBI format differences (e.g. gi|...|gb|ACC| vs ACC).
                preferred_hit = next(
                    (hit for hit in rank_hits if _subject_id_matches(hit.subject_id, preferred_subject_id)),
                    None
                )
                
                if preferred_hit:
                    # Use the preferred representative
                    logger.info(f"Reusing previously selected representative {preferred_subject_id} for rank {rank_key}")
                    representatives.append(preferred_hit)
                    continue
            
            # No preferred representative or it's not in current hits
            # Sort by bit score (best first) and select new representative
            rank_hits.sort(key=lambda h: h.bit_score, reverse=True)
            
            # Take top N representatives
            representatives.extend(rank_hits[:max_per_rank])
        
        logger.info(
            f"Selected {len(representatives)} representative sequences from {len(hits)} hits at rank '{rank}'"
        )
        
        return representatives
    

class SequenceExtractor:
    """
    Class for extracting sequences from BLAST databases.
    """
    
    def __init__(self, blastdb_path: Optional[Path] = None):
        """
        Initialize the SequenceExtractor.
        
        Args:
            blastdb_path: Path to BLAST database directory. If None, uses BLASTDB env var.
        """
        self.blastdb_path = blastdb_path
        if self.blastdb_path is None:
            blastdb_env = os.environ.get('BLASTDB')
            if blastdb_env:
                self.blastdb_path = Path(blastdb_env)
            else:
                self.blastdb_path = Path.home() / "blastdb"
    
    def check_blastdbcmd_available(self) -> bool:
        """
        Check if blastdbcmd is available in the system.
        
        Returns:
            True if blastdbcmd is available, False otherwise
        """
        try:
            result = subprocess.run(
                ['blastdbcmd', '-version'],
                capture_output=True,
                text=True,
                timeout=5
            )
            return result.returncode == 0
        except (subprocess.SubprocessError, FileNotFoundError):
            return False
    
    def extract_sequences(
        self,
        sequence_ids: list[str],
        output_fasta: Path,
        database: str = "core_nt"
    ) -> bool:
        """
        Extract sequences from BLAST database using blastdbcmd.
        
        Args:
            sequence_ids: list of sequence IDs to extract
            output_fasta: Path to output FASTA file
            database: Name of BLAST database (default: "core_nt")
            
        Returns:
            True if extraction was successful, False otherwise
            
        Raises:
            RuntimeError: If blastdbcmd is not available
        """
        if not self.check_blastdbcmd_available():
            raise RuntimeError("blastdbcmd is not available. Please install BLAST+ tools.")
        
        if not sequence_ids:
            logger.warning("No sequence IDs provided for extraction")
            return False
        
        # Build database path
        db_path = self.blastdb_path / database
        
        # Create a temporary file with sequence IDs
        id_file = output_fasta.parent / f"{output_fasta.stem}_ids.txt"
        
        try:
            with open(id_file, 'w') as f:
                for seq_id in sequence_ids:
                    f.write(f"{seq_id}\n")
            
            # Run blastdbcmd
            cmd = [
                'blastdbcmd',
                '-db', str(db_path),
                '-entry_batch', str(id_file),
                '-out', str(output_fasta)
            ]
            
            logger.info(f"Extracting {len(sequence_ids)} sequences from database")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600  # 10 minute timeout
            )
            
            if result.returncode != 0:
                logger.error(f"blastdbcmd failed with error: {result.stderr}")
                return False
            
            logger.info(f"Sequences extracted successfully to {output_fasta}")
            return True
            
        except subprocess.TimeoutExpired:
            logger.error("Sequence extraction timed out")
            return False
        except Exception as e:
            logger.error(f"Error extracting sequences (taxonomy): {e}")
            return False
        finally:
            # Clean up temporary ID file
            if id_file.exists():
                id_file.unlink()
    
    def extract_representatives_for_query(
        self,
        query_id: str,
        representative_hits: list[BlastHit],
        output_dir: Path,
        database: str = "core_nt"
    ) -> Optional[Path]:
        """
        Extract representative sequences for a single query to a FASTA file.
        
        Args:
            query_id: Query sequence identifier
            representative_hits: list of representative BlastHit objects
            output_dir: Output directory for FASTA files
            database: Name of BLAST database
            
        Returns:
            Path to output FASTA file if successful, None otherwise
        """
        if not representative_hits:
            logger.warning(f"No representative hits for query {query_id}")
            return None
        
        # Create output directory if it doesn't exist
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate output filename
        safe_query_id = query_id.replace('|', '_').replace('/', '_')
        output_fasta = output_dir / f"{safe_query_id}_representatives.fasta"
        
        # Extract sequence IDs
        sequence_ids = [hit.subject_id for hit in representative_hits]
        
        # Extract sequences
        success = self.extract_sequences(
            sequence_ids=sequence_ids,
            output_fasta=output_fasta,
            database=database
        )
        
        if success:
            return output_fasta
        else:
            return None


def rewrite_blast_hits(
    blast_hits: List[BlastHit],
    output_file: Path,
    header: bool = True) -> bool:
    """
    Rewrite the blast hits when we have annotated them

    Args:
        blast_hits: list of BlastHit objects
        output_file: the file to write to
        header: whether to include a header line in the file
    
    Returns:
        True on success
    """

    fields = [
        "query_id", "subject_id", "percent_identity", "alignment_length",
        "query_length", "subject_length", "query_start", "query_end",
        "subject_start", "subject_end", "evalue", "bit_score",
        "query_coverage", "subject_taxid", "subject_taxids",
        "subject_taxonomy"
    ]

    with open(output_file, 'w') as out:
        if header:
            print("\t".join(fields), file=out)

        for hit in blast_hits:
            print(
                "\t".join(
                    "" if getattr(hit, f) is None else str(getattr(hit, f))
                    for f in fields
                ),
                file=out
            )
    
    return True


def process_blast_results_for_taxonomy(
    blast_hits: List[BlastHit],
    output_dir: Path,
    rank: str = "genus",
    database: str = "core_nt",
    blastdb_path: Optional[Path] = None,
    custom_header_parser: Optional[Callable[[str], Tuple[str, str, str]]] = None
) -> Dict[str, Optional[Path]]:
    """
    Process BLAST hits to extract representative sequences per taxonomic rank.
    
    When rank='no_rank', skips taxonomy lookup and groups by sequence ID instead.
    This mode is useful for custom databases without taxid information.
    
    Args:
        blast_hits: list of BlastHit objects
        output_dir: Output directory for FASTA files
        rank: Taxonomic rank for representative selection. Use 'no_rank' to bypass
              taxonomy lookup and group by sequence ID (useful for custom databases).
        database: Name of BLAST database
        blastdb_path: Path to BLAST database directory
        custom_header_parser: Optional parser function for custom database headers.
                             Should take a header string and return (seq_id, description, taxid).
        
    Returns:
        dictionary mapping query IDs to output FASTA file paths
    """
    
    if rank not in VALID_RANKS:
        raise ValueError(f"Rank: {rank} is not a valid rank. It must be one of: {VALID_RANKS}")

    tax_extractor = TaxonomyExtractor()
    seq_extractor = SequenceExtractor(blastdb_path)
    
    # If custom header parser is provided, extract taxids from headers and override subject_taxid
    if custom_header_parser:
        for hit in blast_hits:
            # Try to extract taxid from the subject_id using custom parser
            # (subject_id in BLAST output may already be parsed by custom parser)
            # This is a fallback if the header wasn't available during BLAST parsing
            logger.debug(f"Custom header parser provided for hit {hit.subject_id}")
    
    # Skip taxonomy lookup entirely if rank is 'no_rank'
    if rank != 'no_rank':
        # get all the taxonomies
        subject_taxids = {hit.subject_taxid for hit in blast_hits}
        tax_dict = tax_extractor.parse_taxids(list(subject_taxids))

        # add all the ranks to all the hits
        for h in blast_hits:
            h.subject_taxonomy = tax_dict.get(h.subject_taxid)
    else:
        # In no_rank mode, we don't use taxonomy information
        logger.info("Using no_rank mode: skipping taxonomy lookup for custom database")
        for h in blast_hits:
            h.subject_taxonomy = None

    # Group hits by query
    grouped_hits = tax_extractor.group_hits_by_query(blast_hits)
    
    # Track selected representatives across queries to ensure consistency
    # Maps rank_tid -> subject_id of the selected representative
    preferred_representatives = {}
    
    # Process each query
    results = {}
    
    for query_id, query_hits in grouped_hits.items():
        logger.info(f"Processing query {query_id} with {len(query_hits)} hits")
        
        # Select representatives, preferring previously selected ones
        representatives = tax_extractor.select_representatives_by_rank(
            hits=query_hits,
            rank=rank,
            preferred_representatives=preferred_representatives
        )
        if len(representatives) == 0:
            logger.warning(f"Error: No representative sequences for {query_id} at rank {rank}")
            continue
        
        # Update the preferred representatives with newly selected ones
        for rep in representatives:
            if rank == 'no_rank':
                # In no_rank mode, use subject_id as the key
                if rep.subject_id not in preferred_representatives:
                    preferred_representatives[rep.subject_id] = rep.subject_id
                    logger.info(f"Recording {rep.subject_id} as representative for no_rank mode")
            elif (
                rep.subject_taxonomy
                and rank in rep.subject_taxonomy
                and isinstance(rep.subject_taxonomy[rank], tuple)
                and rep.subject_taxonomy[rank][1] not in preferred_representatives
            ):
                preferred_representatives[rep.subject_taxonomy[rank][1]] = rep.subject_id
                logger.info(f"Recording {rep.subject_id} as representative for rank {rep.subject_taxonomy[rank][1]}")
        
        # Create query-specific output directory
        query_output_dir = output_dir / query_id.replace('|', '_').replace('/', '_')
        
        # Extract sequences
        output_fasta = seq_extractor.extract_representatives_for_query(
            query_id=query_id,
            representative_hits=representatives,
            output_dir=query_output_dir,
            database=database
        )
        
        results[query_id] = output_fasta
    
    return results

# This is the new dictionary to help create the look up for the classification output
def make_safe_taxonomic_tree_label(hit: BlastHit, label_rank: str = "species") -> str:
    """
    Create the same safe tree/FASTA label used after trimming in alignment.py.

    This lets classification map a tree leaf label back to the original BlastHit.
    """
    taxid = hit.subject_taxid
    tax_label = hit.subject_id

    if isinstance(hit.subject_taxonomy, dict) and label_rank in hit.subject_taxonomy:
        taxid, tax_label = hit.subject_taxonomy[label_rank]

    safe_label = (
        f"{taxid}_{tax_label}"
        .replace(" ", "_")
        .replace(":", "_")
        .replace("(", "_")
        .replace(")", "_")
        .replace(",", "_")
        .replace(";", "_")
        .replace("|", "_")
        .replace("/", "_")
    )

    return safe_label

def sort_strings_and_numbers(s: str):
    """
    Extract text and numbers from strings for proper sorting.

    Args:
        s: string to extract the number from
    Returns:
        Returns:
            A tuple ``(text_part, num_part)`` that can be used as a sort key. For strings
            matching the pattern ``<non-digits><digits>``, this is the non-digit prefix
            and the trailing integer. For non-matching strings, returns ``(s, 0)``.

    """
    match = re.match(r'(\D+)(\d+)', s)
    if match:
        text_part = match.group(1)
        num_part = int(match.group(2))
        return (text_part, num_part)
    return (s, 0)

# Adding new lines to to update the classification output
DECISION_RANK_ORDER = [
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

MIN_IDENTITY_FOR_RANK = {
    "species": 99.0,
    "genus": 97.0,
    "family": 95.0,
    "order": 90.0,
    "class": 85.0,
    "phylum": 80.0,
    "domain": 0.0,
}

COMPETING_IDENTITY_WINDOW = 1.0
COMPETING_COVERAGE_WINDOW = 5.0


def _safe_float(value):
    try:
        return float(value)
    except Exception:
        return 0.0


def _get_taxon_at_rank(hit: BlastHit, rank: str):
    if hit.subject_taxonomy and rank in hit.subject_taxonomy:
        return hit.subject_taxonomy[rank]
    return ("N/A", "N/A")

def _safe_tree_label(label: str):
    """
    Match the label cleaning used for tree tips.
    """
    if label is None:
        return ""

    return (
        str(label)
        .replace(" ", "_")
        .replace(":", "_")
        .replace("(", "_")
        .replace(")", "_")
        .replace(",", "_")
        .replace(";", "_")
        .replace("|", "_")
        .replace("/", "_")
    )


def _split_newick_top_level(text: str):
    """
    Split a Newick subtree string on commas that are not inside parentheses.
    """
    parts = []
    depth = 0
    start = 0

    for i, char in enumerate(text):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == "," and depth == 0:
            parts.append(text[start:i])
            start = i + 1

    parts.append(text[start:])
    return [p.strip() for p in parts if p.strip()]


def _strip_branch_length(newick_part: str):
    """
    Remove a terminal branch length from a Newick token.
    """
    depth = 0

    for i, char in enumerate(newick_part):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == ":" and depth == 0:
            return newick_part[:i]

    return newick_part


def _collect_tip_labels(newick_part: str):
    """
    Collect tip labels from a simple Newick subtree.
    """
    labels = []
    token = ""

    for char in newick_part:
        if char in "(),:;":
            if token:
                labels.append(token)
                token = ""
        else:
            token += char

    if token:
        labels.append(token)

    # Remove obvious numeric branch support labels.
    cleaned = []
    for label in labels:
        label = label.strip()
        if not label:
            continue
        try:
            float(label)
            continue
        except ValueError:
            cleaned.append(label)

    return cleaned


def _find_query_sister_tip_labels(tree_file: Path, query_id: str):
    """
    Find labels in the smallest immediate sister neighbourhood around query_id.

    This is intentionally simple for Milestone 2. It finds the smallest
    parenthesized Newick subtree containing the query, then returns all other
    tip labels in that subtree.
    """
    if not tree_file or not tree_file.exists():
        return [], "No tree", "No phylogenetic tree was available."

    tree_text = tree_file.read_text().strip()

    query_label = _safe_tree_label(query_id)

    if query_id not in tree_text and query_label not in tree_text:
        return [], "Query not found", "Query was not found in the phylogenetic tree."

    search_label = query_id if query_id in tree_text else query_label
    query_pos = tree_text.find(search_label)

    # Walk left to find the nearest opening parenthesis that contains the query.
    left = query_pos
    while left >= 0:
        if tree_text[left] == "(":
            break
        left -= 1

    if left < 0:
        return [], "No reference sister clade", "Could not identify a parent subtree for the query."

    # Walk right to the matching closing parenthesis.
    depth = 0
    right = None
    for i in range(left, len(tree_text)):
        if tree_text[i] == "(":
            depth += 1
        elif tree_text[i] == ")":
            depth -= 1
            if depth == 0:
                right = i
                break

    if right is None:
        return [], "No reference sister clade", "Could not identify a complete parent subtree for the query."

    subtree = tree_text[left + 1:right]
    tips = _collect_tip_labels(subtree)

    sister_tips = [
        tip for tip in tips
        if tip not in {query_id, query_label, f"_R_{query_label}", f"{query_label}_R"}
    ]

    if not sister_tips:
        return [], "No reference sister clade", "No sister reference tips were found near the query."

    return sister_tips, "Tree neighbourhood found", "A local sister neighbourhood was found for the query."


def _match_tree_label_to_hit(label: str, query_hits: List[BlastHit], tree_label_rank: str):
    """
    Match a tree tip label back to a BlastHit.
    """
    label_options = {label}

    if label.startswith("_R_"):
        label_options.add(label[3:])
    if label.endswith("_R"):
        label_options.add(label[:-2])

    for hit in query_hits:
        possible_labels = {
            hit.subject_id,
            hit.get_accession(),
            _safe_tree_label(hit.subject_id),
            _safe_tree_label(hit.get_accession()),
        }

        if hit.subject_taxonomy and tree_label_rank in hit.subject_taxonomy:
            taxid, name = hit.subject_taxonomy[tree_label_rank]
            possible_labels.add(_safe_tree_label(f"{taxid}_{name}"))

        if label_options & possible_labels:
            return hit

    return None


def make_tree_topology_evidence(
    tree_file: Optional[Path],
    query_id: str,
    query_hits: List[BlastHit],
    tree_label_rank: str,
):
    """
    Summarise topology-based evidence from the local tree neighbourhood.
    """
    evidence = {
        "tree_nearest_neighbor_label": "N/A",
        "tree_sister_reference_count": 0,
        "tree_sister_taxa": "N/A",
        "tree_lowest_consistent_rank": "N/A",
        "tree_lowest_consistent_taxid": "N/A",
        "tree_lowest_consistent_name": "N/A",
        "tree_topology_status": "No tree",
        "tree_topology_reason": "No phylogenetic tree was available.",
    }

    if not tree_file:
        return evidence

    sister_labels, status, reason = _find_query_sister_tip_labels(tree_file, query_id)
    evidence["tree_topology_status"] = status
    evidence["tree_topology_reason"] = reason

    sister_hits = []
    first_reference_label = "N/A"

    for label in sister_labels:
        hit = _match_tree_label_to_hit(label, query_hits, tree_label_rank)
        if hit:
            sister_hits.append(hit)
            if first_reference_label == "N/A":
                first_reference_label = label

    evidence["tree_nearest_neighbor_label"] = first_reference_label

    if not sister_hits:
        evidence["tree_topology_status"] = "No reference sister clade"
        evidence["tree_topology_reason"] = (
            "Tree neighbours were found, but none could be mapped back to retained BLAST reference hits."
        )
        return evidence

    evidence["tree_sister_reference_count"] = len(sister_hits)

    taxa = []
    for hit in sister_hits:
        taxid, name = _get_taxon_at_rank(hit, tree_label_rank)
        if name != "N/A":
            taxa.append(name)
        else:
            taxa.append(hit.subject_id)

    evidence["tree_sister_taxa"] = "; ".join(sorted(set(taxa)))

    lca_rank, lca_taxid, lca_name = _find_lca_rank(sister_hits)
    evidence["tree_lowest_consistent_rank"] = lca_rank
    evidence["tree_lowest_consistent_taxid"] = lca_taxid
    evidence["tree_lowest_consistent_name"] = lca_name

    if lca_rank == "N/A":
        evidence["tree_topology_status"] = "Mixed topology"
        evidence["tree_topology_reason"] = (
            f"The query's sister neighbourhood contains {len(sister_hits)} mapped reference tips, "
            "but they do not share a consistent taxonomic rank."
        )
    elif len(sister_hits) == 1:
        evidence["tree_topology_status"] = "Single reference neighbour"
        evidence["tree_topology_reason"] = (
            f"The query's sister neighbourhood contains one mapped reference tip: "
            f"{evidence['tree_sister_taxa']}."
        )
    else:
        evidence["tree_topology_status"] = "Consistent topology"
        evidence["tree_topology_reason"] = (
            f"The query's sister neighbourhood contains {len(sister_hits)} mapped reference tips "
            f"sharing {lca_rank}: {lca_name}."
        )

    return evidence

def _format_optional_float(value, digits: int = 3):
    """
    Format numeric values for classification.tsv while preserving N/A.
    """
    if value in (None, "N/A", ""):
        return "N/A"

    try:
        return f"{float(value):.{digits}f}"
    except Exception:
        return "N/A"


def _lineage_for_hit(hit: Optional[BlastHit]):
    """
    Return the standard semicolon-delimited lineage string for a BLAST hit.
    """
    if not hit or not hit.subject_taxonomy:
        return ";;;;;"

    return ";".join([
        hit.subject_taxonomy[r][1] if r in hit.subject_taxonomy else ""
        for r in VALID_RANKS
        if r != "no_rank"
    ])


def _annotate_hits_with_taxonomy_if_needed(hits: List[BlastHit]):
    """
    Attach taxonomy dictionaries to raw/unfiltered hits when taxids are available.
    """
    missing_taxids = {
        hit.subject_taxid
        for hit in hits
        if not hit.subject_taxonomy
        and hit.subject_taxid not in (None, "", "0", "N/A")
    }

    if not missing_taxids:
        return

    tax_extractor = TaxonomyExtractor()
    tax_dict = tax_extractor.parse_taxids(list(missing_taxids))

    for hit in hits:
        if not hit.subject_taxonomy:
            hit.subject_taxonomy = tax_dict.get(hit.subject_taxid)


def make_raw_blast_evidence(
    raw_hits: List[BlastHit],
    retained_hits: List[BlastHit],
):
    """
    Summarise raw search evidence for one query.

    This distinguishes true no-hit reads from reads that had weak/subthreshold
    BLAST/MMseqs evidence but no retained hits after filtering.
    """
    evidence = {
        "raw_search_hits": 0,
        "raw_best_percent_identity": "N/A",
        "raw_best_query_coverage": "N/A",
        "raw_best_bit_score": "N/A",
        "raw_best_subject_id": "N/A",
        "raw_best_taxonomy": ";;;;;",
        "raw_evidence_status": "no_detectable_match",
    }

    if not raw_hits:
        return evidence

    best_hit = max(
        raw_hits,
        key=lambda h: (h.bit_score, h.percent_identity, h.query_coverage),
    )

    evidence.update({
        "raw_search_hits": len(raw_hits),
        "raw_best_percent_identity": _format_optional_float(best_hit.percent_identity),
        "raw_best_query_coverage": _format_optional_float(best_hit.query_coverage),
        "raw_best_bit_score": _format_optional_float(best_hit.bit_score),
        "raw_best_subject_id": best_hit.subject_id,
        "raw_best_taxonomy": _lineage_for_hit(best_hit),
        "raw_evidence_status": "retained_match" if retained_hits else "subthreshold_match",
    })

    return evidence


def _find_lca_rank(hits: List[BlastHit]):
    """
    Find the lowest shared taxonomic rank among competing hits.
    """
    if not hits:
        return ("N/A", "N/A", "N/A")

    for rank in reversed(DECISION_RANK_ORDER):
        observed = set()

        for hit in hits:
            taxid, name = _get_taxon_at_rank(hit, rank)

            if taxid == "N/A" or name == "N/A":
                observed.add(("missing", "missing"))
            else:
                observed.add((str(taxid), name))

        if len(observed) == 1:
            taxid, name = next(iter(observed))

            if taxid != "missing":
                return (rank, taxid, name)

    return ("N/A", "N/A", "N/A")


def _degrade_rank_by_identity(best_hit: BlastHit, rank: str):
    """
    If identity is too low for the LCA rank, move up to a safer rank.
    """
    if rank not in DECISION_RANK_ORDER:
        return ("N/A", "N/A", "N/A")

    best_identity = _safe_float(best_hit.percent_identity)
    rank_index = DECISION_RANK_ORDER.index(rank)

    while rank_index >= 0:
        candidate_rank = DECISION_RANK_ORDER[rank_index]
        min_identity = MIN_IDENTITY_FOR_RANK.get(candidate_rank, 0.0)

        if best_identity >= min_identity:
            taxid, name = _get_taxon_at_rank(best_hit, candidate_rank)
            return (candidate_rank, taxid, name)

        rank_index -= 1

    return ("N/A", "N/A", "N/A")


def _get_competing_hits(query_hits: List[BlastHit]):
    """
    Retain hits close enough to the top hit that they should be considered
    competing explanations for the ASV.
    """
    if not query_hits:
        return []

    best_hit = max(query_hits, key=lambda h: h.bit_score)
    best_identity = _safe_float(best_hit.percent_identity)
    best_coverage = _safe_float(best_hit.query_coverage)

    competing = []

    for hit in query_hits:
        identity = _safe_float(hit.percent_identity)
        coverage = _safe_float(hit.query_coverage)

        if (
            identity >= best_identity - COMPETING_IDENTITY_WINDOW
            and coverage >= best_coverage - COMPETING_COVERAGE_WINDOW
        ):
            competing.append(hit)

    return competing


def _count_unique_taxa(hits: List[BlastHit], rank: str):
    observed = set()

    for hit in hits:
        taxid, name = _get_taxon_at_rank(hit, rank)

        if taxid != "N/A":
            observed.add((str(taxid), name))

    return len(observed)


def make_decision_classification(
    query_hits: List[BlastHit],
    tree_best_hit: Optional[BlastHit] = None,
):
    """
    Create a decision-level classification using competing BLAST hits,
    LCA logic, identity thresholds, and tree agreement.
    """
    decision = {
        "best_percent_identity": "N/A",
        "best_query_coverage": "N/A",
        "best_bit_score": "N/A",
        "second_percent_identity": "N/A",
        "identity_gap": "N/A",
        "competing_hit_count": 0,
        "competing_species_count": 0,
        "competing_genus_count": 0,
        "competing_family_count": 0,
        "decision_rank": "N/A",
        "decision_taxid": "N/A",
        "decision_name": "N/A",
        "tree_agrees_with_decision": "N/A",
        "decision_confidence": "No classification",
        "decision_reason": "No retained BLAST hits",
    }

    if not query_hits:
        return decision

    ranked_hits = sorted(
        query_hits,
        key=lambda h: (h.bit_score, h.percent_identity, h.query_coverage),
        reverse=True,
    )

    best_hit = ranked_hits[0]
    second_hit = ranked_hits[1] if len(ranked_hits) > 1 else None

    decision["best_percent_identity"] = f"{best_hit.percent_identity:.3f}"
    decision["best_query_coverage"] = f"{best_hit.query_coverage:.3f}"
    decision["best_bit_score"] = f"{best_hit.bit_score:.3f}"

    if second_hit:
        identity_gap = best_hit.percent_identity - second_hit.percent_identity
        decision["second_percent_identity"] = f"{second_hit.percent_identity:.3f}"
        decision["identity_gap"] = f"{identity_gap:.3f}"
    else:
        decision["second_percent_identity"] = "N/A"
        decision["identity_gap"] = "N/A"

    competing_hits = _get_competing_hits(query_hits)

    decision["competing_hit_count"] = len(competing_hits)
    decision["competing_species_count"] = _count_unique_taxa(competing_hits, "species")
    decision["competing_genus_count"] = _count_unique_taxa(competing_hits, "genus")
    decision["competing_family_count"] = _count_unique_taxa(competing_hits, "family")

    lca_rank, lca_taxid, lca_name = _find_lca_rank(competing_hits)

    if lca_rank == "N/A":
        decision["decision_reason"] = "No shared LCA could be determined among competing hits"
        decision["decision_confidence"] = "Low"
        return decision

    decision_rank, decision_taxid, decision_name = _degrade_rank_by_identity(
        best_hit,
        lca_rank,
    )

    decision["decision_rank"] = decision_rank
    decision["decision_taxid"] = decision_taxid
    decision["decision_name"] = decision_name

    tree_agrees = "N/A"

    if tree_best_hit and decision_rank != "N/A":
        tree_taxid, tree_name = _get_taxon_at_rank(tree_best_hit, decision_rank)

        if str(tree_taxid) == str(decision_taxid):
            tree_agrees = "Yes"
        else:
            tree_agrees = "No"

    decision["tree_agrees_with_decision"] = tree_agrees
    
    best_species_taxid, best_species_name = _get_taxon_at_rank(best_hit, "species")
    best_genus_taxid, best_genus_name = _get_taxon_at_rank(best_hit, "genus")

    tree_species_name = "N/A"
    tree_genus_name = "N/A"

    if tree_best_hit:
        _tree_species_taxid, tree_species_name = _get_taxon_at_rank(tree_best_hit, "species")
        _tree_genus_taxid, tree_genus_name = _get_taxon_at_rank(tree_best_hit, "genus")                                                                                         
    
    reason_parts = []

    if len(competing_hits) == 1:
        reason_parts.append(
            "A single high-quality BLAST match remained after filtering, with no competing taxa meeting the similarity threshold."
        )
    else:
        reason_parts.append(
            f"{len(competing_hits)} competing BLAST matches remained after filtering."
        )
        reason_parts.append(
            f"These competing hits represent {decision['competing_species_count']} species, "
            f"{decision['competing_genus_count']} genera, and "
            f"{decision['competing_family_count']} families."
        )

    reason_parts.append(
        f"The best alignment matched {best_species_name} with "
        f"{best_hit.percent_identity:.2f}% identity and "
        f"{best_hit.query_coverage:.2f}% query coverage."
    )

    if lca_rank == decision_rank:
        reason_parts.append(
            f"The lowest common ancestor of the competing hits supports a "
            f"{decision_rank}-level classification: {decision_name}."
        )
    else:
        threshold = MIN_IDENTITY_FOR_RANK.get(lca_rank, "N/A")
        reason_parts.append(
            f"The competing hits resolved to {lca_rank}, but the best alignment identity "
            f"({best_hit.percent_identity:.2f}%) did not meet the {lca_rank}-level "
            f"confidence threshold ({threshold}%). The classification was therefore "
            f"reported conservatively at {decision_rank}: {decision_name}."
        )

    if tree_agrees == "Yes":
        reason_parts.append(
            f"The phylogenetic nearest neighbour supports this decision at the "
            f"{decision_rank} level."
        )
    elif tree_agrees == "No":
        if tree_best_hit:
            reason_parts.append(
                f"The phylogenetic nearest neighbour was {tree_species_name} "
                f"({tree_genus_name}), which does not support the "
                f"{decision_rank}-level decision of {decision_name}. "
                f"This indicates disagreement between pairwise alignment and "
                f"phylogenetic placement."
            )
        else:
            reason_parts.append(
                "The phylogenetic nearest neighbour did not support this decision."
            )
    else:
        reason_parts.append(
            "No phylogenetic nearest-neighbour support was available."
        )

    decision["decision_reason"] = " ".join(reason_parts)

    if tree_agrees == "No":
        decision["decision_confidence"] = "Low"
    elif (
        decision_rank == "species"
        and decision["competing_species_count"] == 1
        and tree_agrees in {"Yes", "N/A"}
    ):
        decision["decision_confidence"] = "High"
    elif tree_agrees == "Yes":
        decision["decision_confidence"] = "Moderate"
    elif decision_rank in {"genus", "family"}:
        decision["decision_confidence"] = "Moderate"
    else:
        decision["decision_confidence"] = "Low"

    return decision

def generate_classification_summary(
    sequences: dict[str, str],
    blast_hits: List[BlastHit],
    output_file: Path,
    rank: str = "genus",
    group_rank: str = "class",
    tree_label_rank: str = "genus",
    tree_files: Optional[dict[str, Path]] = None,
    raw_blast_hits: Optional[List[BlastHit]] = None,
    placement_plan: Optional[Mapping[str, QueryPlacementPlan]] = None,
) -> bool:
    """
    Generate a classification summary TSV file for each query sequence.
    
    This function creates a TSV file that reports:
    - Query sequence ID
    - Closest organism at the classification rank (--rank)
    - Closest organism at the grouping rank (--group-rank)
    - Closest organism at the tree labeling rank (--tree-label-rank)
    - Whether the sequence appears in multiple groups
    - Whether the sequence has no appropriate classification
    
    When rank='no_rank', taxonomy-based classification is skipped and sequences
    are identified by their sequence ID instead.
    
    The classification is based on the phylogenetically nearest neighbor in the tree
    (if available), otherwise falls back to the best BLAST hit by bit score.
    
    Args:
        sequences: dictionary of sequences that we read from the fasta file
        blast_hits: List of BlastHit objects with taxonomy information
        output_file: Path to output TSV file
        rank: Taxonomic rank for classification (default: genus). Use 'no_rank' to skip taxonomy lookup.
        group_rank: Taxonomic rank for grouping (default: class)
        tree_label_rank: Taxonomic rank for tree labeling (default: genus)
        tree_files: Optional dict mapping query_id to tree file paths for finding nearest neighbors
        
    Returns:
        True if successful, False otherwise
    """
    logger.info(f"Generating classification summary TSV to {output_file}")
    
    # Validate ranks
    for r, r_name in [(rank, 'rank'), (group_rank, 'group_rank'), (tree_label_rank, 'tree_label_rank')]:
        if r not in VALID_RANKS:
            logger.error(f"{r_name}: {r} is not a valid rank. It must be one of: {VALID_RANKS}")
            return False
    
    if raw_blast_hits is None:
        raw_blast_hits = blast_hits

    _annotate_hits_with_taxonomy_if_needed(raw_blast_hits)

    # Group retained hits by query
    query_hits_map = defaultdict(list)
    for hit in blast_hits:
        query_hits_map[hit.query_id].append(hit)

    # Group raw/unfiltered hits by query
    raw_hits_map = defaultdict(list)
    for hit in raw_blast_hits:
        raw_hits_map[hit.query_id].append(hit)
    
    # Collect all query IDs that were searched
    all_query_ids = set(sequences.keys())
    
    # Prepare data for each query
    summary_data = []
    
    for query_id in sorted(all_query_ids, key=sort_strings_and_numbers):
        query_hits = query_hits_map.get(query_id, [])
        raw_query_hits = raw_hits_map.get(query_id, [])
        
        # Initialize classification info
        classification = {
            'query_id': query_id,
            'blast_hits': 0,
            'taxonomy_blast': ';;;;;',
            'blast_classification_rank': rank,
            'blast_classification_taxid': 'N/A',
            'blast_classification_name': 'N/A',
            'blast_group_rank': group_rank,
            'blast_group_taxid': 'N/A',
            'blast_group_name': 'N/A',
            'blast_tree_label_rank': tree_label_rank,
            'blast_tree_label_taxid': 'N/A',
            'blast_tree_label_name': 'N/A',
            'taxonomy_tree': ';;;;;',
            'tree_classification_rank': rank,
            'tree_classification_taxid': 'N/A',
            'tree_classification_name': 'N/A',
            'tree_group_rank': group_rank,
            'tree_group_taxid': 'N/A',
            'tree_group_name': 'N/A',
            'tree_tree_label_rank': tree_label_rank,
            'tree_tree_label_taxid': 'N/A',
            'tree_tree_label_name': 'N/A',
            'tree_based_classification': 'No',
            'best_percent_identity': 'N/A',
            'best_query_coverage': 'N/A',
            'best_bit_score': 'N/A',
            'second_percent_identity': 'N/A',
            'identity_gap': 'N/A',
            'competing_hit_count': 0,
            'competing_species_count': 0,
            'competing_genus_count': 0,
            'competing_family_count': 0,
            'decision_rank': 'N/A',
            'decision_taxid': 'N/A',
            'decision_name': 'N/A',
            'tree_agrees_with_decision': 'N/A',
            'decision_confidence': 'No classification',
            'decision_reason': 'No retained BLAST hits',
            'tree_nearest_neighbor_label': 'N/A',
            'tree_sister_reference_count': 0,
            'tree_sister_taxa': 'N/A',
            'tree_lowest_consistent_rank': 'N/A',
            'tree_lowest_consistent_taxid': 'N/A',
            'tree_lowest_consistent_name': 'N/A',
            'tree_topology_status': 'No tree',
            'tree_topology_reason': 'No phylogenetic tree was available.',
            'raw_search_hits': 0,
            'raw_best_percent_identity': 'N/A',
            'raw_best_query_coverage': 'N/A',
            'raw_best_bit_score': 'N/A',
            'raw_best_subject_id': 'N/A',
            'raw_best_taxonomy': ';;;;;',
            'raw_evidence_status': 'no_detectable_match',
            'placement_route': 'no_evidence',
            'placement_hit_count': 0,
            'placement_best_percent_identity': 'N/A',
            'placement_best_query_coverage': 'N/A',
            'placement_best_bit_score': 'N/A',
            'placement_best_subject_id': 'N/A',
            'placement_best_taxon_name': 'N/A',
            'phylogenetic_placement_attempted': 'No',
            'phylogenetic_placement_reason': 'No placement plan was available.',
            'appears_in_multiple_groups': 'No',
            'has_classification': 'Yes'
        }
        
        raw_evidence = make_raw_blast_evidence(
            raw_hits=raw_query_hits,
            retained_hits=query_hits,
        )
        classification.update(raw_evidence)
        
        placement_query_hits = []
        if placement_plan and query_id in placement_plan:
            query_placement_plan = placement_plan[query_id]
            classification.update(query_placement_plan.to_summary_fields())
            placement_query_hits = query_placement_plan.tree_hits

        # Hits used for tree interpretation. Retained hits remain the evidence
        # for strict BLAST classification, but placement hits can still support
        # phylogenetic neighbourhood reporting.
        tree_context_hits = query_hits if query_hits else placement_query_hits
            
        if not query_hits:
            # No retained hits for this query
            classification['has_classification'] = 'No'

            if raw_query_hits:
                classification['decision_confidence'] = 'No retained classification'
                classification['decision_reason'] = (
                    "No BLAST hits passed the retained-hit filtering thresholds, "
                    f"but {classification['raw_search_hits']} raw search hits were present. "
                    f"The best raw hit was {classification['raw_best_subject_id']} "
                    f"with {classification['raw_best_percent_identity']}% identity and "
                    f"{classification['raw_best_query_coverage']}% query coverage. "
                    "Treat this as weak/subthreshold database evidence rather than a complete absence of signal."
                )

            # Do not continue yet if placement/tree context exists.
            if not tree_context_hits:
                summary_data.append(classification)
                continue

        classification['blast_hits'] = len(query_hits)

        # Extract BLAST-based classification only when retained hits exist.
        # Placement-hit queries may have tree_context_hits but zero retained hits,
        # so do not call max(query_hits) unless query_hits is non-empty.
        blast_best_hit = None
        blast_missing_ranks = []

        if query_hits:
            # Get the best retained BLAST hit (highest bit score) for strict BLAST-based classification
            blast_best_hit = max(query_hits, key=lambda h: h.bit_score)
            
            # Populate BLAST-based classification
            if blast_best_hit.subject_taxonomy:
                classification['taxonomy_blast'] = ';'.join([
                    blast_best_hit.subject_taxonomy[r][1]
                    if r in blast_best_hit.subject_taxonomy else ""
                    for r in VALID_RANKS
                    if r != 'no_rank'
                ])
            elif rank == 'no_rank':
                # In no_rank mode, use subject_id instead of taxonomy
                classification['taxonomy_blast'] = blast_best_hit.subject_id
                classification['blast_classification_taxid'] = blast_best_hit.subject_id
                classification['blast_classification_name'] = blast_best_hit.subject_id
            
            # Skip taxonomy-based classification when in no_rank mode
            if rank != 'no_rank':
                if blast_best_hit.subject_taxonomy and rank in blast_best_hit.subject_taxonomy:
                    taxid, name = blast_best_hit.subject_taxonomy[rank]
                    classification['blast_classification_taxid'] = taxid
                    classification['blast_classification_name'] = name
                else:
                    blast_missing_ranks.append(rank)
                
                if blast_best_hit.subject_taxonomy and group_rank in blast_best_hit.subject_taxonomy:
                    taxid, name = blast_best_hit.subject_taxonomy[group_rank]
                    classification['blast_group_taxid'] = taxid
                    classification['blast_group_name'] = name
                else:
                    blast_missing_ranks.append(group_rank)
                
                if blast_best_hit.subject_taxonomy and tree_label_rank in blast_best_hit.subject_taxonomy:
                    taxid, name = blast_best_hit.subject_taxonomy[tree_label_rank]
                    classification['blast_tree_label_taxid'] = taxid
                    classification['blast_tree_label_name'] = name
                else:
                    blast_missing_ranks.append(tree_label_rank)
        else:
            # No retained BLAST hits. This query may still have placement hits,
            # but it does not have a strict BLAST classification.
            if rank != 'no_rank':
                blast_missing_ranks.extend([rank, group_rank, tree_label_rank])
        
        # Try to get tree-based classification
        tree_best_hit = None
        
        if tree_files and query_id in tree_files:
            # Try to find the nearest neighbor in the phylogenetic tree
            tree_file = tree_files[query_id]
            if tree_file and tree_file.exists():
                # Import here to avoid circular dependency
                from .alignment import find_nearest_neighbor_in_tree
                
                nearest_neighbor = find_nearest_neighbor_in_tree(tree_file, query_id)
                
                tree_topology_evidence = make_tree_topology_evidence(
                    tree_file=tree_file,
                    query_id=query_id,
                    query_hits=tree_context_hits,
                    tree_label_rank=tree_label_rank,
                )
                classification.update(tree_topology_evidence)
                
                if nearest_neighbor:
                    # First try matching the new taxonomic tree labels created in alignment.py.
                    tree_label_map = {}

                    for hit in tree_context_hits:
                        safe_tree_label = make_safe_taxonomic_tree_label(
                            hit,
                            label_rank=tree_label_rank
                        )
                        tree_label_map[safe_tree_label] = hit

                        # Also handle MAFFT reverse-complement labels.
                        tree_label_map[f"_R_{safe_tree_label}"] = hit
                        tree_label_map[f"{safe_tree_label}_R"] = hit

                    tree_best_hit = tree_label_map.get(nearest_neighbor)

                    # Fall back to old subject-id matching for legacy trees or NCBI-style runs.
                    if not tree_best_hit:
                        for hit in tree_context_hits:
                            if _subject_id_matches(hit.subject_id, nearest_neighbor):
                                tree_best_hit = hit
                                break

                    if tree_best_hit:
                        classification['tree_based_classification'] = 'Yes'
                        logger.info(f"Tree-based nearest neighbor for {query_id}: {nearest_neighbor}")
                    else:
                        logger.debug(
                            f"Tree nearest neighbor {nearest_neighbor} not found in BLAST hits for {query_id}"
                        )
        
        # Populate tree-based classification if available
        if tree_best_hit:
            if tree_best_hit.subject_taxonomy:
                classification['taxonomy_tree'] = ';'.join([tree_best_hit.subject_taxonomy[r][1] if r in tree_best_hit.subject_taxonomy else ""
                                                       for r in VALID_RANKS if r != 'no_rank'])
            elif rank == 'no_rank':
                classification['taxonomy_tree'] = tree_best_hit.subject_id
                classification['tree_classification_taxid'] = tree_best_hit.subject_id
                classification['tree_classification_name'] = tree_best_hit.subject_id
            
            tree_missing_ranks = []
            
            if rank != 'no_rank':
                if tree_best_hit.subject_taxonomy and rank in tree_best_hit.subject_taxonomy:
                    taxid, name = tree_best_hit.subject_taxonomy[rank]
                    classification['tree_classification_taxid'] = taxid
                    classification['tree_classification_name'] = name
                else:
                    tree_missing_ranks.append(rank)
                
                if tree_best_hit.subject_taxonomy and group_rank in tree_best_hit.subject_taxonomy:
                    taxid, name = tree_best_hit.subject_taxonomy[group_rank]
                    classification['tree_group_taxid'] = taxid
                    classification['tree_group_name'] = name
                else:
                    tree_missing_ranks.append(group_rank)
                
                if tree_best_hit.subject_taxonomy and tree_label_rank in tree_best_hit.subject_taxonomy:
                    taxid, name = tree_best_hit.subject_taxonomy[tree_label_rank]
                    classification['tree_tree_label_taxid'] = taxid
                    classification['tree_tree_label_name'] = name
                else:
                    tree_missing_ranks.append(tree_label_rank)
        
        # Set classification status based on BLAST missing ranks
        if rank != 'no_rank' and blast_missing_ranks:
            if len(blast_missing_ranks) == 3:  # All three ranks missing
                classification['has_classification'] = 'No'
            else:
                classification['has_classification'] = 'Partial'
        
        # Check if sequence appears in multiple groups at the group_rank level
        if rank != 'no_rank':
            group_names = set()
            group_taxids = set()
            for hit in query_hits:
                if hit.subject_taxonomy and group_rank in hit.subject_taxonomy:
                    taxid, name = hit.subject_taxonomy[group_rank]
                    group_names.add(name)
                    group_taxids.add(str(taxid))

            if len(group_names) > 1:
                classification['appears_in_multiple_groups'] = 'Yes'
                # Update BLAST group names/taxids to show all groups
                classification['blast_group_name'] = '; '.join(sorted(group_names))
                classification['blast_group_taxid'] = '; '.join(sorted(group_taxids))
        
        decision = make_decision_classification(
            query_hits=query_hits,
            tree_best_hit=tree_best_hit,
        )

        classification.update(decision)
        
        summary_data.append(classification)
    
    # Write TSV file
    try:
        with open(output_file, 'w') as f:
            # Write header with both BLAST and tree-based columns
            headers = [
                'query_id',
                'blast_hits',
                'taxonomy_blast',
                'blast_classification_rank',
                'blast_classification_taxid',
                'blast_classification_name',
                'blast_group_rank',
                'blast_group_taxid',
                'blast_group_name',
                'blast_tree_label_rank',
                'blast_tree_label_taxid',
                'blast_tree_label_name',
                'tree_based_classification',
                'taxonomy_tree',
                'tree_classification_rank',
                'tree_classification_taxid',
                'tree_classification_name',
                'tree_group_rank',
                'tree_group_taxid',
                'tree_group_name',
                'tree_tree_label_rank',
                'tree_tree_label_taxid',
                'tree_tree_label_name',
                'appears_in_multiple_groups',
                'best_percent_identity',
                'best_query_coverage',
                'best_bit_score',
                'second_percent_identity',
                'identity_gap',
                'competing_hit_count',
                'competing_species_count',
                'competing_genus_count',
                'competing_family_count',
                'decision_rank',
                'decision_taxid',
                'decision_name',
                'tree_agrees_with_decision',
                'decision_confidence',
                'decision_reason',
                'tree_nearest_neighbor_label',
                'tree_sister_reference_count',
                'tree_sister_taxa',
                'tree_lowest_consistent_rank',
                'tree_lowest_consistent_taxid',
                'tree_lowest_consistent_name',
                'tree_topology_status',
                'tree_topology_reason',
                'raw_search_hits',
                'raw_best_percent_identity',
                'raw_best_query_coverage',
                'raw_best_bit_score',
                'raw_best_subject_id',
                'raw_best_taxonomy',
                'raw_evidence_status',
                'placement_route',
                'placement_hit_count',
                'placement_best_percent_identity',
                'placement_best_query_coverage',
                'placement_best_bit_score',
                'placement_best_subject_id',
                'placement_best_taxon_name',
                'phylogenetic_placement_attempted',
                'phylogenetic_placement_reason',
                'has_classification'
            ]
            f.write('\t'.join(headers) + '\n')
            
            # Write data
            for entry in summary_data:
                row = [
                    entry['query_id'],
                    str(entry['blast_hits']),
                    entry['taxonomy_blast'],
                    entry['blast_classification_rank'],
                    entry['blast_classification_taxid'],
                    entry['blast_classification_name'],
                    entry['blast_group_rank'],
                    entry['blast_group_taxid'],
                    entry['blast_group_name'],
                    entry['blast_tree_label_rank'],
                    entry['blast_tree_label_taxid'],
                    entry['blast_tree_label_name'],
                    entry['tree_based_classification'],
                    entry['taxonomy_tree'],
                    entry['tree_classification_rank'],
                    entry['tree_classification_taxid'],
                    entry['tree_classification_name'],
                    entry['tree_group_rank'],
                    entry['tree_group_taxid'],
                    entry['tree_group_name'],
                    entry['tree_tree_label_rank'],
                    entry['tree_tree_label_taxid'],
                    entry['tree_tree_label_name'],
                    entry['appears_in_multiple_groups'],
                    str(entry['best_percent_identity']),
                    str(entry['best_query_coverage']),
                    str(entry['best_bit_score']),
                    str(entry['second_percent_identity']),
                    str(entry['identity_gap']),
                    str(entry['competing_hit_count']),
                    str(entry['competing_species_count']),
                    str(entry['competing_genus_count']),
                    str(entry['competing_family_count']),
                    entry['decision_rank'],
                    str(entry['decision_taxid']),
                    entry['decision_name'],
                    entry['tree_agrees_with_decision'],
                    entry['decision_confidence'],
                    entry['decision_reason'],
                    entry['tree_nearest_neighbor_label'],
                    str(entry['tree_sister_reference_count']),
                    entry['tree_sister_taxa'],
                    entry['tree_lowest_consistent_rank'],
                    str(entry['tree_lowest_consistent_taxid']),
                    entry['tree_lowest_consistent_name'],
                    entry['tree_topology_status'],
                    entry['tree_topology_reason'],
                    str(entry['raw_search_hits']),
                    entry['raw_best_percent_identity'],
                    entry['raw_best_query_coverage'],
                    entry['raw_best_bit_score'],
                    entry['raw_best_subject_id'],
                    entry['raw_best_taxonomy'],
                    entry['raw_evidence_status'],
                    entry['placement_route'],
                    str(entry['placement_hit_count']),
                    entry['placement_best_percent_identity'],
                    entry['placement_best_query_coverage'],
                    entry['placement_best_bit_score'],
                    entry['placement_best_subject_id'],
                    entry['placement_best_taxon_name'],
                    entry['phylogenetic_placement_attempted'],
                    entry['phylogenetic_placement_reason'],
                    entry['has_classification']
                ]
                f.write('\t'.join(row) + '\n')
        
        logger.info(f"Successfully wrote classification summary for {len(summary_data)} queries to {output_file}")
        return True
        
    except Exception as e:
        logger.exception(f"Error writing classification summary TSV: {e}")
        return False
