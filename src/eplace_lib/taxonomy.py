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
from typing import Optional, List, Dict, Tuple, Callable
from collections import defaultdict

from .blast_analysis import BlastHit, normalize_sequence_id, _parse_nbdl_custom_header

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

def generate_classification_summary(
    sequences: dict[str, str],
    blast_hits: List[BlastHit],
    output_file: Path,
    rank: str = "genus",
    group_rank: str = "class",
    tree_label_rank: str = "genus",
    tree_files: Optional[dict[str, Path]] = None
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
    
    # Group hits by query
    query_hits_map = defaultdict(list)
    for hit in blast_hits:
        query_hits_map[hit.query_id].append(hit)
    
    # Collect all query IDs that were searched
    all_query_ids = set(sequences.keys())
    
    # Prepare data for each query
    summary_data = []
    
    for query_id in sorted(all_query_ids, key=sort_strings_and_numbers):
        query_hits = query_hits_map.get(query_id, [])
        
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
            'appears_in_multiple_groups': 'No',
            'has_classification': 'Yes'
        }
        
        if not query_hits:
            # No hits for this query
            classification['has_classification'] = 'No'
            summary_data.append(classification)
            continue

        classification['blast_hits'] = len(query_hits)

        # Get the best BLAST hit (highest bit score) for BLAST-based classification
        blast_best_hit = max(query_hits, key=lambda h: h.bit_score)
        
        # Populate BLAST-based classification
        if blast_best_hit.subject_taxonomy:
            classification['taxonomy_blast'] = ';'.join([blast_best_hit.subject_taxonomy[r][1] if r in blast_best_hit.subject_taxonomy else ""
                                                   for r in VALID_RANKS if r != 'no_rank'])
        elif rank == 'no_rank':
            # In no_rank mode, use subject_id instead of taxonomy
            classification['taxonomy_blast'] = blast_best_hit.subject_id
            classification['blast_classification_taxid'] = blast_best_hit.subject_id
            classification['blast_classification_name'] = blast_best_hit.subject_id
        
        # Extract BLAST-based classification at different ranks
        blast_missing_ranks = []
        
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
        
        # Try to get tree-based classification
        tree_best_hit = None
        
        if tree_files and query_id in tree_files:
            # Try to find the nearest neighbor in the phylogenetic tree
            tree_file = tree_files[query_id]
            if tree_file and tree_file.exists():
                # Import here to avoid circular dependency
                from .alignment import find_nearest_neighbor_in_tree
                
                nearest_neighbor = find_nearest_neighbor_in_tree(tree_file, query_id)
                
                if nearest_neighbor:
                    # Find the BLAST hit corresponding to the nearest neighbor.
                    # Try exact match first; fall back to normalized comparison to handle
                    # NCBI format differences and MAFFT _R_ markers.
                    for hit in query_hits:
                        if _subject_id_matches(hit.subject_id, nearest_neighbor):
                            tree_best_hit = hit
                            classification['tree_based_classification'] = 'Yes'
                            break
                    
                    if tree_best_hit:
                        logger.info(f"Tree-based nearest neighbor for {query_id}: {nearest_neighbor}")
                    else:
                        logger.debug(f"Tree nearest neighbor {nearest_neighbor} not found in BLAST hits for {query_id}")
        
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
                    entry['has_classification']
                ]
                f.write('\t'.join(row) + '\n')
        
        logger.info(f"Successfully wrote classification summary for {len(summary_data)} queries to {output_file}")
        return True
        
    except Exception as e:
        logger.error(f"Error writing classification summary TSV: {e}")
        return False
