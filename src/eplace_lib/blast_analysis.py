"""
BLAST analysis module for sequence comparison.

This module provides functionality for running BLAST searches and filtering results
based on sequence identity and coverage criteria.
"""

import os
import re
import subprocess
import logging
from typing import Optional, Dict, Tuple, Callable
from pathlib import Path
from dataclasses import dataclass

# Configure module logger
logger = logging.getLogger(__name__)


@dataclass
class BlastHit:
    """
    Represents a single BLAST hit result.
    
    Attributes:
        query_id: Query sequence identifier
        subject_id: Subject (database) sequence identifier
        percent_identity: Percentage of identical matches
        alignment_length: Length of alignment
        query_length: Length of query sequence
        subject_length: Length of subject sequence
        query_start: Start position in query
        query_end: End position in query
        subject_start: Start position in subject
        subject_end: End position in subject
        evalue: Expectation value
        bit_score: Bit score
        query_coverage: Percentage of query covered by alignment
        subject_taxonomy: The subjects taxonomy information. A dictionary with rank as key and a tuple of (taxid, name) as value.
    """
    query_id: str
    subject_id: str
    percent_identity: float
    alignment_length: int
    query_length: int
    subject_length: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    evalue: float
    bit_score: float
    query_coverage: float
    subject_taxid: str
    subject_taxids: str
    subject_taxonomy: Optional[Dict[str, Tuple[str, str]]] = None

    def get_subject_taxonomy(self, rank: str) -> Optional[tuple[str, str]]:
        """
        Return the taxonomy information as a tuple of (taxid, name) for the given rank.
        If the rank is not found, return None.

        Args:
            rank: The rank to return the taxonomy information for.
        Returns:
            tuple of (taxid, name) for the given rank, or None if the rank is not found.
        """

        if not self.subject_taxonomy:
            logger.warning("Taxonomy information is not available. Please parse it before calling this function.")
            return None

        if rank in self.subject_taxonomy:
            return self.subject_taxonomy[rank]
        else:
            return None
    
    def get_accession(self) -> str:
        """
        Extract the accession number from the subject_id.
        
        BLAST IDs can be in various formats:
        - gi|2273658778|gb|MZ387488.1| -> MZ387488.1
        - ref|NZ_CP123456.1| -> NZ_CP123456.1
        - gb|MZ387488.1| -> MZ387488.1
        - MZ387488.1 -> MZ387488.1 (already in accession format)
        
        Note: gnl|database|identifier format is handled by returning the identifier,
        but these may not be standard accessions.
        
        Returns:
            The accession number extracted from subject_id, or the full subject_id
            if no standard format is detected
        """
        return _extract_accession_from_pipe_id(self.subject_id)


def _extract_accession_from_pipe_id(seq_id: str) -> str:
    """
    Extract the accession number from a pipe-delimited NCBI-style sequence ID.

    Handles formats such as:
    - gi|2273658778|gb|MZ387488.1| -> MZ387488.1
    - ref|NZ_CP123456.1|           -> NZ_CP123456.1
    - gb|MZ387488.1|               -> MZ387488.1
    - gnl|BL_ORD_ID|12345          -> 12345
    - MZ387488.1                   -> MZ387488.1 (returned unchanged)
    - sampleA|42                   -> sampleA|42 (no known pattern, returned unchanged)

    No generic fallback is applied: if the ID contains pipes but does not match
    a known NCBI format (gnl or a standard db-prefix like gb/ref/emb/…), the
    original string is returned unchanged.  This prevents unrelated IDs that
    share the same trailing pipe-segment (e.g. ``sampleA|42`` and
    ``sampleB|42``) from being incorrectly considered equivalent.

    Args:
        seq_id: Sequence identifier, which may or may not be pipe-delimited.

    Returns:
        Extracted accession, or *seq_id* itself if no pipe-delimited structure is recognised.
    """
    if '|' in seq_id:
        parts = seq_id.split('|')

        # gnl|database|identifier format
        if len(parts) >= 3 and parts[0] == 'gnl':
            return parts[2]

        # Known database prefixes: gi|123|gb|ACC.1|, ref|ACC.1|, gb|ACC.1|, etc.
        db_identifiers = ['gb', 'ref', 'emb', 'dbj', 'pdb', 'prf', 'sp', 'tr']
        for i, part in enumerate(parts):
            if part in db_identifiers:
                if i + 1 < len(parts) and parts[i + 1]:
                    return parts[i + 1]

    return seq_id


def _parse_nbdl_custom_header(header: str) -> Tuple[str, str, str]:
    """
    Parse a custom NBDL database header to extract key components.
    
    This parser handles headers in the NBDL format:
    >NBDL-HR5SFP7MHKSYMG.v1.mt|13904-15593|+|MT-RNR2 [species=...] [dbxref=...] ... |taxid=80951
    
    Extracts:
    1. Sequence ID: Everything before the first space/bracket (e.g., NBDL-HR5SFP7MHKSYMG.v1.mt|13904-15593|+|MT-RNR2)
    2. Description: The description text between the end of sequence ID and first bracket (e.g., organism name and gene description)
    3. TaxID: The taxid value from the final |taxid=XXXXX segment
    
    Args:
        header: The full FASTA header (including '>').
    
    Returns:
        Tuple of (sequence_id, description, taxid) or (sequence_id, '', '') if parsing fails.
    """
    if not header:
        return (header, '', '')
    
    # Remove leading '>'
    header = header.lstrip('>')
    
    # Extract sequence ID (first token before space or bracket)
    seq_id_match = re.match(r'([^\s\[]+)', header)
    seq_id = seq_id_match.group(1) if seq_id_match else header
    
    # Extract taxid from the header (look for |taxid=XXXXX pattern)
    taxid = ''
    taxid_match = re.search(r'\|taxid=(\d+)', header)
    if taxid_match:
        taxid = taxid_match.group(1)
    
    # Extract description: the text between sequence ID and first bracket
    # This captures organism name and gene description
    desc_match = re.search(r'^' + re.escape(seq_id) + r'\s+([^\[]*)', header)
    description = desc_match.group(1).strip() if desc_match else ''
    
    return (seq_id, description, taxid)


def normalize_sequence_id(seq_id: str, parser: Optional[Callable[[str], str]] = None) -> str:
    """
    Normalize an arbitrary sequence or tree label to a canonical accession-like identifier.

    This is used to compare IDs from different sources (BLAST subject IDs, FASTA headers,
    tree leaf labels) that may be formatted differently but refer to the same sequence.

    Normalization steps:
    1. Strip a leading '>' (FASTA header prefix).
    2. Take only the first whitespace-delimited token.
    3. Remove MAFFT reverse-complement markers: a leading '_R_' prefix or a trailing '_R_' suffix.
    4. If the token contains pipes ('|'), extract the accession via _extract_accession_from_pipe_id()
       (gi|...|gb|ACC|, ref|ACC|, gb|ACC|, etc.).
    5. Otherwise return the token unchanged.
    
    If a custom parser is provided, it will be used instead of the default NCBI-style parsing.

    Args:
        seq_id: Raw sequence identifier from any source.
        parser: Optional custom parser function to apply instead of default parsing.

    Returns:
        Canonical accession string suitable for exact comparison.
    """
    if not seq_id:
        return seq_id

    # If a custom parser is provided, use it
    if parser:
        return parser(seq_id)

    # Step 1: strip leading '>'
    token = seq_id.lstrip('>')

    # Step 2: take first whitespace-delimited token
    parts = token.split()
    token = parts[0] if parts else token

    # Step 3: remove MAFFT reverse-complement markers
    if token.startswith('_R_'):
        token = token[3:]
    if token.endswith('_R_'):
        token = token[:-3]

    # Step 4: parse pipe-delimited NCBI-style identifiers
    return _extract_accession_from_pipe_id(token)


class FastaReader:
    """
    Class for reading FASTA files.
    """
    
    @staticmethod
    def read_fasta(fasta_path: Path) -> dict[str, str]:
        """
        Read sequences from a FASTA file.
        
        Args:
            fasta_path: Path to the FASTA file
            
        Returns:
            dictionary mapping sequence IDs to sequences
            
        Raises:
            FileNotFoundError: If FASTA file doesn't exist
            ValueError: If FASTA file is malformed
        """
        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
        
        sequences = {}
        current_id = None
        current_seq = []
        
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_id is not None:
                        sequences[current_id] = ''.join(current_seq)
                    
                    # Start new sequence
                    # Remove '>' and take first word as ID
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Save last sequence
            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)
        
        if not sequences:
            raise ValueError(f"No sequences found in FASTA file: {fasta_path}")
        
        return sequences
    
    @staticmethod
    def get_sequence_lengths(fasta_path: Path) -> dict[str, int]:
        """
        Get the length of each sequence in a FASTA file.
        
        Args:
            fasta_path: Path to the FASTA file
            
        Returns:
            dictionary mapping sequence IDs to their lengths
        """
        sequences = FastaReader.read_fasta(fasta_path)
        return {seq_id: len(seq) for seq_id, seq in sequences.items()}


class BlastRunner:
    """
    Class for running BLAST searches and parsing results.
    """
    
    def __init__(self, blastdb_path: Optional[Path] = None):
        """
        Initialize the BlastRunner.
        
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
    
    def check_blastn_available(self) -> bool:
        """
        Check if blastn is available in the system.
        
        Returns:
            True if blastn is available, False otherwise
        """
        try:
            result = subprocess.run(
                ['blastn', '-version'],
                capture_output=True,
                text=True,
                timeout=5
            )
            return result.returncode == 0
        except (subprocess.SubprocessError, FileNotFoundError):
            return False
    
    def run_blastn(
        self,
        query_fasta: Path,
        output_file: Path,
        database: str = "core_nt",
        num_threads: int = 1,
        max_target_seqs: int = 100,
        evalue: float = 1e-5,
        outfmt: str = "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore staxid staxids"
    ) -> bool:
        """
        Run blastn search.
        
        Args:
            query_fasta: Path to query FASTA file
            output_file: Path to output file
            database: Name of BLAST database (default: "core_nt")
            num_threads: Number of threads to use
            max_target_seqs: Maximum number of target sequences to report
            evalue: E-value threshold
            outfmt: Output format string
            
        Returns:
            True if BLAST ran successfully, False otherwise
            
        Raises:
            FileNotFoundError: If query file doesn't exist
            RuntimeError: If blastn is not available
        """
        if not query_fasta.exists():
            raise FileNotFoundError(f"Query FASTA file not found: {query_fasta}")
        
        if not self.check_blastn_available():
            raise RuntimeError("blastn is not available. Please install BLAST+ tools.")
        
        # Build database path
        db_path = self.blastdb_path / database
        
        # Build blastn command
        cmd = [
            'blastn',
            '-query', str(query_fasta),
            '-db', str(db_path),
            '-out', str(output_file),
            '-outfmt', outfmt,
            '-num_threads', str(num_threads),
            '-max_target_seqs', str(max_target_seqs),
            '-evalue', str(evalue)
        ]
        
        logger.info(f"Running BLAST search: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout
            )
            
            if result.returncode != 0:
                logger.error(f"BLAST failed with error: {result.stderr}")
                return False
            
            logger.info(f"BLAST search completed successfully. Output: {output_file}")
            return True
            
        except subprocess.TimeoutExpired:
            logger.error("BLAST search timed out")
            return False
        except Exception as e:
            logger.error(f"Error running BLAST: {e}")
            return False
    
    def parse_blast_results(
        self,
        blast_output: Path,
        query_lengths: Optional[dict[str, int]] = None
    ) -> list[BlastHit]:
        """
        Parse BLAST tabular output.
        
        Args:
            blast_output: Path to BLAST output file (tabular format)
            query_lengths: dictionary of query sequence lengths. If None, uses qlen from results.
            
        Returns:
            list of BlastHit objects
            
        Raises:
            FileNotFoundError: If BLAST output file doesn't exist
            ValueError: If BLAST output is malformed
        """
        if not blast_output.exists():
            raise FileNotFoundError(f"BLAST output file not found: {blast_output}")
        
        hits = []
        
        with open(blast_output, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                fields = line.split('\t')
                if len(fields) < 12:
                    raise ValueError(
                        f"Invalid BLAST output format at line {line_num}: "
                        f"expected at least 12 fields, got {len(fields)}"
                    )
                
                try:
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
                    staxid = fields[12]
                    staxids = fields[13]
                    
                    # Calculate query coverage
                    query_coverage = (abs(query_end - query_start) + 1) / query_length * 100
                    
                    hit = BlastHit(
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
                        subject_taxids=staxids
                    )
                    hits.append(hit)
                    
                except (ValueError, IndexError) as e:
                    raise ValueError(
                        f"Error parsing BLAST output at line {line_num}: {e}"
                    )
        
        return hits
    
    def filter_blast_hits(
        self,
        hits: list[BlastHit],
        min_identity: float = 90.0,
        min_coverage: float = 80.0,
        min_alignment_length: Optional[int] = None
    ) -> list[BlastHit]:
        """
        Filter BLAST hits based on identity and coverage thresholds.
        
        Args:
            hits: list of BlastHit objects
            min_identity: Minimum percent identity (default: 90.0)
            min_coverage: Minimum query coverage percentage (default: 80.0)
            min_alignment_length: Minimum alignment length (optional)
            
        Returns:
            Filtered list of BlastHit objects
        """
        filtered_hits = []
        
        for hit in hits:
            if hit.percent_identity < min_identity:
                continue
            if hit.query_coverage < min_coverage:
                continue
            if min_alignment_length and hit.alignment_length < min_alignment_length:
                continue
            
            filtered_hits.append(hit)
        
        logger.info(
            f"Filtered {len(hits)} hits to {len(filtered_hits)} hits "
            f"(min_identity={min_identity}%, min_coverage={min_coverage}%)"
        )
        
        return filtered_hits


def run_blast_search(
    query_fasta: Path,
    output_file: Path,
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
    database: str = "core_nt",
    blastdb_path: Optional[Path] = None,
    num_threads: int = 1,
    skip_existing: bool = True
) -> tuple[bool, list[BlastHit]]:
    """
    Convenience function to run BLAST search and return filtered hits.
    
    Args:
        query_fasta: Path to query FASTA file
        output_file: Path to output file
        min_identity: Minimum percent identity (default: 90.0)
        min_coverage: Minimum query coverage percentage (default: 80.0)
        database: Name of BLAST database (default: "core_nt")
        blastdb_path: Path to BLAST database directory
        num_threads: Number of threads to use
        skip_existing: Skip search if output file already exists (default: True)
        
    Returns:
        Tuple of (success: bool, filtered_hits: list[BlastHit])
    """

    runner = BlastRunner(blastdb_path)

    if os.path.exists(output_file) and skip_existing:
        logger.info(f"The blast output file {output_file} already exists. Skipping and using these results")
    else:    
        # Run BLAST
        success = runner.run_blastn(
            query_fasta=query_fasta,
            output_file=output_file,
            database=database,
            num_threads=num_threads
        )
        
        if not success:
            return False, []
        
    # Parse results
    try:
        hits = runner.parse_blast_results(output_file)
        
        # Filter results
        filtered_hits = runner.filter_blast_hits(
            hits,
            min_identity=min_identity,
            min_coverage=min_coverage
        )
        
        return True, filtered_hits
        
    except Exception as e:
        logger.error(f"Error processing BLAST results: {e}")
        return False, []


class MMseqs2Runner:
    """
    Class for running MMseqs2 searches and parsing results.

    MMseqs2 (Many-against-Many sequence searching) is an alternative to BLAST
    for sequence similarity searching, offering improved speed and sensitivity.
    Results are parsed into BlastHit objects for compatibility with the rest of
    the ePLACE pipeline.

    The target database can be either a pre-built MMseqs2 database (created with
    ``mmseqs createdb``) or a FASTA file that MMseqs2 indexes automatically.
    Taxonomy fields (taxid) are populated only when the database was built with
    taxonomy information (``mmseqs createtaxdb``); otherwise they default to "0".

    **Database selection**: To keep results comparable with the BLAST workflow
    (which uses NCBI ``core_nt``), the recommended MMseqs2 database should be
    built from the same underlying sequence collection as ``core_nt``.  This
    means creating an MMseqs2 database from the FASTA sequences that make up
    NCBI ``core_nt`` (e.g. by exporting them with ``blastdbcmd -db core_nt
    -entry all`` and then running ``mmseqs createdb``).  Using a different
    nucleotide collection will change the search space and may produce
    classification differences that reflect the database rather than the search
    algorithm.  There is no official pre-built MMseqs2 ``core_nt`` database;
    users must provide their own.
    """

    def __init__(self, db_path: Optional[Path] = None):
        """
        Initialize the MMseqs2Runner.

        Args:
            db_path: Path to the MMseqs2 database directory. If None the
                ``MMSEQS_DB_DIR`` environment variable is used; if unset,
                ``MMSEQS2DB`` is used as a legacy fallback; if both are unset,
                the directory ``~/mmseqs2db`` is used.
        """
        self.db_path = db_path
        if self.db_path is None:
            mmseqs_env = os.environ.get('MMSEQS_DB_DIR') or os.environ.get('MMSEQS2DB')
            if mmseqs_env:
                self.db_path = Path(mmseqs_env)
            else:
                self.db_path = Path.home() / "mmseqs2db"

    def check_mmseqs_available(self) -> bool:
        """
        Check if mmseqs is available in the system PATH.

        Returns:
            True if mmseqs is available, False otherwise
        """
        try:
            result = subprocess.run(
                ['mmseqs', 'version'],
                capture_output=True,
                text=True,
                timeout=5
            )
            return result.returncode == 0
        except (subprocess.SubprocessError, FileNotFoundError):
            return False

    def run_easy_search(
        self,
        query_fasta: Path,
        output_file: Path,
        database: str = "core_nt",
        num_threads: int = 1,
        max_target_seqs: int = 100,
        evalue: float = 1e-5,
        sensitivity: float = 5.7,
        tmp_dir: Optional[Path] = None,
        search_type: int = 3,
        split_memory_limit: Optional[str] = None,
        timeout: int = 3600
    ) -> bool:
        """
        Run an MMseqs2 easy-search.

        The output is written in a tab-separated format with the following
        columns (in order):
        query, target, pident, alnlen, qlen, tlen, qstart, qend, tstart, tend,
        evalue, bits, taxid, taxlineage

        Args:
            query_fasta: Path to query FASTA file
            output_file: Path to output file
            database: Name of the MMseqs2 database inside ``db_path``
                (default: "core_nt").  There is no official pre-built MMseqs2
                ``core_nt`` database; users must build their own from the same
                sequence collection as BLAST ``core_nt`` to keep results
                comparable across backends.
            num_threads: Number of threads to use
            max_target_seqs: Maximum number of target sequences to report
            evalue: E-value threshold
            sensitivity: MMseqs2 sensitivity (1–7.5, default: 5.7)
            tmp_dir: Temporary directory for MMseqs2 intermediate files.
                Defaults to a ``mmseqs_tmp`` subdirectory next to ``output_file``.
            search_type: MMseqs2 search type passed as ``--search-type`` to
                ``easy-search``. Commonly used values: 2 (translated),
                3 (nucleotide), 4 (translated nucleotide backtrace).
                Default is 3 (nucleotide). See MMseqs2 documentation for all
                valid values.
            split_memory_limit: Maximum RAM for the MMseqs2 prefilter/index
                step, passed as ``--split-memory-limit`` to ``easy-search``
                (e.g. ``"400G"``).  When ``None`` the flag is omitted and
                MMseqs2 uses its own default.
            timeout: Maximum runtime for the MMseqs2 search in seconds
                (default: 3600).

        Returns:
            True if MMseqs2 ran successfully, False otherwise

        Raises:
            FileNotFoundError: If query file doesn't exist
            RuntimeError: If mmseqs is not available
        """
        if not query_fasta.exists():
            raise FileNotFoundError(f"Query FASTA file not found: {query_fasta}")

        if not self.check_mmseqs_available():
            raise RuntimeError("mmseqs is not available. Please install MMseqs2.")

        # Build database path
        db_path = self.db_path / database

        # Set up tmp directory
        if tmp_dir is None:
            tmp_dir = output_file.parent / "mmseqs_tmp"
        tmp_dir.mkdir(parents=True, exist_ok=True)

        # Custom output format: matches the field order used in parse_mmseqs_results
        format_output = (
            "query,target,pident,alnlen,qlen,tlen,"
            "qstart,qend,tstart,tend,evalue,bits,taxid,taxlineage"
        )

        cmd = [
            'mmseqs', 'easy-search',
            str(query_fasta),
            str(db_path),
            str(output_file),
            str(tmp_dir),
            '--format-output', format_output,
            '--threads', str(num_threads),
            '--max-seqs', str(max_target_seqs),
            '-e', str(evalue),
            '-s', str(sensitivity),
            '--search-type', str(search_type)
        ]

        if split_memory_limit is not None:
            cmd.extend(['--split-memory-limit', split_memory_limit])

        logger.info(f"Running MMseqs2 search: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout
            )

            if result.returncode != 0:
                logger.error(f"MMseqs2 failed with error: {result.stderr}")
                return False

            logger.info(f"MMseqs2 search completed successfully. Output: {output_file}")
            return True

        except subprocess.TimeoutExpired:
            logger.error("MMseqs2 search timed out")
            return False
        except Exception as e:
            logger.error(f"Error running MMseqs2: {e}")
            return False

    def parse_mmseqs_results(
        self,
        mmseqs_output: Path,
        query_lengths: Optional[dict[str, int]] = None
    ) -> list[BlastHit]:
        """
        Parse MMseqs2 tabular output into BlastHit objects.

        Expects output generated with ``--format-output`` set to:
        ``query,target,pident,alnlen,qlen,tlen,qstart,qend,tstart,tend,evalue,bits,taxid,taxlineage``

        The ``taxid`` and ``taxlineage`` columns are optional; if absent or
        set to "N/A" / "0", ``subject_taxid`` will be stored as "0".

        Args:
            mmseqs_output: Path to MMseqs2 output file
            query_lengths: Unused; kept for API compatibility with
                BlastRunner.parse_blast_results.

        Returns:
            list of BlastHit objects

        Raises:
            FileNotFoundError: If the output file doesn't exist
            ValueError: If the output is malformed
        """
        if not mmseqs_output.exists():
            raise FileNotFoundError(f"MMseqs2 output file not found: {mmseqs_output}")

        hits = []

        with open(mmseqs_output, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                fields = line.split('\t')
                if len(fields) < 12:
                    raise ValueError(
                        f"Invalid MMseqs2 output format at line {line_num}: "
                        f"expected at least 12 fields, got {len(fields)}"
                    )

                try:
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

                    # taxid field (may be absent or "N/A" when database has no taxonomy)
                    taxid = "0"
                    if len(fields) > 12 and fields[12] not in ("", "N/A", "0"):
                        taxid = fields[12]

                    # Calculate query coverage
                    query_coverage = (abs(query_end - query_start) + 1) / query_length * 100

                    hit = BlastHit(
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
                        subject_taxid=taxid,
                        subject_taxids=taxid
                    )
                    hits.append(hit)

                except (ValueError, IndexError) as e:
                    raise ValueError(
                        f"Error parsing MMseqs2 output at line {line_num}: {e}"
                    )

        return hits

    def filter_hits(
        self,
        hits: list[BlastHit],
        min_identity: float = 90.0,
        min_coverage: float = 80.0,
        min_alignment_length: Optional[int] = None
    ) -> list[BlastHit]:
        """
        Filter MMseqs2 hits based on identity and coverage thresholds.

        Args:
            hits: list of BlastHit objects
            min_identity: Minimum percent identity (default: 90.0)
            min_coverage: Minimum query coverage percentage (default: 80.0)
            min_alignment_length: Minimum alignment length (optional)

        Returns:
            Filtered list of BlastHit objects
        """
        filtered_hits = []

        for hit in hits:
            if hit.percent_identity < min_identity:
                continue
            if hit.query_coverage < min_coverage:
                continue
            if min_alignment_length and hit.alignment_length < min_alignment_length:
                continue
            filtered_hits.append(hit)

        logger.info(
            f"Filtered {len(hits)} hits to {len(filtered_hits)} hits "
            f"(min_identity={min_identity}%, min_coverage={min_coverage}%)"
        )

        return filtered_hits


def validate_mmseqs_memory_limit(value: str) -> str:
    """Validate a MMseqs2-style memory limit string.

    Accepts a positive integer (no leading zeros) followed by a single unit
    suffix ``K``, ``M``, ``G``, or ``T`` (case-sensitive, no space).  Examples
    of valid values::

        64G   128G   400G   1T   512M

    Args:
        value: The memory limit string to validate.

    Returns:
        The validated string unchanged.

    Raises:
        ValueError: If the string is empty, missing units, has an invalid
            unit suffix, or is otherwise malformed.
    """
    if not value:
        raise ValueError("MMseqs2 memory limit must not be empty.")
    if not re.fullmatch(r'[1-9][0-9]*[KMGT]', value):
        raise ValueError(
            f"Invalid MMseqs2 memory limit: '{value}'. "
            "Expected a positive integer (no leading zeros) followed by a unit (K, M, G, or T), "
            "e.g. '400G' or '1T'."
        )
    return value


def run_mmseqs_search(
    query_fasta: Path,
    output_file: Path,
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
    database: str = "core_nt",
    db_path: Optional[Path] = None,
    num_threads: int = 1,
    sensitivity: float = 5.7,
    skip_existing: bool = True,
    search_type: int = 3,
    memory_limit: Optional[str] = None,
    timeout: int = 3600
) -> tuple[bool, list[BlastHit]]:
    """
    Convenience function to run an MMseqs2 search and return filtered hits.

    To keep results comparable with the BLAST workflow (which searches NCBI
    ``core_nt``), the MMseqs2 database should be built from the same underlying
    sequence collection as ``core_nt``.  There is no official pre-built MMseqs2
    ``core_nt`` database; users must create one from the relevant FASTA
    sequences (e.g. exported from BLAST ``core_nt`` with ``blastdbcmd``).
    Using a different nucleotide collection changes the search space and may
    produce classification differences unrelated to the choice of search engine.

    Args:
        query_fasta: Path to query FASTA file
        output_file: Path to output file
        min_identity: Minimum percent identity (default: 90.0)
        min_coverage: Minimum query coverage percentage (default: 80.0)
        database: Name of MMseqs2 database inside ``db_path`` (default: "core_nt")
        db_path: Path to the MMseqs2 database directory
        num_threads: Number of threads to use
        sensitivity: MMseqs2 sensitivity (1–7.5, default: 5.7)
        skip_existing: Skip search if output file already exists (default: True)
        search_type: MMseqs2 search type passed as ``--search-type`` to
            ``easy-search``. Commonly used values: 2 (translated),
            3 (nucleotide), 4 (translated nucleotide backtrace).
            Default is 3 (nucleotide). See MMseqs2 documentation for all
            valid values.
        memory_limit: Maximum RAM for the MMseqs2 prefilter/index step,
            passed as ``--split-memory-limit`` to ``easy-search``
            (e.g. ``"400G"``).  When ``None`` the flag is omitted.
        timeout: Maximum runtime for the MMseqs2 search in seconds
            (default: 3600).

    Returns:
        Tuple of (success: bool, filtered_hits: list[BlastHit])

    Raises:
        ValueError: If sensitivity is outside the valid range (1–7.5)
    """
    if not 1.0 <= sensitivity <= 7.5:
        raise ValueError(
            f"MMseqs2 sensitivity must be between 1.0 and 7.5, got {sensitivity}"
        )
    runner = MMseqs2Runner(db_path)

    if os.path.exists(output_file) and skip_existing:
        logger.info(
            f"The MMseqs2 output file {output_file} already exists. "
            "Skipping and using these results"
        )
    else:
        success = runner.run_easy_search(
            query_fasta=query_fasta,
            output_file=output_file,
            database=database,
            num_threads=num_threads,
            sensitivity=sensitivity,
            search_type=search_type,
            split_memory_limit=memory_limit,
            timeout=timeout
        )

        if not success:
            return False, []

    # Parse results
    try:
        hits = runner.parse_mmseqs_results(output_file)

        # Filter results
        filtered_hits = runner.filter_hits(
            hits,
            min_identity=min_identity,
            min_coverage=min_coverage
        )

        return True, filtered_hits

    except Exception as e:
        logger.error(f"Error processing MMseqs2 results: {e}")
        return False, []
