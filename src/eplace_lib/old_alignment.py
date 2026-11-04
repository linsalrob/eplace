"""
Sequence alignment and phylogenetic tree building module.

This module provides functionality for trimming sequences based on BLAST alignments,
aligning sequences using MAFFT, and building phylogenetic trees using IQTree.
"""

import os
import sys
import subprocess
import logging
from pathlib import Path
from typing import Optional, List, Dict, Tuple
from collections import defaultdict

from .blast_analysis import BlastHit, FastaReader
from .taxonomy import SequenceExtractor

# Configure module logger
logger = logging.getLogger(__name__)


class SequenceTrimmer:
    """
    Class for trimming sequences based on BLAST alignment coordinates.
    """
    
    @staticmethod
    def trim_sequence_by_coordinates(
        sequence: str,
        start: int,
        end: int
    ) -> str:
        """
        Trim a sequence to extract the region between start and end coordinates.
        
        BLAST coordinates are 1-indexed, so we need to adjust for Python's 0-indexing.
        
        Args:
            sequence: The full sequence string
            start: Start position (1-indexed, inclusive)
            end: End position (1-indexed, inclusive)
            
        Returns:
            Trimmed sequence string
        """
        # Convert 1-indexed BLAST coordinates to 0-indexed Python
        # Also handle reverse complement alignments (start > end)
        if start > end:
            # Reverse strand alignment
            python_start = end - 1
            python_end = start
        else:
            # Forward strand alignment
            python_start = start - 1
            python_end = end
        
        # Ensure coordinates are within bounds
        python_start = max(0, python_start)
        python_end = min(len(sequence), python_end)
        
        return sequence[python_start:python_end]
    
    @staticmethod
    def trim_sequences_from_blast_hits(
        fasta_path: Path,
        blast_hits: List[BlastHit],
        output_fasta: Path,
        query_id: str,
        taxonomic_rank: str
    ) -> bool:
        """
        Trim sequences in a FASTA file based on BLAST hit coordinates.
        
        This reads the representative sequences, trims them to the aligned regions,
        and writes them to a new FASTA file along with the query sequence.
        
        Args:
            fasta_path: Path to input FASTA file with full-length sequences
            blast_hits: List of BlastHit objects for this query
            output_fasta: Path to output FASTA file with trimmed sequences
            query_id: The query sequence ID to include in output
            taxonomic_rank: the taxonomic rank to use for taxonomic labels (e.g., "genus")

        Returns:
            True if successful, False otherwise
        """
        try:
            # Read all sequences from the input FASTA
            sequences = FastaReader.read_fasta(fasta_path)
            
            # Create a mapping of subject accession to blast hits for quick lookup
            # The FASTA file will have accessions (e.g., MZ387488.1) not full IDs
            hit_map = {hit.get_accession(): hit for hit in blast_hits}
            
            # Open output file
            with open(output_fasta, 'w') as out:
                # First, write the query sequence if it exists
                if query_id in sequences:
                    query_seq = sequences[query_id]
                    out.write(f">{query_id}\n")
                    # Write sequence in lines of 60 characters
                    for i in range(0, len(query_seq), 60):
                        out.write(query_seq[i:i+60] + "\n")
                    logger.info(f"Added query sequence {query_id} ({len(query_seq)} bp)")
                
                # Now process subject sequences
                for seq_id, sequence in sequences.items():
                    if seq_id == query_id:
                        continue  # Skip query, already written
                    
                    # Find the corresponding BLAST hit by accession
                    # The seq_id from FASTA should match the accession from the hit
                    hit = hit_map.get(seq_id)
                    
                    if hit is None:
                        # This sequence doesn't have a BLAST hit, skip it
                        logger.warning(f"No BLAST hit found for sequence {seq_id}, skipping")
                        continue
                    
                    # Trim the sequence based on subject coordinates
                    trimmed_seq = SequenceTrimmer.trim_sequence_by_coordinates(
                        sequence,
                        hit.subject_start,
                        hit.subject_end
                    )
                    
                    # Write trimmed sequence with taxonomic information in header
                    header = seq_id
                    if isinstance(hit.subject_taxonomy, dict) and taxonomic_rank in hit.subject_taxonomy:
                        header = f"{seq_id} {hit.subject_taxonomy[taxonomic_rank][1]}"
                    
                    out.write(f">{header}\n")
                    # Write sequence in lines of 60 characters
                    for i in range(0, len(trimmed_seq), 60):
                        out.write(trimmed_seq[i:i+60] + "\n")
                    
                    logger.info(
                        f"Trimmed {seq_id} from {len(sequence)} bp to {len(trimmed_seq)} bp "
                        f"(coords: {hit.subject_start}-{hit.subject_end})"
                    )
            
            return True
            
        except Exception as e:
            logger.error(f"Error trimming sequences: {e}")
            return False


class MAFFTAligner:
    """
    Class for running MAFFT sequence alignments.
    """
    
    @staticmethod
    def check_mafft_available() -> bool:
        """
        Check if MAFFT is available in the system.
        
        Returns:
            True if MAFFT is available, False otherwise
        """
        try:
            result = subprocess.run(
                ['mafft', '--version'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                logger.info(f"Using mafft to build the alignment ({result.stdout})")
            return result.returncode == 0
        except (subprocess.SubprocessError, FileNotFoundError):
            return False
    
    @staticmethod
    def align_sequences(
        input_fasta: Path,
        output_fasta: Path,
        auto_orient: bool = True,
        num_threads: int = 1,
        strategy: str = 'default'
    ) -> bool:
        """
        Align sequences using MAFFT.
        
        Args:
            input_fasta: Path to input FASTA file with sequences to align
            output_fasta: Path to output aligned FASTA file
            auto_orient: Use MAFFT's auto-orient feature (default: True)
            num_threads: Number of threads to use
            strategy: MAFFT alignment strategy (default: 'default')
                      Options: 'default', 'auto', 'retree2', 'fftns'
                      'auto': Let MAFFT choose the best strategy automatically
                      'retree2': Fast progressive method, good for large datasets
                      'fftns': Fastest method for very large datasets
            
        Returns:
            True if alignment was successful, False otherwise
        """
        if not MAFFTAligner.check_mafft_available():
            logger.error("MAFFT is not available. Please install MAFFT.")
            return False
        
        if not input_fasta.exists():
            logger.error(f"Input FASTA file not found: {input_fasta}")
            return False
        
        # Build MAFFT command
        cmd = ['mafft']
        
        # Add strategy-specific options
        if strategy == 'auto':
            cmd.append('--auto')
        elif strategy == 'retree2':
            cmd.append('--retree')
            cmd.append('2')
        elif strategy == 'fftns':
            cmd.append('--retree')
            cmd.append('1')
        # 'default' uses MAFFT's default strategy (no special flag)
        
        # Add auto-orient option
        if auto_orient:
            cmd.append('--adjustdirection')
        
        # Add threading
        cmd.extend(['--thread', str(num_threads)])
        
        # Add input file
        cmd.append(str(input_fasta))
        
        logger.info(f"Running MAFFT alignment: {' '.join(cmd)}")
        
        try:
            with open(output_fasta, 'w') as out:
                result = subprocess.run(
                    cmd,
                    stdout=out,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=3600  # 1 hour timeout
                )
            
            if result.returncode != 0:
                logger.error(f"MAFFT failed with error: {result.stderr}")
                return False
            
            logger.info(f"MAFFT alignment completed successfully. Output: {output_fasta}")
            return True
            
        except subprocess.TimeoutExpired:
            logger.error("MAFFT alignment timed out")
            return False
        except Exception as e:
            logger.error(f"Error running MAFFT: {e}")
            return False


class IQTreeBuilder:
    """
    Class for building phylogenetic trees using IQTree.
    """
    
    @staticmethod
    def check_iqtree_available() -> Tuple[bool, Optional[str]]:
        """
        Check if IQTree is available in the system.
        
        Returns:
            Tuple of (available: bool, command: str or None)
        """
        # Try 'iqtree', 'iqtree2', and 'iqtree3' commands
        for cmd in ['iqtree3', 'iqtree2', 'iqtree']:
            try:
                result = subprocess.run(
                    [cmd, '--version'],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                if result.returncode == 0:
                    logger.info(f"Using {cmd} to build the trees ({result.stdout})")
                    return True, cmd
            except (subprocess.SubprocessError, FileNotFoundError):
                continue
        
        return False, None
    
    @staticmethod
    def build_tree(
        alignment_fasta: Path,
        output_prefix: Path,
        model: str = "MFP",
        num_threads: int = None
    ) -> bool:
        """
        Build a phylogenetic tree using IQTree.
        
        Args:
            alignment_fasta: Path to aligned FASTA file
            output_prefix: Prefix for output files
            model: Substitution model (default: "MFP" for automatic ModelFinder Plus selection)
            num_threads: Number of threads to use (default: None, which uses AUTO)
            
        Returns:
            True if tree building was successful, False otherwise
        """
        available, iqtree_cmd = IQTreeBuilder.check_iqtree_available()
        if not available:
            logger.error("IQTree is not available. Please install IQTree or IQTree2.")
            return False
        
        if not alignment_fasta.exists():
            logger.error(f"Alignment file not found: {alignment_fasta}")
            return False
        
        # Build IQTree command
        cmd = [
            iqtree_cmd,
            '-s', str(alignment_fasta),
            '-pre', str(output_prefix),
            '-m', model,
            '-T', str(num_threads) if num_threads else "AUTO"
        ]
        
        logger.info(f"Running IQTree: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=14400  # 4 hour timeout (increased from 1 hour to handle larger datasets)
            )
            
            if result.returncode != 0:
                logger.error(f"IQTree failed with error: {result.stderr}")
                return False
            
            # Check if tree file was created
            tree_file = Path(str(output_prefix) + ".treefile")
            if not tree_file.exists():
                logger.error("IQTree did not produce a tree file")
                return False
            
            logger.info(f"IQTree completed successfully. Tree: {tree_file}")
            return True
            
        except subprocess.TimeoutExpired:
            logger.error("IQTree timed out")
            return False
        except Exception as e:
            logger.error(f"Error running IQTree: {e}")
            return False
    
    @staticmethod
    def build_tree_background(
        alignment_fasta: Path,
        output_prefix: Path,
        model: str = "MFP"
    ) -> Optional[Dict]:
        """
        Start building a phylogenetic tree using IQTree in the background.
        
        This method starts IQTree as a background process and returns immediately,
        allowing multiple trees to be built in parallel.
        
        Args:
            alignment_fasta: Path to aligned FASTA file
            output_prefix: Prefix for output files
            model: Substitution model (default: "MFP" for automatic ModelFinder Plus selection)
            
        Returns:
            Dictionary with process information if successful, None otherwise:
                - 'process': subprocess.Popen object
                - 'output_prefix': output prefix path
                - 'alignment_fasta': input alignment file path
                - 'tree_file': expected tree file path
        """
        available, iqtree_cmd = IQTreeBuilder.check_iqtree_available()
        if not available:
            logger.error("IQTree is not available. Please install IQTree or IQTree2.")
            return None
        
        if not alignment_fasta.exists():
            logger.error(f"Alignment file not found: {alignment_fasta}")
            return None
        
        # Build IQTree command
        cmd = [
            iqtree_cmd,
            '-s', str(alignment_fasta),
            '-pre', str(output_prefix),
            '-m', model,
            '-T', "AUTO"
        ]
        
        logger.info(f"Starting IQTree in background: {' '.join(cmd)}")
        
        try:
            # Start process in background
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            
            tree_file = Path(str(output_prefix) + ".treefile")
            
            return {
                'process': process,
                'output_prefix': output_prefix,
                'alignment_fasta': alignment_fasta,
                'tree_file': tree_file,
                'cmd': ' '.join(cmd)
            }
            
        except Exception as e:
            logger.error(f"Error starting IQTree: {e}")
            return None
    
    @staticmethod
    def wait_for_tree_jobs(
        jobs: List[Dict],
        timeout: int = 14400
    ) -> Dict[str, bool]:
        """
        Wait for multiple IQTree jobs to complete.
        
        This method polls all running processes and waits for them to complete.
        Since the processes are already running in parallel (started with Popen),
        this method just collects their results as they finish.
        
        Args:
            jobs: List of job dictionaries returned by build_tree_background()
            timeout: Maximum time to wait for each individual job in seconds (default: 14400 = 4 hours)
                     Increased because of mega tree created at the end!
            
        Returns:
            Dictionary mapping tree_file path to success status (True/False)
        """
        results = {}
        
        logger.info(f"Waiting for {len(jobs)} IQTree jobs to complete...")
        
        for job in jobs:
            process = job['process']
            output_prefix = job['output_prefix']
            tree_file = job['tree_file']
            cmd = job.get('cmd', 'IQTree')
            
            try:
                # Wait for process to complete
                # Note: This waits for THIS process, but other processes continue running in parallel
                stdout, stderr = process.communicate(timeout=timeout)
                
                if process.returncode != 0:
                    logger.error(f"IQTree failed for {output_prefix} with error: {stderr}")
                    results[str(tree_file)] = False
                    continue
                
                # Check if tree file was created
                if not tree_file.exists():
                    logger.error(f"IQTree did not produce a tree file for {output_prefix}")
                    results[str(tree_file)] = False
                    continue
                
                logger.info(f"IQTree completed successfully for {output_prefix}. Tree: {tree_file}")
                results[str(tree_file)] = True
                
            except subprocess.TimeoutExpired:
                logger.error(f"IQTree timed out for {output_prefix} after {timeout} seconds")
                process.kill()
                # Clean up the killed process
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    logger.error(f"Failed to terminate IQTree process for {output_prefix}")
                results[str(tree_file)] = False
            except Exception as e:
                logger.error(f"Error waiting for IQTree job {output_prefix}: {e}")
                # Ensure process is cleaned up
                if process.poll() is None:  # Process still running
                    process.kill()
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        pass
                results[str(tree_file)] = False
        
        successful = sum(1 for success in results.values() if success)
        logger.info(f"Completed {successful}/{len(jobs)} IQTree jobs successfully")
        
        return results
    
    @staticmethod
    def relabel_tree_with_taxonomy(
        tree_file: Path,
        blast_hits: List[BlastHit],
        output_tree: Path,
        taxonomic_rank: str,
    ) -> bool:
        """
        Relabel tree nodes with taxonomic names.
        
        This reads a Newick tree file and replaces sequence IDs with taxonomic names
        from the BLAST hits.
        
        Args:
            tree_file: Path to input tree file (Newick format)
            blast_hits: List of BlastHit objects with taxonomic information
            output_tree: Path to output tree file with relabeled nodes
            taxonomic_rank: the taxonomic rank to use for relabeling (e.g., "genus")
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Create mapping of sequence accession to taxonomic name
            # Trees will have accessions (e.g., MZ387488.1) not full IDs
            label_map = {}
            for hit in blast_hits:
                label = "unknown"
                if isinstance(hit.subject_taxonomy, dict) and taxonomic_rank in hit.subject_taxonomy:
                    label = hit.subject_taxonomy[taxonomic_rank][1]
                label_map[hit.subject_id] = label
                if label:
                    # Clean up the label for tree format (Newick format constraints)
                    # Replace spaces, colons, parentheses, commas, and semicolons
                    clean_label = (label.replace(' ', '_')
                                  .replace(':', '_')
                                  .replace('(', '_')
                                  .replace(')', '_')
                                  .replace(',', '_')
                                  .replace(';', '_'))
                    # Use accession for mapping since that's what appears in trees
                    accession = hit.get_accession()
                    label_map[accession] = clean_label

            # Read the tree file
            with open(tree_file, 'r') as f:
                tree_string = f.read()
            
            # Replace sequence IDs with taxonomic names
            for seq_id, tax_name in label_map.items():
                # Handle normal sequences (not reversed)
                tree_string = tree_string.replace(f"({seq_id}:", f"({tax_name}:")
                tree_string = tree_string.replace(f",{seq_id}:", f",{tax_name}:")
                tree_string = tree_string.replace(f" {seq_id}:", f" {tax_name}:")
                
                # Handle sequences with _R_ prefix (reversed by MAFFT)
                # MAFFT prepends _R_ to sequence IDs when it adjusts direction
                # We need to remove _R_ to find the correct ID, then append "_R" to the label
                # (using underscore to maintain Newick format compliance, representing " R")
                reversed_seq_id = f"_R_{seq_id}"
                reversed_label = f"{tax_name}_R"
                tree_string = tree_string.replace(f"({reversed_seq_id}:", f"({reversed_label}:")
                tree_string = tree_string.replace(f",{reversed_seq_id}:", f",{reversed_label}:")
                tree_string = tree_string.replace(f" {reversed_seq_id}:", f" {reversed_label}:")
            
            # Write the relabeled tree
            with open(output_tree, 'w') as f:
                f.write(tree_string)
            
            logger.info(f"Tree relabeled with {len(label_map)} taxonomic names")
            logger.info(f"Relabeled tree saved to: {output_tree}")
            return True
            
        except Exception as e:
            logger.error(f"Error relabeling tree: {e}")
            return False


def process_query_alignment_and_tree(
    query_id: str,
    query_dir: Path,
    blast_hits: List[BlastHit],
    query_fasta: Path,
    taxonomic_rank: str,
    num_threads: int = 1
) -> Dict[str, Optional[Path]]:
    """
    Complete pipeline for a single query: trim, align, and build tree.
    
    Args:
        query_id: Query sequence identifier
        query_dir: Directory containing query-specific files
        blast_hits: List of BlastHit objects for this query (with taxonomy info)
        query_fasta: Path to original query FASTA file
        taxonomic_rank: The taxonomic rank to use for relabeling the tree
        num_threads: Number of threads to use
        
    Returns:
        Dictionary with paths to generated files:
            - 'trimmed_fasta': Trimmed sequences
            - 'alignment': Aligned sequences
            - 'tree': Phylogenetic tree
            - 'labeled_tree': Tree with taxonomic labels
    """
    results = {
        'trimmed_fasta': None,
        'alignment': None,
        'tree': None,
        'labeled_tree': None
    }
    
    # File paths
    safe_query_id = query_id.replace('|', '_').replace('/', '_')
    representatives_fasta = query_dir / f"{safe_query_id}_representatives.fasta"
    trimmed_fasta = query_dir / f"{safe_query_id}_trimmed.fasta"
    alignment_fasta = query_dir / f"{safe_query_id}_aligned.fasta"
    tree_prefix = query_dir / f"{safe_query_id}_tree"
    tree_file = Path(str(tree_prefix) + ".treefile")
    labeled_tree = query_dir / f"{safe_query_id}_tree_labeled.treefile"
    
    # Step 1: Read query sequence and add it to the representatives file
    try:
        query_sequences = FastaReader.read_fasta(query_fasta)
        if query_id not in query_sequences:
            logger.error(f"Query {query_id} not found in {query_fasta}")
            return results
        
        # Read representatives and combine with query
        if not representatives_fasta.exists():
            logger.error(f"Representatives file not found: {representatives_fasta}")
            return results
        
        # Create combined FASTA with query + representatives
        combined_fasta = query_dir / f"{safe_query_id}_with_query.fasta"
        with open(combined_fasta, 'w') as out:
            # Write query first
            query_seq = query_sequences[query_id]
            out.write(f">{query_id}\n")
            for i in range(0, len(query_seq), 60):
                out.write(query_seq[i:i+60] + "\n")
            
            # Append representatives
            with open(representatives_fasta, 'r') as rep:
                out.write(rep.read())
        
        logger.info(f"Combined query with representatives: {combined_fasta}")
        
    except Exception as e:
        logger.error(f"Error preparing sequences: {e}")
        return results
    
    # Step 2: Trim sequences based on BLAST coordinates
    logger.info(f"Trimming sequences for {query_id}...")
    if SequenceTrimmer.trim_sequences_from_blast_hits(
        fasta_path=combined_fasta,
        blast_hits=blast_hits,
        output_fasta=trimmed_fasta,
        taxonomic_rank=taxonomic_rank,
        query_id=query_id
    ):
        results['trimmed_fasta'] = trimmed_fasta
        logger.info(f"Trimmed sequences saved to: {trimmed_fasta}")
    else:
        logger.error(f"Failed to trim sequences for {query_id}")
        return results
    
    # Step 3: Align sequences with MAFFT
    logger.info(f"Aligning sequences for {query_id}...")
    if MAFFTAligner.align_sequences(
        input_fasta=trimmed_fasta,
        output_fasta=alignment_fasta,
        auto_orient=True,
        num_threads=num_threads
    ):
        results['alignment'] = alignment_fasta
        logger.info(f"Alignment saved to: {alignment_fasta}")
    else:
        logger.error(f"Failed to align sequences for {query_id}")
        return results
    
    # Step 4: Build phylogenetic tree with IQTree
    logger.info(f"Building phylogenetic tree for {query_id}...")
    if IQTreeBuilder.build_tree(
        alignment_fasta=alignment_fasta,
        output_prefix=tree_prefix,
    ):
        results['tree'] = tree_file
        logger.info(f"Tree saved to: {tree_file}")
    else:
        logger.error(f"Failed to build tree for {query_id}")
        return results
    
    # Step 5: Relabel tree with taxonomic names
    logger.info(f"Relabeling tree with taxonomic names for {query_id}...")
    if IQTreeBuilder.relabel_tree_with_taxonomy(
        tree_file=tree_file,
        blast_hits=blast_hits,
        output_tree=labeled_tree,
        taxonomic_rank=taxonomic_rank
    ):
        results['labeled_tree'] = labeled_tree
        logger.info(f"Labeled tree saved to: {labeled_tree}")
    else:
        logger.warning(f"Failed to relabel tree for {query_id}, but unlabeled tree is available")
    
    return results


def process_query_alignment_and_tree_parallel(
    query_id: str,
    query_dir: Path,
    blast_hits: List[BlastHit],
    query_fasta: Path,
    taxonomic_rank: str,
    num_threads: int = 1,
    background_tree: bool = False
) -> Dict[str, Optional[Path]]:
    """
    Complete pipeline for a single query: trim, align, and optionally build tree in background.
    
    This is similar to process_query_alignment_and_tree, but with an option to start
    tree building in the background and return immediately without waiting for completion.
    
    Args:
        query_id: Query sequence identifier
        query_dir: Directory containing query-specific files
        blast_hits: List of BlastHit objects for this query (with taxonomy info)
        query_fasta: Path to original query FASTA file
        taxonomic_rank: The taxonomic rank to use for relabeling the tree
        num_threads: Number of threads to use
        background_tree: If True, start tree building in background and return immediately
        
    Returns:
        Dictionary with paths to generated files:
            - 'trimmed_fasta': Trimmed sequences
            - 'alignment': Aligned sequences
            - 'tree_job': Background job info if background_tree=True, None otherwise
            - 'tree_file': Expected tree file path
            - 'blast_hits': BLAST hits for later tree relabeling
            - 'taxonomic_rank': Taxonomic rank for later tree relabeling
    """
    results = {
        'trimmed_fasta': None,
        'alignment': None,
        'tree_job': None,
        'tree_file': None,
        'labeled_tree_path': None,
        'blast_hits': blast_hits,
        'taxonomic_rank': taxonomic_rank
    }
    
    # File paths
    safe_query_id = query_id.replace('|', '_').replace('/', '_')
    representatives_fasta = query_dir / f"{safe_query_id}_representatives.fasta"
    trimmed_fasta = query_dir / f"{safe_query_id}_trimmed.fasta"
    alignment_fasta = query_dir / f"{safe_query_id}_aligned.fasta"
    tree_prefix = query_dir / f"{safe_query_id}_tree"
    tree_file = Path(str(tree_prefix) + ".treefile")
    labeled_tree = query_dir / f"{safe_query_id}_tree_labeled.treefile"
    
    # Step 1: Read query sequence and add it to the representatives file
    try:
        query_sequences = FastaReader.read_fasta(query_fasta)
        if query_id not in query_sequences:
            logger.error(f"Query {query_id} not found in {query_fasta}")
            return results
        
        # Read representatives and combine with query
        if not representatives_fasta.exists():
            logger.error(f"Representatives file not found: {representatives_fasta}")
            return results
        
        # Create combined FASTA with query + representatives
        combined_fasta = query_dir / f"{safe_query_id}_with_query.fasta"
        with open(combined_fasta, 'w') as out:
            # Write query first
            query_seq = query_sequences[query_id]
            out.write(f">{query_id}\n")
            for i in range(0, len(query_seq), 60):
                out.write(query_seq[i:i+60] + "\n")
            
            # Append representatives
            with open(representatives_fasta, 'r') as rep:
                out.write(rep.read())
        
        logger.info(f"Combined query with representatives: {combined_fasta}")
        
    except Exception as e:
        logger.error(f"Error preparing sequences: {e}")
        return results
    
    # Step 2: Trim sequences based on BLAST coordinates
    logger.info(f"Trimming sequences for {query_id}...")
    if SequenceTrimmer.trim_sequences_from_blast_hits(
        fasta_path=combined_fasta,
        blast_hits=blast_hits,
        output_fasta=trimmed_fasta,
        taxonomic_rank=taxonomic_rank,
        query_id=query_id
    ):
        results['trimmed_fasta'] = trimmed_fasta
        logger.info(f"Trimmed sequences saved to: {trimmed_fasta}")
    else:
        logger.error(f"Failed to trim sequences for {query_id}")
        return results
    
    # Step 3: Align sequences with MAFFT
    logger.info(f"Aligning sequences for {query_id}...")
    if MAFFTAligner.align_sequences(
        input_fasta=trimmed_fasta,
        output_fasta=alignment_fasta,
        auto_orient=True,
        num_threads=num_threads
    ):
        results['alignment'] = alignment_fasta
        logger.info(f"Alignment saved to: {alignment_fasta}")
    else:
        logger.error(f"Failed to align sequences for {query_id}")
        return results
    
    # Step 4: Build phylogenetic tree with IQTree
    if background_tree:
        # Start tree building in background
        logger.info(f"Starting phylogenetic tree building in background for {query_id}...")
        tree_job = IQTreeBuilder.build_tree_background(
            alignment_fasta=alignment_fasta,
            output_prefix=tree_prefix,
        )
        if tree_job:
            results['tree_job'] = tree_job
            results['tree_file'] = tree_file
            results['labeled_tree_path'] = labeled_tree
            logger.info(f"Tree building started in background for {query_id}")
        else:
            logger.error(f"Failed to start tree building for {query_id}")
    else:
        # Build tree synchronously
        logger.info(f"Building phylogenetic tree for {query_id}...")
        if IQTreeBuilder.build_tree(
            alignment_fasta=alignment_fasta,
            output_prefix=tree_prefix,
        ):
            results['tree_file'] = tree_file
            logger.info(f"Tree saved to: {tree_file}")
            
            # Step 5: Relabel tree with taxonomic names
            logger.info(f"Relabeling tree with taxonomic names for {query_id}...")
            if IQTreeBuilder.relabel_tree_with_taxonomy(
                tree_file=tree_file,
                blast_hits=blast_hits,
                output_tree=labeled_tree,
                taxonomic_rank=taxonomic_rank
            ):
                results['labeled_tree_path'] = labeled_tree
                logger.info(f"Labeled tree saved to: {labeled_tree}")
            else:
                logger.warning(f"Failed to relabel tree for {query_id}, but unlabeled tree is available")
        else:
            logger.error(f"Failed to build tree for {query_id}")
    
    return results


def check_alignment_consistency(blast_hits: List[BlastHit], tolerance: int = 50) -> Dict[str, bool]:
    """
    Check if BLAST hits align to similar locations on reference sequences.
    
    For each reference sequence that appears in multiple hits, check if the alignment
    coordinates are consistent (within tolerance).
    
    Args:
        blast_hits: List of BlastHit objects to check
        tolerance: Maximum allowed difference in coordinates (default: 50 bp)
        
    Returns:
        Dictionary mapping subject_id to consistency status (True if consistent)
    """
    # Group hits by subject sequence
    hits_by_subject = defaultdict(list)
    for hit in blast_hits:
        hits_by_subject[hit.subject_id].append(hit)
    
    consistency_status = {}
    
    for subject_id, subject_hits in hits_by_subject.items():
        if len(subject_hits) == 1:
            # Only one hit, so it's consistent by definition
            consistency_status[subject_id] = True
            continue
        
        # Check if all hits have similar start and end coordinates
        starts = [hit.subject_start for hit in subject_hits]
        ends = [hit.subject_end for hit in subject_hits]
        
        start_range = max(starts) - min(starts)
        end_range = max(ends) - min(ends)
        
        is_consistent = start_range <= tolerance and end_range <= tolerance
        consistency_status[subject_id] = is_consistent
        
        if not is_consistent:
            logger.warning(
                f"Subject {subject_id} has inconsistent alignments: "
                f"start range={start_range}, end range={end_range}"
            )
        else:
            logger.info(
                f"Subject {subject_id} has consistent alignments across {len(subject_hits)} hits"
            )
    
    return consistency_status


def group_hits_by_group_rank(
    blast_hits: List[BlastHit],
    group_rank: str,
) -> Dict[str, Dict[str, List[BlastHit]]]:
    """
    Group BLAST hits by group_rank across all queries.
    
    Args:
        blast_hits: List of BlastHit objects with group taxonomy information
        
    Returns:
        Dictionary mapping group_rank_name (taxonomy name) to another dict mapping
        query_id to list of hits.
        Format: {group_rank_name: {query_id: [hits]}}
    """
    grouped = defaultdict(lambda: defaultdict(list))
    
    for hit in blast_hits:
        if hit.subject_taxonomy and group_rank in hit.subject_taxonomy:
            grouped[hit.subject_taxonomy[group_rank][1]][hit.query_id].append(hit)
        else:
            logger.warning(
                f"Hit {hit.subject_id} for query {hit.query_id} has no group taxonomy information"
            )
    
    logger.info(f"Grouped hits into {len(grouped)} taxonomic groups")
    for group_name, queries in grouped.items():
        logger.info(
            f"  Group {group_name}: {len(queries)} queries, "
            f"{sum(len(hits) for hits in queries.values())} total hits"
        )
    
    return dict(grouped)


def create_grouped_fasta_with_queries(
    group_tid: str,
    group_name: str,
    query_hits_map: Dict[str, List[BlastHit]],
    labeling_rank: str,
    query_fasta: Path,
    output_fasta: Path,
    database: str = "core_nt",
    blastdb_path: Optional[Path] = None
) -> bool:
    """
    Create a FASTA file for a taxonomic group containing all queries and unique references.
    
    Args:
        group_tid: Taxonomy ID of the group
        group_name: Name of the taxonomic group
        query_hits_map: Dictionary mapping query_id to list of BlastHit objects
        labeling_rank: Taxonomic rank to use for labeling (e.g., "genus")
        query_fasta: Path to original query FASTA file
        output_fasta: Path to output grouped FASTA file
        database: Name of BLAST database
        blastdb_path: Path to BLAST database directory
        
    Returns:
        True if successful, False otherwise
    """
    logger.info(f"Creating grouped FASTA for {group_name} ({group_tid})")
    
    # Read all query sequences
    try:
        query_sequences = FastaReader.read_fasta(query_fasta)
    except Exception as e:
        logger.error(f"Error reading query FASTA: {e}")
        return False
    
    # Collect unique reference sequences (by labeling rank) so we only get one example
    # For each unique reference, keep the hit with the best bit score
    unique_references = {}
    unique_labels = {}
    for query_id, hits in query_hits_map.items():
        for hit in hits:
            label = hit.subject_id
            if isinstance(hit.subject_taxonomy, dict) and labeling_rank in hit.subject_taxonomy:
                label = hit.subject_taxonomy[labeling_rank][1]
            if label not in unique_labels:
                unique_labels[label] = hit
                unique_references[hit.subject_id] = hit
            else:
                # Keep the hit with better bit score
                if hit.bit_score > unique_labels[label].bit_score:
                    # Remove the old reference for this label, if present
                    old_hit = unique_labels[label]
                    old_subject_id = old_hit.subject_id
                    if old_subject_id in unique_references:
                        del unique_references[old_subject_id]
                    # Add the new, better-scoring hit
                    unique_references[hit.subject_id] = hit
                    unique_labels[label] = hit

    logger.info(f"Found {len(unique_references)} unique reference sequences")
    unique_label_keys = list(unique_labels.keys())
    logger.info(
        "Unique labels: count=%d, example_labels=%s",
        len(unique_labels),
        unique_label_keys[:10],
    )
    
    # Extract reference sequences
    seq_extractor = SequenceExtractor(blastdb_path)
    temp_ref_fasta = output_fasta.parent / f"{output_fasta.stem}_temp_refs.fasta"
    
    try:
        success = seq_extractor.extract_sequences(
            sequence_ids=list(unique_references.keys()),
            output_fasta=temp_ref_fasta,
            database=database
        )
        
        if not success:
            logger.error("Failed to extract reference sequences")
            sys.exit(1)
        
        # Read extracted references
        ref_sequences = FastaReader.read_fasta(temp_ref_fasta)
        
        # Write combined FASTA file
        with open(output_fasta, 'w') as out:
            # Write all query sequences first
            for query_id in query_hits_map.keys():
                if query_id in query_sequences:
                    query_seq = query_sequences[query_id]
                    out.write(f">{query_id}\n")
                    for i in range(0, len(query_seq), 60):
                        out.write(query_seq[i:i+60] + "\n")
                    logger.info(f"Added query {query_id} ({len(query_seq)} bp)")
                else:
                    logger.warning(f"Query {query_id} not found in query FASTA file, skipping")
            
            # Write reference sequences with taxonomic labels
            for subject_id, hit in unique_references.items():
                # Get accession for lookup
                accession = hit.get_accession()
                if accession in ref_sequences:
                    ref_seq = ref_sequences[accession]
                    header = accession
                    if isinstance(hit.subject_taxonomy, dict) and labeling_rank in hit.subject_taxonomy:
                        header = f"{accession} {hit.subject_taxonomy[labeling_rank][1]}"
                    
                    out.write(f">{header}\n")
                    for i in range(0, len(ref_seq), 60):
                        out.write(ref_seq[i:i+60] + "\n")
                    logger.info(f"Added reference {accession} ({len(ref_seq)} bp)")
                else:
                    logger.warning(f"Reference {accession} not found in extracted sequences, skipping")
        
        logger.info(f"Created grouped FASTA file: {output_fasta}")
        return True
        
    except Exception as e:
        logger.error(f"Error creating grouped FASTA: {e}")
        return False
    finally:
        # Clean up temporary file
        if temp_ref_fasta.exists():
            temp_ref_fasta.unlink()


def trim_grouped_sequences(
    input_fasta: Path,
    blast_hits: List[BlastHit],
    output_fasta: Path,
    query_ids: List[str]
) -> bool:
    """
    Trim sequences in a grouped FASTA file based on BLAST hit coordinates.
    
    This is similar to trim_sequences_from_blast_hits but handles multiple queries.
    
    Args:
        input_fasta: Path to input FASTA file with full-length sequences
        blast_hits: List of BlastHit objects for all queries in the group
        output_fasta: Path to output FASTA file with trimmed sequences
        query_ids: List of query sequence IDs to include (untrimmed)
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Read all sequences from the input FASTA
        sequences = FastaReader.read_fasta(input_fasta)
        
        # Create a mapping of subject accession to blast hits
        hit_map = {}
        for hit in blast_hits:
            accession = hit.get_accession()
            if accession not in hit_map:
                hit_map[accession] = []
            hit_map[accession].append(hit)
        
        # For sequences with multiple hits, use the one with the best bit score
        # This is consistent with the deduplication logic in create_grouped_fasta_with_queries
        best_hits = {}
        for accession, hits in hit_map.items():
            if len(hits) == 1:
                best_hits[accession] = hits[0]
            else:
                # Use the hit with the best bit score for consistency
                best_hits[accession] = max(hits, key=lambda h: h.bit_score)
        
        # Open output file
        with open(output_fasta, 'w') as out:
            # First, write all query sequences (untrimmed)
            for query_id in query_ids:
                if query_id in sequences:
                    query_seq = sequences[query_id]
                    out.write(f">{query_id}\n")
                    for i in range(0, len(query_seq), 60):
                        out.write(query_seq[i:i+60] + "\n")
                    logger.info(f"Added query sequence {query_id} ({len(query_seq)} bp)")
            
            # Now process subject sequences (trimmed)
            for seq_id, sequence in sequences.items():
                if seq_id in query_ids:
                    continue  # Skip queries, already written
                
                # Extract just the accession from the header (might have taxonomy info)
                accession = seq_id.split()[0]
                hit = best_hits.get(accession)
                
                if hit is None:
                    logger.error(f"No BLAST hit found for sequence {accession}, data consistency issue")
                    continue
                
                # Trim the sequence based on subject coordinates
                trimmed_seq = SequenceTrimmer.trim_sequence_by_coordinates(
                    sequence,
                    hit.subject_start,
                    hit.subject_end
                )
                
                # Write trimmed sequence
                out.write(f">{seq_id}\n")
                for i in range(0, len(trimmed_seq), 60):
                    out.write(trimmed_seq[i:i+60] + "\n")
                
                logger.info(
                    f"Trimmed {accession} from {len(sequence)} bp to {len(trimmed_seq)} bp "
                    f"(coords: {hit.subject_start}-{hit.subject_end})"
                )
        
        return True
        
    except Exception as e:
        logger.error(f"Error trimming grouped sequences: {e}")
        return False


def process_grouped_alignment_and_tree(
    group_name: str,
    group_dir: Path,
    taxonomic_rank: str,
    blast_hits: List[BlastHit],
    query_ids: List[str],
    num_threads: int = 1
) -> Dict[str, Optional[Path]]:
    """
    Complete pipeline for a taxonomic group: trim, align, and build tree.
    
    Args:
        group_name: The name of the group, used for file naming
        group_dir: Directory containing group-specific files
        taxonomic_rank: Taxonomic rank to use for labeling the tree
        blast_hits: List of BlastHit objects for all queries in the group
        query_ids: List of query sequence IDs in this group
        num_threads: Number of threads to use
        
    Returns:
        Dictionary with paths to generated files:
            - 'combined_fasta': Combined sequences (queries + references)
            - 'trimmed_fasta': Trimmed sequences
            - 'alignment': Aligned sequences
            - 'tree': Phylogenetic tree
            - 'labeled_tree': Tree with taxonomic labels
    """
    results = {
        'combined_fasta': None,
        'trimmed_fasta': None,
        'alignment': None,
        'tree': None,
        'labeled_tree': None
    }

    # File paths
    safe_group_name = group_name.replace(' ', '_').replace('/', '_').replace('|', '_')
    combined_fasta = group_dir / f"{safe_group_name}_combined.fasta"
    trimmed_fasta = group_dir / f"{safe_group_name}_trimmed.fasta"
    alignment_fasta = group_dir / f"{safe_group_name}_aligned.fasta"
    tree_prefix = group_dir / f"{safe_group_name}_tree"
    tree_file = Path(str(tree_prefix) + ".treefile")
    labeled_tree = group_dir / f"{safe_group_name}_tree_labeled.treefile"
    
    # Check if combined FASTA exists
    if not combined_fasta.exists():
        logger.error(f"Combined FASTA file not found: {combined_fasta}")
        return results
    
    results['combined_fasta'] = combined_fasta
    
    # Step 1: Trim sequences based on BLAST coordinates
    logger.info(f"Trimming sequences for group {group_name}...")
    if trim_grouped_sequences(
        input_fasta=combined_fasta,
        blast_hits=blast_hits,
        output_fasta=trimmed_fasta,
        query_ids=query_ids
    ):
        results['trimmed_fasta'] = trimmed_fasta
        logger.info(f"Trimmed sequences saved to: {trimmed_fasta}")
    else:
        logger.error(f"Failed to trim sequences for group {group_name}")
        return results
    
    # Step 2: Align sequences with MAFFT
    logger.info(f"Aligning sequences for group {group_name}...")
    if MAFFTAligner.align_sequences(
        input_fasta=trimmed_fasta,
        output_fasta=alignment_fasta,
        auto_orient=True,
        num_threads=num_threads
    ):
        results['alignment'] = alignment_fasta
        logger.info(f"Alignment saved to: {alignment_fasta}")
    else:
        logger.error(f"Failed to align sequences for group {group_name}")
        return results
    
    # Step 3: Build phylogenetic tree with IQTree
    logger.info(f"Building phylogenetic tree for group {group_name}...")
    if IQTreeBuilder.build_tree(
        alignment_fasta=alignment_fasta,
        output_prefix=tree_prefix,
    ):
        results['tree'] = tree_file
        logger.info(f"Tree saved to: {tree_file}")
    else:
        logger.error(f"Failed to build tree for group {group_name}")
        return results
    
    # Step 4: Relabel tree with taxonomic names
    logger.info(f"Relabeling tree with taxonomic names for group {group_name}...")
    if IQTreeBuilder.relabel_tree_with_taxonomy(
        tree_file=tree_file,
        blast_hits=blast_hits,
        output_tree=labeled_tree,
        taxonomic_rank=taxonomic_rank
    ):
        results['labeled_tree'] = labeled_tree
        logger.info(f"Labeled tree saved to: {labeled_tree}")
    else:
        logger.warning(f"Failed to relabel tree for group {group_name}, but unlabeled tree is available")
    
    return results


def process_grouped_alignment_and_tree_parallel(
    group_name: str,
    group_dir: Path,
    taxonomic_rank: str,
    blast_hits: List[BlastHit],
    query_ids: List[str],
    num_threads: int = 1,
    background_tree: bool = False
) -> Dict[str, Optional[Path]]:
    """
    Complete pipeline for a taxonomic group: trim, align, and optionally build tree in background.
    
    This is similar to process_grouped_alignment_and_tree, but with an option to start
    tree building in the background and return immediately without waiting for completion.
    
    Args:
        group_name: The name of the group, used for file naming
        group_dir: Directory containing group-specific files
        taxonomic_rank: Taxonomic rank to use for labeling the tree
        blast_hits: List of BlastHit objects for all queries in the group
        query_ids: List of query sequence IDs in this group
        num_threads: Number of threads to use
        background_tree: If True, start tree building in background and return immediately
        
    Returns:
        Dictionary with paths to generated files:
            - 'combined_fasta': Combined sequences (queries + references)
            - 'trimmed_fasta': Trimmed sequences
            - 'alignment': Aligned sequences
            - 'tree_job': Background job info if background_tree=True, None otherwise
            - 'tree_file': Expected tree file path
            - 'blast_hits': BLAST hits for later tree relabeling
            - 'taxonomic_rank': Taxonomic rank for later tree relabeling
    """
    results = {
        'combined_fasta': None,
        'trimmed_fasta': None,
        'alignment': None,
        'tree_job': None,
        'tree_file': None,
        'labeled_tree_path': None,
        'blast_hits': blast_hits,
        'taxonomic_rank': taxonomic_rank
    }

    # File paths
    safe_group_name = group_name.replace(' ', '_').replace('/', '_').replace('|', '_')
    combined_fasta = group_dir / f"{safe_group_name}_combined.fasta"
    trimmed_fasta = group_dir / f"{safe_group_name}_trimmed.fasta"
    alignment_fasta = group_dir / f"{safe_group_name}_aligned.fasta"
    tree_prefix = group_dir / f"{safe_group_name}_tree"
    tree_file = Path(str(tree_prefix) + ".treefile")
    labeled_tree = group_dir / f"{safe_group_name}_tree_labeled.treefile"
    
    # Check if combined FASTA exists
    if not combined_fasta.exists():
        logger.error(f"Combined FASTA file not found: {combined_fasta}")
        return results
    
    results['combined_fasta'] = combined_fasta
    
    # Step 1: Trim sequences based on BLAST coordinates
    logger.info(f"Trimming sequences for group {group_name}...")
    if trim_grouped_sequences(
        input_fasta=combined_fasta,
        blast_hits=blast_hits,
        output_fasta=trimmed_fasta,
        query_ids=query_ids
    ):
        results['trimmed_fasta'] = trimmed_fasta
        logger.info(f"Trimmed sequences saved to: {trimmed_fasta}")
    else:
        logger.error(f"Failed to trim sequences for group {group_name}")
        return results
    
    # Step 2: Align sequences with MAFFT
    logger.info(f"Aligning sequences for group {group_name}...")
    if MAFFTAligner.align_sequences(
        input_fasta=trimmed_fasta,
        output_fasta=alignment_fasta,
        auto_orient=True,
        num_threads=num_threads
    ):
        results['alignment'] = alignment_fasta
        logger.info(f"Alignment saved to: {alignment_fasta}")
    else:
        logger.error(f"Failed to align sequences for group {group_name}")
        return results
    
    # Step 3: Build phylogenetic tree with IQTree
    if background_tree:
        # Start tree building in background
        logger.info(f"Starting phylogenetic tree building in background for group {group_name}...")
        tree_job = IQTreeBuilder.build_tree_background(
            alignment_fasta=alignment_fasta,
            output_prefix=tree_prefix,
        )
        if tree_job:
            results['tree_job'] = tree_job
            results['tree_file'] = tree_file
            results['labeled_tree_path'] = labeled_tree
            logger.info(f"Tree building started in background for {group_name}")
        else:
            logger.error(f"Failed to start tree building for group {group_name}")
    else:
        # Build tree synchronously
        logger.info(f"Building phylogenetic tree for group {group_name}...")
        if IQTreeBuilder.build_tree(
            alignment_fasta=alignment_fasta,
            output_prefix=tree_prefix,
        ):
            results['tree_file'] = tree_file
            logger.info(f"Tree saved to: {tree_file}")
            
            # Step 4: Relabel tree with taxonomic names
            logger.info(f"Relabeling tree with taxonomic names for group {group_name}...")
            if IQTreeBuilder.relabel_tree_with_taxonomy(
                tree_file=tree_file,
                blast_hits=blast_hits,
                output_tree=labeled_tree,
                taxonomic_rank=taxonomic_rank
            ):
                results['labeled_tree_path'] = labeled_tree
                logger.info(f"Labeled tree saved to: {labeled_tree}")
            else:
                logger.warning(f"Failed to relabel tree for group {group_name}, but unlabeled tree is available")
        else:
            logger.error(f"Failed to build tree for group {group_name}")
    
    return results


def concatenate_all_groups_and_build_tree(
    output_dir: Path,
    query_fasta: Path,
    classification_file: Path,
    blast_hits: List[BlastHit],
    combined_tree_label_rank: str = "genus",
    num_threads: int = 1,
    alignment_strategy: str = "auto"
) -> Dict[str, Optional[Path]]:
    """
    Concatenate all group _trimmed.fasta files, add queries with 0 blast hits,
    build a final alignment and tree.
    
    This function:
    1. Finds all *_trimmed.fasta files in group directories
    2. Reads the classification file to identify queries with 0 blast hits
    3. Concatenates all sequences into a single file
    4. Uses MAFFT to build an alignment (with optimal parameters for many sequences)
    5. Uses IQTree to build a phylogenetic tree
    6. Relabels tree nodes with taxonomic names
    
    Args:
        output_dir: Output directory containing group subdirectories
        query_fasta: Original query FASTA file
        classification_file: Path to classifications.tsv file
        blast_hits: List of all BlastHit objects with taxonomy information
        combined_tree_label_rank: Taxonomic rank for tree labeling (default: genus)
        num_threads: Number of threads for alignment and tree building (default: 1)
        alignment_strategy: MAFFT alignment strategy (default: 'auto')
                           Options: 'default', 'auto', 'retree2', 'fftns'
        
    Returns:
        Dictionary with paths to generated files:
            - 'combined_fasta': Combined sequences from all groups + zero-hit queries
            - 'alignment': Aligned sequences
            - 'tree': Phylogenetic tree
            - 'labeled_tree': Tree with taxonomic labels
    """
    results = {
        'combined_fasta': None,
        'alignment': None,
        'tree': None,
        'labeled_tree': None
    }
    
    # Define total steps for consistent logging
    TOTAL_STEPS = 5
    
    logger.info("\n" + "=" * 60)
    logger.info("Building combined tree from all groups")
    logger.info("=" * 60)
    
    # Output file paths
    combined_fasta = output_dir / "all_groups_combined.fasta"
    alignment_fasta = output_dir / "all_groups_aligned.fasta"
    tree_prefix = output_dir / "all_groups_tree"
    tree_file = Path(str(tree_prefix) + ".treefile")
    labeled_tree = output_dir / "all_groups_tree_labeled.treefile"
    
    try:
        # Step 1: Find all group directories and their _trimmed.fasta files
        logger.info(f"\n[Step 1/{TOTAL_STEPS}] Finding all group trimmed FASTA files...")
        trimmed_files = []
        
        # Look for directories in output_dir
        for item in output_dir.iterdir():
            if item.is_dir():
                # Look for *_trimmed.fasta files in this directory
                for fasta_file in item.glob("*_trimmed.fasta"):
                    trimmed_files.append(fasta_file)
                    logger.info(f"  Found: {fasta_file}")
        
        if not trimmed_files:
            logger.error("No _trimmed.fasta files found in group directories")
            return results
        
        logger.info(f"Found {len(trimmed_files)} trimmed FASTA files")
        
        # Step 2: Read classification file to find queries with 0 blast hits
        logger.info(f"\n[Step 2/{TOTAL_STEPS}] Identifying queries with 0 blast hits...")
        zero_hit_queries = []
        
        if classification_file.exists():
            with open(classification_file, 'r') as f:
                # Skip header
                next(f, None)
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        query_id = parts[0]
                        try:
                            blast_hits_count = int(parts[1])
                            if blast_hits_count == 0:
                                zero_hit_queries.append(query_id)
                                logger.info(f"  Query with 0 hits: {query_id}")
                        except ValueError:
                            # Skip lines where blast_hits_count is not a valid integer (e.g., 'N/A')
                            logger.debug(f"Skipping query {query_id} with non-numeric hit count: {parts[1]}")
                            continue
        else:
            logger.warning(f"Classification file not found: {classification_file}")
        
        if zero_hit_queries:
            logger.info(f"Found {len(zero_hit_queries)} queries with 0 blast hits")
        else:
            logger.info("No queries with 0 blast hits")
        
        # Step 3: Concatenate all sequences
        logger.info(f"\n[Step 3/{TOTAL_STEPS}] Concatenating all sequences...")
        
        # Read original query sequences
        query_sequences = FastaReader.read_fasta(query_fasta)
        
        # Track which sequences we've already written to avoid duplicates
        written_sequences = set()
        
        with open(combined_fasta, 'w') as out:
            # First, write all sequences from trimmed files
            for trimmed_file in sorted(trimmed_files):
                logger.info(f"  Reading: {trimmed_file}")
                sequences = FastaReader.read_fasta(trimmed_file)
                for seq_id, sequence in sequences.items():
                    if seq_id not in written_sequences:
                        out.write(f">{seq_id}\n")
                        for i in range(0, len(sequence), 60):
                            out.write(sequence[i:i+60] + "\n")
                        written_sequences.add(seq_id)
            
            # Then, add queries with 0 blast hits
            for query_id in zero_hit_queries:
                if query_id in query_sequences and query_id not in written_sequences:
                    query_seq = query_sequences[query_id]
                    out.write(f">{query_id}\n")
                    for i in range(0, len(query_seq), 60):
                        out.write(query_seq[i:i+60] + "\n")
                    written_sequences.add(query_id)
                    logger.info(f"  Added zero-hit query: {query_id}")
        
        results['combined_fasta'] = combined_fasta
        logger.info(f"Combined {len(written_sequences)} sequences into: {combined_fasta}")
        
        # Step 4: Align sequences with MAFFT (using optimal parameters for many sequences)
        logger.info(f"\n[Step 4/{TOTAL_STEPS}] Aligning sequences with MAFFT...")
        logger.info(f"Using '{alignment_strategy}' strategy for alignment")
        
        if MAFFTAligner.align_sequences(
            input_fasta=combined_fasta,
            output_fasta=alignment_fasta,
            auto_orient=True,
            num_threads=num_threads,
            strategy=alignment_strategy
        ):
            results['alignment'] = alignment_fasta
            logger.info(f"Alignment saved to: {alignment_fasta}")
        else:
            logger.error("Failed to align combined sequences")
            return results
        
        # Step 5: Build phylogenetic tree with IQTree
        logger.info(f"\n[Step 5/{TOTAL_STEPS}] Building phylogenetic tree...")
        if IQTreeBuilder.build_tree(
            alignment_fasta=alignment_fasta,
            output_prefix=tree_prefix,
            num_threads=num_threads
        ):
            results['tree'] = tree_file
            logger.info(f"Tree saved to: {tree_file}")
            
            # Relabel tree with taxonomic names
            logger.info("Relabeling tree with taxonomic names...")
            if IQTreeBuilder.relabel_tree_with_taxonomy(
                tree_file=tree_file,
                blast_hits=blast_hits,
                output_tree=labeled_tree,
                taxonomic_rank=combined_tree_label_rank
            ):
                results['labeled_tree'] = labeled_tree
                logger.info(f"Labeled tree saved to: {labeled_tree}")
            else:
                logger.warning("Failed to relabel tree, but unlabeled tree is available")
        else:
            logger.error("Failed to build phylogenetic tree")
            return results
        
        logger.info("\n" + "=" * 60)
        logger.info("Combined tree building completed successfully!")
        logger.info("=" * 60)
        
    except Exception as e:
        logger.error(f"Error in concatenate_all_groups_and_build_tree: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return results
    
    return results


class SimpleNewickNode:
    """
    Simple Newick tree node representation for finding nearest neighbors.
    """
    def __init__(self, name: str = "", distance: float = 0.0):
        self.name = name
        self.distance = distance
        self.children: List['SimpleNewickNode'] = []
        self.parent: Optional['SimpleNewickNode'] = None
    
    def is_leaf(self) -> bool:
        """Check if this node is a leaf."""
        return len(self.children) == 0
    
    def get_leaves(self) -> List['SimpleNewickNode']:
        """Get all leaf nodes under this node."""
        if self.is_leaf():
            return [self]
        leaves = []
        for child in self.children:
            leaves.extend(child.get_leaves())
        return leaves


def parse_simple_newick(newick_str: str) -> Optional[SimpleNewickNode]:
    """
    Parse a simple Newick tree string into a tree structure.
    
    This is a lightweight parser that handles basic Newick format with branch lengths.
    Format: ((A:0.1,B:0.2):0.3,C:0.4);
    
    Args:
        newick_str: Newick format tree string
        
    Returns:
        Root node of the parsed tree, or None if parsing fails
    """
    try:
        # Remove trailing semicolon and whitespace
        newick_str = newick_str.strip().rstrip(';').strip()
        
        # Stack to track parent nodes
        stack = []
        current_node = None
        
        i = 0
        token = ""
        
        while i < len(newick_str):
            char = newick_str[i]
            
            if char == '(':
                # Start a new internal node
                new_node = SimpleNewickNode()
                
                if current_node is not None:
                    # Add as child to current node
                    current_node.children.append(new_node)
                    new_node.parent = current_node
                
                stack.append(new_node)
                current_node = new_node
                token = ""
                
            elif char == ',':
                # Process accumulated token as a leaf
                if token:
                    name, distance = _parse_node_info(token)
                    leaf = SimpleNewickNode(name=name, distance=distance)
                    leaf.parent = current_node
                    current_node.children.append(leaf)
                    token = ""
                
            elif char == ')':
                # Process any remaining token
                if token:
                    name, distance = _parse_node_info(token)
                    leaf = SimpleNewickNode(name=name, distance=distance)
                    leaf.parent = current_node
                    current_node.children.append(leaf)
                    token = ""
                
                # Look for node label and distance after ')'
                i += 1
                label_token = ""
                while i < len(newick_str) and newick_str[i] not in '(),;':
                    label_token += newick_str[i]
                    i += 1
                
                if label_token and current_node:
                    name, distance = _parse_node_info(label_token)
                    if name:
                        current_node.name = name
                    current_node.distance = distance
                
                # Move up to parent
                if len(stack) > 1:
                    stack.pop()
                    current_node = stack[-1]
                
                i -= 1  # Adjust because we'll increment at the end
                
            else:
                # Accumulate characters for node name/distance
                token += char
            
            i += 1
        
        # Process any final token
        if token and current_node:
            name, distance = _parse_node_info(token)
            if name:
                current_node.name = name
            current_node.distance = distance
        
        # Return the root (first node in stack)
        return stack[0] if stack else None
        
    except Exception as e:
        logger.error(f"Error parsing Newick tree: {e}")
        return None


def _parse_node_info(token: str) -> Tuple[str, float]:
    """
    Parse node name and distance from a token like 'name:0.123' or just 'name' or ':0.123'.
    
    Args:
        token: Token string to parse
        
    Returns:
        Tuple of (name, distance)
    """
    token = token.strip()
    if ':' in token:
        parts = token.split(':', 1)
        name = parts[0].strip()
        try:
            distance = float(parts[1].strip())
        except (ValueError, IndexError):
            distance = 0.0
        return name, distance
    else:
        return token, 0.0


def find_nearest_neighbor_in_tree(
    tree_file: Path,
    query_id: str
) -> Optional[str]:
    """
    Find the nearest neighbor (closest leaf) to a query sequence in a phylogenetic tree.
    
    This function parses the Newick tree and finds the leaf node that is phylogenetically
    closest to the query sequence based on tree topology and branch lengths.
    
    Args:
        tree_file: Path to the Newick tree file (.treefile)
        query_id: Query sequence identifier to find neighbors for
        
    Returns:
        Name of the nearest neighbor leaf node, or None if not found or error
    """
    try:
        if not tree_file.exists():
            logger.warning(f"Tree file not found: {tree_file}")
            return None
        
        # Read the tree file
        tree_content = tree_file.read_text().strip()
        if not tree_content:
            logger.warning(f"Empty tree file: {tree_file}")
            return None
        
        # Parse the tree
        root = parse_simple_newick(tree_content)
        if not root:
            logger.warning(f"Failed to parse tree file: {tree_file}")
            return None
        
        # Find the query node
        query_node = _find_node_by_name(root, query_id)
        if not query_node:
            logger.warning(f"Query {query_id} not found in tree {tree_file}")
            return None
        
        # Get all leaves in the tree
        all_leaves = root.get_leaves()
        
        # Find the closest leaf (excluding the query itself)
        min_distance = float('inf')
        nearest_neighbor = None
        
        for leaf in all_leaves:
            # Skip the query itself
            if leaf.name == query_id:
                continue
            
            # Calculate distance between query and this leaf
            distance = _calculate_tree_distance(query_node, leaf)
            
            if distance < min_distance:
                min_distance = distance
                nearest_neighbor = leaf.name
        
        return nearest_neighbor
        
    except Exception as e:
        logger.error(f"Error finding nearest neighbor in tree {tree_file}: {e}")
        return None


def _find_node_by_name(node: SimpleNewickNode, name: str) -> Optional[SimpleNewickNode]:
    """
    Find a node with a given name in the tree.
    
    Args:
        node: Node to start search from
        name: Name to search for
        
    Returns:
        Node with matching name, or None if not found
    """
    if node.name == name:
        return node
    
    for child in node.children:
        result = _find_node_by_name(child, name)
        if result:
            return result
    
    return None


def _calculate_tree_distance(node1: SimpleNewickNode, node2: SimpleNewickNode) -> float:
    """
    Calculate the phylogenetic distance between two nodes in the tree.
    
    This finds the path from node1 to their most recent common ancestor (MRCA),
    then from MRCA to node2, summing the branch lengths.
    
    Args:
        node1: First node
        node2: Second node
        
    Returns:
        Total distance between the nodes
    """
    # Get paths from root to each node
    path1 = _get_path_to_root(node1)
    path2 = _get_path_to_root(node2)
    
    # Find the most recent common ancestor (MRCA)
    # The paths are from node to root, so reverse them (without modifying originals)
    path1_reversed = path1[::-1]
    path2_reversed = path2[::-1]
    
    # Find where paths diverge
    mrca_index = 0
    for i in range(min(len(path1_reversed), len(path2_reversed))):
        if path1_reversed[i] is path2_reversed[i]:
            mrca_index = i
        else:
            break
    
    # Calculate distance: sum from node1 to MRCA + sum from MRCA to node2
    distance = 0.0
    
    # Distance from node1 to MRCA (excluding MRCA itself to avoid double counting)
    for i in range(len(path1) - mrca_index - 1):
        distance += path1[i].distance
    
    # Distance from node2 to MRCA (excluding MRCA itself)
    for i in range(len(path2) - mrca_index - 1):
        distance += path2[i].distance
    
    return distance


def _get_path_to_root(node: SimpleNewickNode) -> List[SimpleNewickNode]:
    """
    Get the path from a node to the root.
    
    Args:
        node: Starting node
        
    Returns:
        List of nodes from the starting node to root
    """
    path = []
    current = node
    while current:
        path.append(current)
        current = current.parent
    return path
