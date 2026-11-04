"""
Tests for taxonomy extraction functionality.

Note: in this test suite we have "taxonomy IDs" at various locations. These are completely fictional and tend to be
1495400 or so. They are designed so that if the suite has an error it is much easier to locate where the error is
rather than using real IDs.

Also note that the mock return methods taxonomy IDs _must_ match the fake taxonomy IDs we used in the Blast Hits.
"""

import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
import pytest

from eplace_lib.taxonomy import (
    TaxonomyExtractor,
    SequenceExtractor,
    process_blast_results_for_taxonomy,
    rewrite_blast_hits,
    generate_classification_summary,
    _subject_id_matches,
)
from eplace_lib.blast_analysis import BlastHit


class TestTaxonomyExtractor:
    """Test cases for TaxonomyExtractor class."""
    
    def setup_method(self):
        self.taxonomy_extractor = TaxonomyExtractor()
        self.salmonella_taxonomy = {
            'domain': ('2', 'Bacteria'),
            'phylum': ('1224', 'Pseudomonadota'),
            'class': ('1236', 'Gammaproteobacteria'),
            'order': ('91347', 'Enterobacterales'),
            'family': ('543', 'Enterobacteriaceae'),
            'genus': ('590', 'Salmonella')
        }
        self.human_taxonomy = {
            'domain': ('2759', 'Eukaryota'),
            'phylum': ('7711', 'Chordata'),
            'class': ('40674', 'Mammalia'),
            'order': ('9443', 'Primates'),
            'family': ('9604', 'Hominidae'),
            'genus': ('9605', 'Homo'),
            'species': ('9606', 'Homo sapiens')
        }
        self.pan_taxonomy = {
            'domain': ('2759', 'Eukaryota'),
            'phylum': ('7711', 'Chordata'),
            'class': ('40674', 'Mammalia'),
            'order': ('9443', 'Primates'),
            'family': ('9604', 'Hominidae'),
            'genus': ('9596', 'Pan'),
            'species': ('9597', 'Pan paniscus')
        }
        self.hits = [
            BlastHit(
                query_id='seq1', subject_id='gi|156763568|gb|EU014687.1|',
                percent_identity=100.000, 
                alignment_length=540,
                query_length=540,
                subject_length=1432,
                query_start=1, 
                query_end=540,
                subject_start=1,
                subject_end=540,
                evalue=0.0,
                bit_score=998,
                query_coverage=100,
                subject_taxid="590",
                subject_taxids="590",
                subject_taxonomy = self.salmonella_taxonomy
            ),
            BlastHit(
                query_id='seq2', subject_id='gi|34190046|gb|BC014593.2|',
                percent_identity=100.000, 
                alignment_length=420,
                query_length=420,
                subject_length=784,
                query_start=1, 
                query_end=420,
                subject_start=1,
                subject_end=420,
                evalue=0.0,
                bit_score=776,
                query_coverage=100,
                subject_taxid="9606",
                subject_taxids="9606;9605",
                subject_taxonomy = self.human_taxonomy
            ),
            BlastHit(
                query_id='seq3', subject_id='gi|2694387494|ref|XM_055113774.3|',
                percent_identity=91.304, 
                alignment_length=115,
                query_length=420,
                subject_length=1626,
                query_start=218, 
                query_end=331,
                subject_start=173,
                subject_end=286,
                evalue=3.56e-33,
                bit_score=156,
                query_coverage=27.3809523809524,
                subject_taxid="9597",
                subject_taxids="9597;9596;9604",
                subject_taxonomy=self.pan_taxonomy
            )
        ]

    def test_parse_taxids(self):
        tax_info = self.taxonomy_extractor.parse_taxids([9606, 590])

        assert tax_info["9606"] == self.human_taxonomy
        assert tax_info["590"] == self.salmonella_taxonomy
    
    def test_group_hits_by_query(self):
        """Test grouping BLAST hits by query."""
        hits = self.hits
        
        extractor = self.taxonomy_extractor
        grouped = extractor.group_hits_by_query(hits)
        
        assert len(grouped) == 3
        assert len(grouped['seq1']) == 1
        assert len(grouped['seq2']) == 1
    
    def test_select_representatives_by_rank_valid(self):
        """Test selecting representatives by valid rank."""
        hits = self.hits
        
        extractor = self.taxonomy_extractor
        representatives = extractor.select_representatives_by_rank(
            hits=hits,
            rank='genus'
        )
        
        # Should return at least one representative
        assert len(representatives) >= 1
        # Best hit should be included (highest bit score)
        assert representatives[0].bit_score == 998
    
    def test_select_representatives_empty_hits(self):
        """Test selecting representatives from empty hits list."""
        extractor = self.taxonomy_extractor
        representatives = extractor.select_representatives_by_rank(
            hits=[],
            rank='genus'
        )
        
        assert len(representatives) == 0

    def test_select_representatives_no_taxonomy_fallback(self):
        """Hits with no subject_taxonomy (e.g. MMseqs2 without tax DB) are
        still returned as representatives, grouped by subject_id."""
        extractor = self.taxonomy_extractor

        hits = [
            BlastHit(
                query_id='query1', subject_id='subjectA',
                percent_identity=95.0,
                alignment_length=500,
                query_length=500,
                subject_length=1000,
                query_start=1,
                query_end=500,
                subject_start=1,
                subject_end=500,
                evalue=1e-50,
                bit_score=900,
                query_coverage=100,
                subject_taxid="0",
                subject_taxids="0",
                subject_taxonomy=None
            ),
            BlastHit(
                query_id='query1', subject_id='subjectA',
                percent_identity=90.0,
                alignment_length=500,
                query_length=500,
                subject_length=1000,
                query_start=1,
                query_end=500,
                subject_start=1,
                subject_end=500,
                evalue=1e-40,
                bit_score=800,
                query_coverage=100,
                subject_taxid="0",
                subject_taxids="0",
                subject_taxonomy=None
            ),
            BlastHit(
                query_id='query1', subject_id='subjectB',
                percent_identity=85.0,
                alignment_length=500,
                query_length=500,
                subject_length=1000,
                query_start=1,
                query_end=500,
                subject_start=1,
                subject_end=500,
                evalue=1e-30,
                bit_score=700,
                query_coverage=100,
                subject_taxid="0",
                subject_taxids="0",
                subject_taxonomy=None
            ),
        ]

        representatives = extractor.select_representatives_by_rank(
            hits=hits,
            rank='genus'
        )

        # Two unique subject_ids → two groups → two representatives
        assert len(representatives) == 2
        rep_ids = {r.subject_id for r in representatives}
        assert rep_ids == {'subjectA', 'subjectB'}
        # Best hit per group (highest bit_score) should be selected
        subjectA_rep = next(r for r in representatives if r.subject_id == 'subjectA')
        assert subjectA_rep.bit_score == 900

    def test_select_representatives_with_preferred(self):
        """Test that preferred representatives are reused when available."""
        extractor = self.taxonomy_extractor
        
        # First query with two hits from the same genus
        first_query_hits = [
            BlastHit(
                query_id='query1', subject_id='subject1',
                percent_identity=95.0, 
                alignment_length=500,
                query_length=500,
                subject_length=1000,
                query_start=1, 
                query_end=500,
                subject_start=1,
                subject_end=500,
                evalue=0.0,
                bit_score=900,
                query_coverage=100,
                subject_taxid="590",
                subject_taxids="590",
                subject_taxonomy=self.salmonella_taxonomy
            ),
            BlastHit(
                query_id='query1', subject_id='subject2',
                percent_identity=98.0, 
                alignment_length=500,
                query_length=500,
                subject_length=1000,
                query_start=1, 
                query_end=500,
                subject_start=1,
                subject_end=500,
                evalue=0.0,
                bit_score=950,  # Higher score - should be selected
                query_coverage=100,
                subject_taxid="590",
                subject_taxids="590",
                subject_taxonomy=self.salmonella_taxonomy
            ),
        ]
        
        # Select representatives for first query (no preferred yet)
        reps1 = extractor.select_representatives_by_rank(
            hits=first_query_hits,
            rank='genus'
        )
        
        # Should select subject2 (higher bit score)
        assert len(reps1) == 1
        assert reps1[0].subject_id == 'subject2'
        
        # Build preferred representatives dict
        preferred = {reps1[0].subject_taxonomy['genus'][1]: reps1[0].subject_id}
        
        # Second query with same genus but different scores
        second_query_hits = [
            BlastHit(
                query_id='query2', subject_id='subject1',
                percent_identity=95.0, 
                alignment_length=500,
                query_length=500,
                subject_length=1000,
                query_start=1, 
                query_end=500,
                subject_start=1,
                subject_end=500,
                evalue=0.0,
                bit_score=850,
                query_coverage=100,
                subject_taxid="590",
                subject_taxids="590",
                subject_taxonomy=self.salmonella_taxonomy
            ),
            BlastHit(
                query_id='query2', subject_id='subject2',
                percent_identity=92.0, 
                alignment_length=500,
                query_length=500,
                subject_length=1000,
                query_start=1, 
                query_end=500,
                subject_start=1,
                subject_end=500,
                evalue=0.0,
                bit_score=800,  # Lower score than subject1
                query_coverage=100,
                subject_taxid="590",
                subject_taxids="590",
                subject_taxonomy=self.salmonella_taxonomy
            ),
        ]
        
        # Select representatives for second query with preferred
        reps2 = extractor.select_representatives_by_rank(
            hits=second_query_hits,
            rank='genus',
            preferred_representatives=preferred
        )
        
        # Should still select subject2 (from preferred), even though subject1 has higher score
        assert len(reps2) == 1
        assert reps2[0].subject_id == 'subject2'



class TestSequenceExtractor:
    """Test cases for SequenceExtractor class."""
    
    def setup_method(self):
        self.salmonella_taxonomy = {
            'domain': ('2', 'Bacteria'),
            'phylum': ('1224', 'Pseudomonadota'),
            'class': ('1236', 'Gammaproteobacteria'),
            'order': ('91347', 'Enterobacterales'),
            'family': ('543', 'Enterobacteriaceae'),
            'genus': ('590', 'Salmonella')
        }
        self.human_taxonomy = {
            'domain': ('2759', 'Eukaryota'),
            'phylum': ('7711', 'Chordata'),
            'class': ('40674', 'Mammalia'),
            'order': ('9443', 'Primates'),
            'family': ('9604', 'Hominidae'),
            'genus': ('9605', 'Homo'),
            'species': ('9606', 'Homo sapiens')
        }
        self.pan_taxonomy = {
            'domain': ('2759', 'Eukaryota'),
            'phylum': ('7711', 'Chordata'),
            'class': ('40674', 'Mammalia'),
            'order': ('9443', 'Primates'),
            'family': ('9604', 'Hominidae'),
            'genus': ('9596', 'Pan'),
            'species': ('9597', 'Pan paniscus')
        }
        self.hits = [
            BlastHit(
                query_id='seq1', subject_id='gi|156763568|gb|EU014687.1|',
                percent_identity=100.000, 
                alignment_length=540,
                query_length=540,
                subject_length=1432,
                query_start=1, 
                query_end=540,
                subject_start=1,
                subject_end=540,
                evalue=0.0,
                bit_score=998,
                query_coverage=100,
                subject_taxid="590",
                subject_taxids="590",
                subject_taxonomy=self.salmonella_taxonomy
            ),
            BlastHit(
                query_id='seq2', subject_id='gi|34190046|gb|BC014593.2|',
                percent_identity=100.000, 
                alignment_length=420,
                query_length=420,
                subject_length=784,
                query_start=1, 
                query_end=420,
                subject_start=1,
                subject_end=420,
                evalue=0.0,
                bit_score=776,
                query_coverage=100,
                subject_taxid="9606",
                subject_taxids="9606;9605",
                subject_taxonomy=self.human_taxonomy
            ),
            BlastHit(
                query_id='seq3', subject_id='gi|2694387494|ref|XM_055113774.3|',
                percent_identity=91.304, 
                alignment_length=115,
                query_length=420,
                subject_length=1626,
                query_start=218, 
                query_end=331,
                subject_start=173,
                subject_end=286,
                evalue=3.56e-33,
                bit_score=156,
                query_coverage=27.3809523809524,
                subject_taxid="9597",
                subject_taxids="9597;9596;9604",
                subject_taxonomy=self.pan_taxonomy
            )
        ]

    
    def test_init_default(self):
        """Test SequenceExtractor initialization with defaults."""
        with patch.dict('os.environ', {}, clear=True):
            extractor = SequenceExtractor()
            assert extractor.blastdb_path == Path.home() / "blastdb"
    
    def test_init_with_env_var(self):
        """Test SequenceExtractor initialization with BLASTDB env var."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.dict('os.environ', {'BLASTDB': tmpdir}):
                extractor = SequenceExtractor()
                assert extractor.blastdb_path == Path(tmpdir)
    
    def test_init_with_path(self):
        """Test SequenceExtractor initialization with explicit path."""
        test_path = Path("/test/path")
        extractor = SequenceExtractor(blastdb_path=test_path)
        assert extractor.blastdb_path == test_path
    
    @patch('subprocess.run')
    def test_check_blastdbcmd_available_true(self, mock_run):
        """Test checking blastdbcmd availability when available."""
        mock_run.return_value = MagicMock(returncode=0)
        
        extractor = SequenceExtractor()
        assert extractor.check_blastdbcmd_available() is True
    
    @patch('subprocess.run')
    def test_check_blastdbcmd_available_false(self, mock_run):
        """Test checking blastdbcmd availability when not available."""
        mock_run.side_effect = FileNotFoundError()
        
        extractor = SequenceExtractor()
        assert extractor.check_blastdbcmd_available() is False
    
    @patch('eplace_lib.taxonomy.SequenceExtractor.check_blastdbcmd_available')
    @patch('subprocess.run')
    def test_extract_sequences_success(self, mock_run, mock_check):
        """Test extracting sequences successfully."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_fasta = tmppath / "output.fasta"
            
            mock_check.return_value = True
            mock_run.return_value = MagicMock(returncode=0, stderr="")
            
            extractor = SequenceExtractor(blastdb_path=tmppath)
            success = extractor.extract_sequences(
                sequence_ids=['seq1', 'seq2'],
                output_fasta=output_fasta
            )
            
            assert success is True
            mock_run.assert_called_once()
    
    @patch('eplace_lib.taxonomy.SequenceExtractor.check_blastdbcmd_available')
    def test_extract_sequences_not_available(self, mock_check):
        """Test extracting sequences when blastdbcmd not available."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_fasta = tmppath / "output.fasta"
            
            mock_check.return_value = False
            
            extractor = SequenceExtractor()
            with pytest.raises(RuntimeError, match="blastdbcmd is not available"):
                extractor.extract_sequences(
                    sequence_ids=['seq1', 'seq2'],
                    output_fasta=output_fasta
                )
    
    @patch('eplace_lib.taxonomy.SequenceExtractor.check_blastdbcmd_available')
    def test_extract_sequences_empty_list(self, mock_check):
        """Test extracting sequences with empty ID list."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_fasta = tmppath / "output.fasta"
            
            mock_check.return_value = True
            
            extractor = SequenceExtractor()
            success = extractor.extract_sequences(
                sequence_ids=[],
                output_fasta=output_fasta
            )
            
            assert success is False
    
    @patch('eplace_lib.taxonomy.SequenceExtractor.extract_sequences')
    def test_extract_representatives_for_query(self, mock_extract):
        """Test extracting representatives for a single query."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            hits = self.hits
            mock_extract.return_value = True
            
            extractor = SequenceExtractor()
            output_fasta = extractor.extract_representatives_for_query(
                query_id='seq1',
                representative_hits=hits,
                output_dir=tmppath
            )
            
            assert output_fasta is not None
            assert output_fasta.parent.exists()
    
    def test_extract_representatives_for_query_empty_hits(self):
        """Test extracting representatives with no hits."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            extractor = SequenceExtractor()
            output_fasta = extractor.extract_representatives_for_query(
                query_id='seq1',
                representative_hits=[],
                output_dir=tmppath
            )
            
            assert output_fasta is None


class TestProcessBlastResultsForTaxonomy:
    """Test cases for process_blast_results_for_taxonomy function."""

    def setup_method(self):
        self.salmonella_taxonomy = {
            'domain': ('2', 'Bacteria'),
            'phylum': ('1224', 'Pseudomonadota'),
            'class': ('1236', 'Gammaproteobacteria'),
            'order': ('91347', 'Enterobacterales'),
            'family': ('543', 'Enterobacteriaceae'),
            'genus': ('590', 'Salmonella')
        }
        self.human_taxonomy = {
            'domain': ('2759', 'Eukaryota'),
            'phylum': ('7711', 'Chordata'),
            'class': ('40674', 'Mammalia'),
            'order': ('9443', 'Primates'),
            'family': ('9604', 'Hominidae'),
            'genus': ('9605', 'Homo'),
            'species': ('9606', 'Homo sapiens')
        }
        self.pan_taxonomy = {
            'domain': ('2759', 'Eukaryota'),
            'phylum': ('7711', 'Chordata'),
            'class': ('40674', 'Mammalia'),
            'order': ('9443', 'Primates'),
            'family': ('9604', 'Hominidae'),
            'genus': ('9596', 'Pan'),
            'species': ('9597', 'Pan paniscus')
        }

    @patch('eplace_lib.taxonomy.TaxonomyExtractor.parse_taxids')
    @patch('eplace_lib.taxonomy.SequenceExtractor.extract_representatives_for_query')
    def test_process_blast_results_for_taxonomy(self, mock_extract, mock_parse_taxids):
        """Test processing BLAST results for taxonomy."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            hits = [
                BlastHit(
                    query_id='seq1', subject_id='gi|156763568|gb|EU014687.1|',
                    percent_identity=100.000, 
                    alignment_length=540,
                    query_length=540,
                    subject_length=1432,
                    query_start=1, 
                    query_end=540,
                    subject_start=1,
                    subject_end=540,
                    evalue=0.0,
                    bit_score=998,
                    query_coverage=100,
                    subject_taxid="1495401",
                    subject_taxids="1495402",
                    subject_taxonomy=self.salmonella_taxonomy
                ),
                BlastHit(
                    query_id='seq2', subject_id='gi|34190046|gb|BC014593.2|',
                    percent_identity=100.000, 
                    alignment_length=420,
                    query_length=420,
                    subject_length=784,
                    query_start=1, 
                    query_end=420,
                    subject_start=1,
                    subject_end=420,
                    evalue=0.0,
                    bit_score=776,
                    query_coverage=100,
                    subject_taxid="1495403",
                    subject_taxids="1495405;1495404",
                    subject_taxonomy=self.human_taxonomy
                ),
                BlastHit(
                    query_id='seq3', subject_id='gi|2694387494|ref|XM_055113774.3|',
                    percent_identity=91.304, 
                    alignment_length=115,
                    query_length=420,
                    subject_length=1626,
                    query_start=218, 
                    query_end=331,
                    subject_start=173,
                    subject_end=286,
                    evalue=3.56e-33,
                    bit_score=156,
                    query_coverage=27.3809523809524,
                    subject_taxid="1495406",
                    subject_taxids="1495407;1495408;1495409",
                    subject_taxonomy=self.pan_taxonomy
                )
            ]
            
            # Mock parse_taxids to return taxonomy and phylum info
            mock_parse_taxids.return_value = {
                    "1495401":self.salmonella_taxonomy,
                    "1495403": self.human_taxonomy,
                    "1495406": self.pan_taxonomy
            }
            
            mock_extract.return_value = tmppath / "output.fasta"
            
            results = process_blast_results_for_taxonomy(
                blast_hits=hits,
                output_dir=tmppath,
                rank='genus'
            )
            
            assert len(results) == 3
            assert 'seq1' in results
            assert 'seq2' in results
            
            # Verify phylum information is correctly set on blast hits
            # seq1 has taxid 1495401 (Salmonella) which should map to Pseudomonadota phylum
            assert hits[0].subject_taxonomy == self.salmonella_taxonomy
            
            # seq2 has taxid 9606 (Homo sapiens) which should map to Chordata phylum
            assert hits[1].subject_taxonomy == self.human_taxonomy
            
            # seq3 has taxid 9597 (Pan) which should also map to Chordata phylum
            assert hits[2].subject_taxonomy == self.pan_taxonomy

class TestRewriteBlastHits:
    """Test cases for rewrite_blast_hits function."""
    
    def setup_method(self):
        self.salmonella_taxonomy = {
            'domain': ('2', 'Bacteria'),
            'phylum': ('1224', 'Pseudomonadota'),
            'class': ('1236', 'Gammaproteobacteria'),
            'order': ('91347', 'Enterobacterales'),
            'family': ('543', 'Enterobacteriaceae'),
            'genus': ('590', 'Salmonella')
        }
        self.human_taxonomy = {
            'domain': ('2759', 'Eukaryota'),
            'phylum': ('7711', 'Chordata'),
            'class': ('40674', 'Mammalia'),
            'order': ('9443', 'Primates'),
            'family': ('9604', 'Hominidae'),
            'genus': ('9605', 'Homo'),
            'species': ('9606', 'Homo sapiens')
        }
        self.expected_fields = [
            "query_id", "subject_id", "percent_identity", "alignment_length",
            "query_length", "subject_length", "query_start", "query_end",
            "subject_start", "subject_end", "evalue", "bit_score",
            "query_coverage", "subject_taxid", "subject_taxids",
            "subject_taxonomy"
        ]

    def test_rewrite_blast_hits_with_complete_annotations(self):
        """Test writing blast hits with all fields populated."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_file = tmppath / "blast_hits.tsv"
            
            hits = [
                BlastHit(
                    query_id='seq1',
                    subject_id='gi|156763568|gb|EU014687.1|',
                    percent_identity=100.0,
                    alignment_length=540,
                    query_length=540,
                    subject_length=1432,
                    query_start=1,
                    query_end=540,
                    subject_start=1,
                    subject_end=540,
                    evalue=0.0,
                    bit_score=998.0,
                    query_coverage=100.0,
                    subject_taxid="149539",
                    subject_taxids="149539",
                    subject_taxonomy=self.salmonella_taxonomy
                )
            ]
            
            result = rewrite_blast_hits(hits, output_file, header=True)
            
            assert result is True
            assert output_file.exists()
            
            # Read and verify content
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            assert len(lines) == 2  # header + 1 data line
            
            # Check header
            header = lines[0].strip().split('\t')
            assert header == self.expected_fields
            
            # Check data line
            data = lines[1].strip().split('\t')
            assert data[0] == "seq1"
            assert data[1] == "gi|156763568|gb|EU014687.1|"
            assert data[15] == str(self.salmonella_taxonomy)
    
    def test_rewrite_blast_hits_with_none_values(self):
        """Test writing blast hits with None values in optional fields."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_file = tmppath / "blast_hits.tsv"
            
            hits = [
                BlastHit(
                    query_id='seq1',
                    subject_id='gi|156763568|gb|EU014687.1|',
                    percent_identity=100.0,
                    alignment_length=540,
                    query_length=540,
                    subject_length=1432,
                    query_start=1,
                    query_end=540,
                    subject_start=1,
                    subject_end=540,
                    evalue=0.0,
                    bit_score=998.0,
                    query_coverage=100.0,
                    subject_taxid="149539",
                    subject_taxids="149539",
                    subject_taxonomy=None
                )
            ]
            
            result = rewrite_blast_hits(hits, output_file, header=True)
            
            assert result is True
            assert output_file.exists()
            
            # Read and verify content
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            assert len(lines) == 2  # header + 1 data line
            
            # Check data line - None values should be empty strings
            # Remove only newline to preserve tab-delimited format
            data = lines[1].rstrip('\n').split('\t')
            assert len(data) == 16  # Verify all fields present
            assert data[0] == "seq1"
            assert data[15] == ""  # subject_taxonomy
    
    def test_rewrite_blast_hits_without_header(self):
        """Test writing blast hits without header."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_file = tmppath / "blast_hits.tsv"
            
            hits = [
                BlastHit(
                    query_id='seq1',
                    subject_id='gi|156763568|gb|EU014687.1|',
                    percent_identity=100.0,
                    alignment_length=540,
                    query_length=540,
                    subject_length=1432,
                    query_start=1,
                    query_end=540,
                    subject_start=1,
                    subject_end=540,
                    evalue=0.0,
                    bit_score=998.0,
                    query_coverage=100.0,
                    subject_taxid="149539",
                    subject_taxids="149539",
                    subject_taxonomy=self.salmonella_taxonomy
                )
            ]
            
            result = rewrite_blast_hits(hits, output_file, header=False)
            
            assert result is True
            assert output_file.exists()
            
            # Read and verify content
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            # No header, just data line
            assert len(lines) == 1
            
            # Check data line
            data = lines[0].strip().split('\t')
            assert data[0] == "seq1"
            assert data[1] == "gi|156763568|gb|EU014687.1|"
    
    @patch('eplace_lib.taxonomy.SequenceExtractor.extract_representatives_for_query')
    @patch('eplace_lib.taxonomy.TaxonomyExtractor.parse_taxids')
    def test_process_blast_results_reuses_representatives(self, mock_parse_taxids, mock_extract):
        """Test that the same representative sequence is used across multiple queries."""
        with (tempfile.TemporaryDirectory() as tmpdir):
            tmppath = Path(tmpdir)
            
            # Mock the taxonomy parsing to return taxonomy info
            # Map taxid to (rank_tid, rank_name) for genus level
            mock_parse_taxids.return_value = {
                    "1495501": self.salmonella_taxonomy,
                    "1495502": self.human_taxonomy
            }
            
            # Create hits where two queries match the same genus (Salmonella)
            hits = [
                # Query 1 hits - Salmonella genus
                BlastHit(
                    query_id='query1', subject_id='subject_a',
                    percent_identity=95.0, 
                    alignment_length=500,
                    query_length=500,
                    subject_length=1000,
                    query_start=1, 
                    query_end=500,
                    subject_start=1,
                    subject_end=500,
                    evalue=0.0,
                    bit_score=900,
                    query_coverage=100,
                    subject_taxid="1495501",
                    subject_taxids="1495501",
                    subject_taxonomy=self.salmonella_taxonomy
                ),
                # Query 2 hits - Also Salmonella genus
                BlastHit(
                    query_id='query1', subject_id='subject_b',
                    percent_identity=97.0, 
                    alignment_length=500,
                    query_length=500,
                    subject_length=1000,
                    query_start=1, 
                    query_end=500,
                    subject_start=1,
                    subject_end=500,
                    evalue=0.0,
                    bit_score=920,  # Higher than subject_b for this query
                    query_coverage=100,
                    subject_taxid="1495501",
                    subject_taxids="1495501",
                    subject_taxonomy=self.salmonella_taxonomy
                ),
                BlastHit(
                    query_id='query2', subject_id='subject_a',
                    percent_identity=94.0, 
                    alignment_length=500,
                    query_length=500,
                    subject_length=1000,
                    query_start=1, 
                    query_end=500,
                    subject_start=1,
                    subject_end=500,
                    evalue=0.0,
                    bit_score=880,  # Lower than subject_a for this query
                    query_coverage=100,
                    subject_taxid="1495501",
                    subject_taxids="1495501",
                    subject_taxonomy=self.salmonella_taxonomy
                ),
                BlastHit(
                    query_id='query2', subject_id='subject_b',
                    percent_identity=94.0,
                    alignment_length=500,
                    query_length=500,
                    subject_length=1000,
                    query_start=1,
                    query_end=500,
                    subject_start=1,
                    subject_end=500,
                    evalue=0.0,
                    bit_score=880,  # Lower than subject_a for this query
                    query_coverage=100,
                    subject_taxid="1495501",
                    subject_taxids="1495501",
                    subject_taxonomy=self.salmonella_taxonomy
                )
            ]
            
            # Track what representatives were extracted for each query
            extracted_representatives = {}
            
            def track_extraction(query_id: str, representative_hits: list, output_dir: Path, database: str = 'core_nt') -> Path:
                # Store the subject IDs that were extracted
                extracted_representatives[query_id] = [hit.subject_id for hit in representative_hits]
                return tmppath / f"{query_id}_output.fasta"
            
            mock_extract.side_effect = track_extraction
            
            results = process_blast_results_for_taxonomy(
                blast_hits=hits,
                output_dir=tmppath,
                rank='genus'
            )

            # All queries should have results
            assert len(results) == 2
            assert 'query1' in results
            assert 'query2' in results
            
            # Query1 should select subject_b (highest bit score)
            assert 'subject_b' in extracted_representatives['query1']
            
            # Query2 should ALSO use subject_b (reused from query1),
            # even though subject_a has a higher bit score for query2
            assert 'subject_b' in extracted_representatives['query2']
            
            # Both queries should use the same representative
            assert extracted_representatives['query1'] == extracted_representatives['query2']

    
    def test_rewrite_blast_hits_empty_list(self):
        """Test writing blast hits with empty list."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_file = tmppath / "blast_hits.tsv"
            
            hits = []
            
            result = rewrite_blast_hits(hits, output_file, header=True)
            
            assert result is True
            assert output_file.exists()
            
            # Read and verify content
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            # Should only have header, no data lines
            assert len(lines) == 1
            
            # Check header exists with all expected fields
            header = lines[0].strip().split('\t')
            assert header == self.expected_fields
    
    def test_rewrite_blast_hits_multiple_records(self):
        """Test writing multiple blast hits."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_file = tmppath / "blast_hits.tsv"
            
            hits = [
                BlastHit(
                    query_id='seq1',
                    subject_id='gi|156763568|gb|EU014687.1|',
                    percent_identity=100.0,
                    alignment_length=540,
                    query_length=540,
                    subject_length=1432,
                    query_start=1,
                    query_end=540,
                    subject_start=1,
                    subject_end=540,
                    evalue=0.0,
                    bit_score=998.0,
                    query_coverage=100.0,
                    subject_taxid="1495393",
                    subject_taxids="1495394",
                    subject_taxonomy=self.salmonella_taxonomy
                ),
                BlastHit(
                    query_id='seq2',
                    subject_id='gi|34190046|gb|BC014593.2|',
                    percent_identity=95.5,
                    alignment_length=420,
                    query_length=420,
                    subject_length=784,
                    query_start=1,
                    query_end=420,
                    subject_start=1,
                    subject_end=420,
                    evalue=1e-100,
                    bit_score=776.0,
                    query_coverage=100.0,
                    subject_taxid="9606",
                    subject_taxids="9606",
                    subject_taxonomy=self.human_taxonomy
                )
            ]
            
            result = rewrite_blast_hits(hits, output_file, header=True)
            
            assert result is True
            assert output_file.exists()
            
            # Read and verify content
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            assert len(lines) == 3  # header + 2 data lines
            
            # Check first data line
            data1 = lines[1].strip().split('\t')
            assert data1[0] == "seq1"
            assert data1[15] == str(self.salmonella_taxonomy)
            
            # Check second data line
            data2 = lines[2].strip().split('\t')
            assert data2[0] == "seq2"
            assert data2[15] == str(self.human_taxonomy)


class TestGenerateClassificationSummary:
    """Test cases for generate_classification_summary function."""

    def setup_method(self):
        self.salmonella_taxonomy = {
            'domain': ('2', 'Bacteria'),
            'phylum': ('1224', 'Pseudomonadota'),
            'class': ('1236', 'Gammaproteobacteria'),
            'order': ('91347', 'Enterobacterales'),
            'family': ('543', 'Enterobacteriaceae'),
            'genus': ('590', 'Salmonella')
        }
        self.ecoli_taxonomy = {
            'domain': ('2', 'Bacteria'),
            'phylum': ('1224', 'Pseudomonadota'),
            'class': ('1236', 'Gammaproteobacteria'),
            'order': ('91347', 'Enterobacterales'),
            'family': ('543', 'Enterobacteriaceae'),
            'genus': ('561', 'Escherichia'),
            'species': ('562', 'Escherichia coli')
        }
        # in this test taxonomy we deliberartely delete order to confirm that 
        # the missing taxonomies are handeled correctly
        self.human_taxonomy = {
            'domain': ('2759', 'Eukaryota'),
            'phylum': ('7711', 'Chordata'),
            'class': ('40674', 'Mammalia'),
            'family': ('9604', 'Hominidae'),
            'genus': ('9605', 'Homo'),
            'species': ('9606', 'Homo sapiens')
        }
        self.pan_taxonomy = {
            'domain': ('2759', 'Eukaryota'),
            'phylum': ('7711', 'Chordata'),
            'class': ('40674', 'Mammalia'),
            'order': ('9443', 'Primates'),
            'family': ('9604', 'Hominidae'),
            'genus': ('9596', 'Pan'),
            'species': ('9597', 'Pan paniscus')
        }
    
    def _get_column_indices(self, header: list) -> dict:
        """Helper method to get column indices from header row."""
        return {
            'query_id': header.index('query_id'),
            'blast_hits': header.index('blast_hits'),
            'taxonomy_blast': header.index('taxonomy_blast'),
            'blast_classification_name': header.index('blast_classification_name'),
            'blast_group_name': header.index('blast_group_name'),
            'blast_tree_label_name': header.index('blast_tree_label_name'),
            'tree_based_classification': header.index('tree_based_classification'),
            'taxonomy_tree': header.index('taxonomy_tree'),
            'tree_classification_name': header.index('tree_classification_name'),
            'tree_group_name': header.index('tree_group_name'),
            'tree_tree_label_name': header.index('tree_tree_label_name'),
            'appears_in_multiple_groups': header.index('appears_in_multiple_groups'),
            'has_classification': header.index('has_classification')
        }
    
    def test_generate_classification_summary_basic(self):
        """Test generating basic classification summary."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_file = tmppath / "classification.tsv"
            
            hits = [
                BlastHit(
                    query_id='seq1', subject_id='gi|156763568|gb|EU014687.1|',
                    percent_identity=100.0, 
                    alignment_length=540,
                    query_length=540,
                    subject_length=1432,
                    query_start=1, 
                    query_end=540,
                    subject_start=1,
                    subject_end=540,
                    evalue=0.0,
                    bit_score=998,
                    query_coverage=100,
                    subject_taxid="590",
                    subject_taxids="590",
                    subject_taxonomy=self.salmonella_taxonomy
                ),
                BlastHit(
                    query_id='seq2', subject_id='gi|34190046|gb|BC014593.2|',
                    percent_identity=100.0, 
                    alignment_length=420,
                    query_length=420,
                    subject_length=784,
                    query_start=1, 
                    query_end=420,
                    subject_start=1,
                    subject_end=420,
                    evalue=0.0,
                    bit_score=776,
                    query_coverage=100,
                    subject_taxid="9606",
                    subject_taxids="9606",
                    subject_taxonomy=self.human_taxonomy
                )
            ]

            seqs = {'seq1': 'this has a hit', 'seq2': 'so does this', 'seq3': 'a missing sequence'}
            
            result = generate_classification_summary(
                sequences=seqs,
                blast_hits=hits,
                output_file=output_file,
                rank='genus',
                group_rank='class',
                tree_label_rank='family'
            )
            
            assert result is True
            assert output_file.exists()
            
            # Read and verify content
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            # Should have header + 3 data lines
            assert len(lines) == 4
            
            # Check header and get column indices
            header = lines[0].strip().split('\t')
            assert 'query_id' in header
            assert 'blast_hits' in header
            assert 'blast_classification_name' in header
            assert 'blast_group_name' in header
            assert 'blast_tree_label_name' in header
            assert 'tree_based_classification' in header
            assert 'tree_classification_name' in header
            assert 'appears_in_multiple_groups' in header
            assert 'has_classification' in header
            
            # Get column indices using helper method
            col_idx = self._get_column_indices(header)

            # Check seq1 data (BLAST-based, no tree)
            data1 = lines[1].strip().split('\t')
            assert data1[col_idx['query_id']] == 'seq1'
            assert 'Salmonella' in data1[col_idx['blast_classification_name']]
            assert 'Gammaproteobacteria' in data1[col_idx['blast_group_name']]
            assert 'Enterobacteriaceae' in data1[col_idx['blast_tree_label_name']]
            assert data1[col_idx['taxonomy_blast']] == 'Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Salmonella;'
            assert data1[col_idx['tree_based_classification']] == 'No'
            assert data1[col_idx['appears_in_multiple_groups']] == 'No'
            assert data1[col_idx['has_classification']] == 'Yes'
            
            # Check seq2 data (BLAST-based, no tree)
            data2 = lines[2].strip().split('\t')
            assert data2[col_idx['query_id']] == 'seq2'
            assert 'Homo' in data2[col_idx['blast_classification_name']]
            assert 'Mammalia' in data2[col_idx['blast_group_name']]
            assert 'Hominidae' in data2[col_idx['blast_tree_label_name']]
            assert data2[col_idx['taxonomy_blast']] == 'Eukaryota;Chordata;Mammalia;;Hominidae;Homo;Homo sapiens'
            assert data2[col_idx['tree_based_classification']] == 'No'
            assert data2[col_idx['appears_in_multiple_groups']] == 'No'
            assert data2[col_idx['has_classification']] == 'Yes'

            data3 = lines[3].strip().split('\t')
            assert data3[col_idx['query_id']] == 'seq3'
            assert 'N/A' in data3[col_idx['blast_classification_name']]
            assert 'N/A' in data3[col_idx['blast_group_name']]
            assert 'N/A' in data3[col_idx['blast_tree_label_name']]
            assert data3[col_idx['taxonomy_blast']] == ';;;;;;'
            assert data3[col_idx['tree_based_classification']] == 'No'
            assert data3[col_idx['appears_in_multiple_groups']] == 'No'
            assert data3[col_idx['has_classification']] == 'No'

    
    def test_generate_classification_summary_multiple_groups(self):
        """Test detecting sequences appearing in multiple groups."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_file = tmppath / "classification.tsv"
            
            # Create hits where one query has multiple hits from different classes
            hits = [
                BlastHit(
                    query_id='seq1', subject_id='hit1',
                    percent_identity=95.0, 
                    alignment_length=500,
                    query_length=500,
                    subject_length=1000,
                    query_start=1, 
                    query_end=500,
                    subject_start=1,
                    subject_end=500,
                    evalue=0.0,
                    bit_score=900,
                    query_coverage=100,
                    subject_taxid="590",
                    subject_taxids="590",
                    subject_taxonomy=self.salmonella_taxonomy
                ),
                BlastHit(
                    query_id='seq1', subject_id='hit2',
                    percent_identity=93.0, 
                    alignment_length=500,
                    query_length=500,
                    subject_length=1000,
                    query_start=1, 
                    query_end=500,
                    subject_start=1,
                    subject_end=500,
                    evalue=0.0,
                    bit_score=850,
                    query_coverage=100,
                    subject_taxid="9606",
                    subject_taxids="9606",
                    subject_taxonomy=self.human_taxonomy
                )
            ]

            seqs = {'seq1': 'this has a hit', 'seq2': 'so does this', 'seq3': 'a missing sequence'}
            
            result = generate_classification_summary(
                sequences=seqs,
                blast_hits=hits,
                output_file=output_file,
                rank='genus',
                group_rank='class',
                tree_label_rank='genus'
            )
            
            assert result is True
            
            # Read and verify content
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            # Should have header + 3 data lines (one per sequence)
            assert len(lines) == 4
            
            # Check data using helper method for column indices
            header = lines[0].strip().split('\t')
            col_idx = self._get_column_indices(header)
            
            data = lines[1].strip().split('\t')
            assert data[col_idx['query_id']] == 'seq1'
            assert data[col_idx['blast_hits']] == '2'
            assert data[col_idx['appears_in_multiple_groups']] == 'Yes'
            # Group name should be semicolon-separated and sorted
            group_names = data[col_idx['blast_group_name']].split('; ')
            assert len(group_names) == 2
            assert 'Gammaproteobacteria' in group_names
            assert 'Mammalia' in group_names
            # Verify they are sorted
            assert group_names == sorted(group_names)

    
    def test_generate_classification_summary_invalid_rank(self):
        """Test with invalid taxonomic rank."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_file = tmppath / "classification.tsv"
            
            hits = [
                BlastHit(
                    query_id='seq1', subject_id='hit1',
                    percent_identity=95.0, 
                    alignment_length=500,
                    query_length=500,
                    subject_length=1000,
                    query_start=1, 
                    query_end=500,
                    subject_start=1,
                    subject_end=500,
                    evalue=0.0,
                    bit_score=900,
                    query_coverage=100,
                    subject_taxid="590",
                    subject_taxids="590",
                    subject_taxonomy=self.salmonella_taxonomy
                )
            ]
            seqs = {'seq1': 'this has a hit', 'seq2': 'so does this', 'seq3': 'a missing sequence'}

            result = generate_classification_summary(
                sequences=seqs,
                blast_hits=hits,
                output_file=output_file,
                rank='invalid_rank',
                group_rank='class',
                tree_label_rank='genus'
            )
            
            assert result is False
    
    def test_generate_classification_summary_with_tree(self):
        """Test classification using tree-based nearest neighbor alongside BLAST."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_file = tmppath / "classification.tsv"
            
            # Create mock BLAST hits
            # subject1 has higher bit score (900) but subject2 is closer in tree
            hits = [
                BlastHit(
                    query_id='seq1', subject_id='subject1',
                    percent_identity=98.0, 
                    alignment_length=540,
                    query_length=540,
                    subject_length=540,
                    query_start=1, 
                    query_end=540,
                    subject_start=1,
                    subject_end=540,
                    evalue=0.0,
                    bit_score=900,
                    query_coverage=100,
                    subject_taxid="590",
                    subject_taxids="590",
                    subject_taxonomy=self.salmonella_taxonomy
                ),
                BlastHit(
                    query_id='seq1', subject_id='subject2',
                    percent_identity=95.0, 
                    alignment_length=540,
                    query_length=540,
                    subject_length=540,
                    query_start=1, 
                    query_end=540,
                    subject_start=1,
                    subject_end=540,
                    evalue=1e-100,
                    bit_score=800,
                    query_coverage=100,
                    subject_taxid="562",
                    subject_taxids="562",
                    subject_taxonomy=self.ecoli_taxonomy
                )
            ]
            
            seqs = {'seq1': 'this has a hit'}
            
            # Create a tree where subject2 is phylogenetically closer to seq1
            # Tree: (seq1:0.1,(subject2:0.05,subject1:0.5):0.05);
            tree_file = tmppath / "seq1.treefile"
            tree_content = "(seq1:0.1,(subject2:0.05,subject1:0.5):0.05);"
            tree_file.write_text(tree_content)
            
            # Generate classification with tree info
            result = generate_classification_summary(
                sequences=seqs,
                blast_hits=hits,
                output_file=output_file,
                rank='genus',
                group_rank='class',
                tree_label_rank='family',
                tree_files={'seq1': tree_file}
            )
            
            assert result is True
            assert output_file.exists()
            
            # Read and verify content
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            # Should have header + 1 data line
            assert len(lines) == 2
            
            # Check that we have both BLAST and tree-based classifications
            header = lines[0].strip().split('\t')
            col_idx = self._get_column_indices(header)
            data = lines[1].strip().split('\t')
            
            # BLAST-based classification should be Salmonella (subject1 - best bit score)
            assert 'Salmonella' in data[col_idx['blast_classification_name']]
            assert 'Gammaproteobacteria' in data[col_idx['blast_group_name']]
            assert 'Enterobacteriaceae' in data[col_idx['blast_tree_label_name']]
            
            # Tree-based classification should be E. coli (subject2 - phylogenetically closest)
            assert data[col_idx['tree_based_classification']] == 'Yes'
            assert 'Escherichia' in data[col_idx['tree_classification_name']]
            assert 'Gammaproteobacteria' in data[col_idx['tree_group_name']]
            assert 'Enterobacteriaceae' in data[col_idx['tree_tree_label_name']]
    
    def test_generate_classification_summary_with_missing_tree(self):
        """Test classification falls back gracefully when tree is not available."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_file = tmppath / "classification.tsv"
            
            hits = [
                BlastHit(
                    query_id='seq1', subject_id='subject1',
                    percent_identity=100.0, 
                    alignment_length=540,
                    query_length=540,
                    subject_length=540,
                    query_start=1, 
                    query_end=540,
                    subject_start=1,
                    subject_end=540,
                    evalue=0.0,
                    bit_score=998,
                    query_coverage=100,
                    subject_taxid="590",
                    subject_taxids="590",
                    subject_taxonomy=self.salmonella_taxonomy
                )
            ]
            
            seqs = {'seq1': 'this has a hit'}
            
            # Pass tree_files but with non-existent file
            non_existent_tree = tmppath / "nonexistent.treefile"
            
            result = generate_classification_summary(
                sequences=seqs,
                blast_hits=hits,
                output_file=output_file,
                rank='genus',
                group_rank='class',
                tree_label_rank='genus',
                tree_files={'seq1': non_existent_tree}
            )
            
            assert result is True
            assert output_file.exists()
            
            # Should have BLAST-based classification, but no tree-based
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            header = lines[0].strip().split('\t')
            col_idx = self._get_column_indices(header)
            data = lines[1].strip().split('\t')
            
            # BLAST-based classification should have Salmonella
            assert 'Salmonella' in data[col_idx['blast_classification_name']]
            assert 'Gammaproteobacteria' in data[col_idx['blast_group_name']]
            
            # Tree-based classification should be unavailable
            assert data[col_idx['tree_based_classification']] == 'No'
            assert data[col_idx['tree_classification_name']] == 'N/A'
            assert data[col_idx['tree_group_name']] == 'N/A'


class TestSubjectIdMatches:
    """Test cases for the _subject_id_matches helper."""

    def test_exact_match(self):
        """Two identical IDs always match."""
        assert _subject_id_matches("HQ641676.1", "HQ641676.1")

    def test_exact_custom_pipe_id(self):
        """Custom pipe-delimited IDs with the same value match exactly."""
        assert _subject_id_matches("sampleA|42", "sampleA|42")

    def test_different_custom_pipe_ids_do_not_match(self):
        """Custom IDs sharing only a trailing segment must not be conflated (Codex P2 fix)."""
        assert not _subject_id_matches("sampleA|42", "sampleB|42")

    def test_ncbi_gi_format_vs_plain_accession(self):
        """gi|...|gb|ACC| matches plain ACC via normalized comparison."""
        assert _subject_id_matches("gi|336317909|gb|HQ641676.1|", "HQ641676.1")

    def test_plain_accession_vs_ncbi_gi_format(self):
        """Plain ACC matches gi|...|gb|ACC| (commutative)."""
        assert _subject_id_matches("HQ641676.1", "gi|336317909|gb|HQ641676.1|")

    def test_mafft_leading_r_prefix(self):
        """_R_ACC in tree labels matches the plain ACC in BLAST hits."""
        assert _subject_id_matches("HQ641676.1", "_R_HQ641676.1")

    def test_mafft_trailing_r_suffix(self):
        """ACC_R_ in tree labels matches the plain ACC in BLAST hits."""
        assert _subject_id_matches("HQ641676.1", "HQ641676.1_R_")

    def test_different_plain_accessions_do_not_match(self):
        """Two different plain accessions are not equal."""
        assert not _subject_id_matches("HQ641676.1", "MZ387488.1")


