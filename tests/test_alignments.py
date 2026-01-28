"""Tests for alignment parsing and processing functions."""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from flexidot.utils.alignments import (
    alignments_to_coordinates,
    filter_redundant_alignments,
    get_sequence_ids_from_alignments,
    load_alignments,
    parse_blast6,
    parse_paf,
)


# Test data paths
TEST_DATA_DIR = Path(__file__).parent / 'test-data'
BLAST6_FILE = TEST_DATA_DIR / 'test_alignments.blast6'
PAF_FILE = TEST_DATA_DIR / 'test_alignments.paf'


class TestParseBlast6:
    """Tests for parse_blast6 function."""

    def test_parse_basic_blast6(self):
        """Test parsing a valid BLAST6 file."""
        alignments = parse_blast6(BLAST6_FILE)
        assert len(alignments) > 0
        # Check that required fields are present
        for aln in alignments:
            assert 'query_id' in aln
            assert 'subject_id' in aln
            assert 'query_start' in aln
            assert 'query_end' in aln
            assert 'subject_start' in aln
            assert 'subject_end' in aln
            assert 'identity' in aln
            assert 'length' in aln
            assert 'strand' in aln

    def test_parse_blast6_identity_filter(self):
        """Test filtering alignments by minimum identity."""
        all_alignments = parse_blast6(BLAST6_FILE, min_identity=0)
        filtered_alignments = parse_blast6(BLAST6_FILE, min_identity=92)
        assert len(filtered_alignments) < len(all_alignments)
        for aln in filtered_alignments:
            assert aln['identity'] >= 92

    def test_parse_blast6_length_filter(self):
        """Test filtering alignments by minimum length."""
        all_alignments = parse_blast6(BLAST6_FILE, min_identity=0)
        filtered_alignments = parse_blast6(BLAST6_FILE, min_length=75)
        assert len(filtered_alignments) < len(all_alignments)
        for aln in filtered_alignments:
            assert aln['length'] >= 75

    def test_parse_blast6_file_not_found(self):
        """Test that FileNotFoundError is raised for missing file."""
        with pytest.raises(FileNotFoundError):
            parse_blast6('nonexistent_file.blast6')

    def test_parse_blast6_strand_detection(self):
        """Test that strand is correctly determined from subject coordinates."""
        # Create a test file with reverse strand alignment
        with tempfile.NamedTemporaryFile(mode='w', suffix='.blast6', delete=False) as f:
            # Forward strand: sstart < send
            f.write('seq1\tseq2\t95\t100\t5\t0\t1\t100\t1\t100\t1e-50\t180\n')
            # Reverse strand: sstart > send
            f.write('seq1\tseq2\t90\t100\t10\t0\t1\t100\t200\t101\t1e-40\t160\n')
            f.flush()

            alignments = parse_blast6(f.name)
            assert len(alignments) == 2
            assert alignments[0]['strand'] == '+'
            assert alignments[1]['strand'] == '-'
            # Check that coordinates are normalized for reverse strand
            assert alignments[1]['subject_start'] < alignments[1]['subject_end']


class TestParsePaf:
    """Tests for parse_paf function."""

    def test_parse_basic_paf(self):
        """Test parsing a valid PAF file."""
        alignments = parse_paf(PAF_FILE)
        assert len(alignments) > 0
        # Check that required fields are present
        for aln in alignments:
            assert 'query_id' in aln
            assert 'subject_id' in aln
            assert 'query_start' in aln
            assert 'query_end' in aln
            assert 'subject_start' in aln
            assert 'subject_end' in aln
            assert 'identity' in aln
            assert 'length' in aln
            assert 'strand' in aln

    def test_parse_paf_coordinate_conversion(self):
        """Test that PAF 0-based coordinates are converted to 1-based."""
        # PAF uses 0-based half-open coordinates
        with tempfile.NamedTemporaryFile(mode='w', suffix='.paf', delete=False) as f:
            # 0-based: [0, 100) should become 1-based: [1, 100]
            f.write('seq1\t1000\t0\t100\t+\tseq2\t1200\t0\t100\t95\t100\t60\n')
            f.flush()

            alignments = parse_paf(f.name)
            assert len(alignments) == 1
            assert alignments[0]['query_start'] == 1
            assert alignments[0]['query_end'] == 100
            assert alignments[0]['subject_start'] == 1
            assert alignments[0]['subject_end'] == 100

    def test_parse_paf_identity_filter(self):
        """Test filtering alignments by minimum identity."""
        all_alignments = parse_paf(PAF_FILE, min_identity=0)
        filtered_alignments = parse_paf(PAF_FILE, min_identity=92)
        assert len(filtered_alignments) < len(all_alignments)
        for aln in filtered_alignments:
            assert aln['identity'] >= 92

    def test_parse_paf_length_filter(self):
        """Test filtering alignments by minimum length."""
        all_alignments = parse_paf(PAF_FILE, min_identity=0)
        filtered_alignments = parse_paf(PAF_FILE, min_length=75)
        assert len(filtered_alignments) < len(all_alignments)
        for aln in filtered_alignments:
            assert aln['length'] >= 75

    def test_parse_paf_file_not_found(self):
        """Test that FileNotFoundError is raised for missing file."""
        with pytest.raises(FileNotFoundError):
            parse_paf('nonexistent_file.paf')

    def test_parse_paf_strand(self):
        """Test that strand is correctly parsed from PAF."""
        alignments = parse_paf(PAF_FILE)
        strands = [aln['strand'] for aln in alignments]
        assert '+' in strands
        assert '-' in strands


class TestFilterRedundantAlignments:
    """Tests for filter_redundant_alignments function."""

    def test_filter_removes_duplicates(self):
        """Test that redundant alignments are removed."""
        alignments = [
            {
                'query_id': 'seqA',
                'subject_id': 'seqB',
                'query_start': 1,
                'query_end': 100,
                'subject_start': 1,
                'subject_end': 100,
                'strand': '+',
                'identity': 95,
                'length': 100,
            },
            {
                'query_id': 'seqB',
                'subject_id': 'seqA',
                'query_start': 1,
                'query_end': 100,
                'subject_start': 1,
                'subject_end': 100,
                'strand': '+',
                'identity': 95,
                'length': 100,
            },
        ]
        filtered = filter_redundant_alignments(alignments)
        assert len(filtered) == 1

    def test_filter_keeps_self_alignments(self):
        """Test that self-alignments are kept."""
        alignments = [
            {
                'query_id': 'seqA',
                'subject_id': 'seqA',
                'query_start': 1,
                'query_end': 50,
                'subject_start': 51,
                'subject_end': 100,
                'strand': '+',
                'identity': 90,
                'length': 50,
            },
        ]
        filtered = filter_redundant_alignments(alignments)
        assert len(filtered) == 1

    def test_filter_keeps_different_coordinates(self):
        """Test that alignments with different coordinates are kept."""
        alignments = [
            {
                'query_id': 'seqA',
                'subject_id': 'seqB',
                'query_start': 1,
                'query_end': 100,
                'subject_start': 1,
                'subject_end': 100,
                'strand': '+',
                'identity': 95,
                'length': 100,
            },
            {
                'query_id': 'seqA',
                'subject_id': 'seqB',
                'query_start': 200,
                'query_end': 300,
                'subject_start': 200,
                'subject_end': 300,
                'strand': '+',
                'identity': 95,
                'length': 100,
            },
        ]
        filtered = filter_redundant_alignments(alignments)
        assert len(filtered) == 2


class TestAlignmentsToCoordinates:
    """Tests for alignments_to_coordinates function."""

    def test_forward_alignment_coordinates(self):
        """Test coordinate conversion for forward strand alignment."""
        alignments = [
            {
                'query_id': 'seq1',
                'subject_id': 'seq2',
                'query_start': 10,
                'query_end': 50,
                'subject_start': 20,
                'subject_end': 60,
                'strand': '+',
                'identity': 95,
                'length': 41,
            },
        ]
        x1, y1, x2, y2 = alignments_to_coordinates(alignments, 'seq1', 'seq2', 100)
        assert len(x1) == 1
        assert len(y1) == 1
        assert len(x2) == 0
        assert len(y2) == 0
        np.testing.assert_array_equal(x1[0], [10, 50])
        np.testing.assert_array_equal(y1[0], [20, 60])

    def test_reverse_alignment_coordinates(self):
        """Test coordinate conversion for reverse strand alignment."""
        alignments = [
            {
                'query_id': 'seq1',
                'subject_id': 'seq2',
                'query_start': 10,
                'query_end': 50,
                'subject_start': 60,
                'subject_end': 20,
                'strand': '-',
                'identity': 95,
                'length': 41,
            },
        ]
        x1, y1, x2, y2 = alignments_to_coordinates(alignments, 'seq1', 'seq2', 100)
        assert len(x1) == 0
        assert len(y1) == 0
        assert len(x2) == 1
        assert len(y2) == 1

    def test_swapped_query_subject(self):
        """Test that alignments work when query and subject are swapped."""
        alignments = [
            {
                'query_id': 'seq2',
                'subject_id': 'seq1',
                'query_start': 20,
                'query_end': 60,
                'subject_start': 10,
                'subject_end': 50,
                'strand': '+',
                'identity': 95,
                'length': 41,
            },
        ]
        # Request seq1 vs seq2, but alignment is stored as seq2 vs seq1
        x1, y1, x2, y2 = alignments_to_coordinates(alignments, 'seq1', 'seq2', 100)
        assert len(x1) == 1
        # Coordinates should be swapped
        np.testing.assert_array_equal(x1[0], [10, 50])
        np.testing.assert_array_equal(y1[0], [20, 60])

    def test_no_matching_alignments(self):
        """Test when no alignments match the requested sequences."""
        alignments = [
            {
                'query_id': 'seq1',
                'subject_id': 'seq2',
                'query_start': 10,
                'query_end': 50,
                'subject_start': 20,
                'subject_end': 60,
                'strand': '+',
                'identity': 95,
                'length': 41,
            },
        ]
        x1, y1, x2, y2 = alignments_to_coordinates(alignments, 'seq3', 'seq4', 100)
        assert len(x1) == 0
        assert len(y1) == 0
        assert len(x2) == 0
        assert len(y2) == 0


class TestLoadAlignments:
    """Tests for load_alignments function."""

    def test_load_blast6_by_extension(self):
        """Test loading BLAST6 file by extension."""
        alignments = load_alignments(BLAST6_FILE)
        assert len(alignments) > 0

    def test_load_paf_by_extension(self):
        """Test loading PAF file by extension."""
        alignments = load_alignments(PAF_FILE)
        assert len(alignments) > 0

    def test_load_with_explicit_format(self):
        """Test loading with explicit format specification."""
        alignments = load_alignments(BLAST6_FILE, file_format='blast6')
        assert len(alignments) > 0

    def test_load_invalid_format(self):
        """Test that ValueError is raised for invalid format."""
        with pytest.raises(ValueError):
            load_alignments(BLAST6_FILE, file_format='invalid')

    def test_load_with_filters(self):
        """Test loading with identity and length filters."""
        all_alignments = load_alignments(BLAST6_FILE, filter_redundant=False)
        filtered = load_alignments(
            BLAST6_FILE, min_identity=95, min_length=80, filter_redundant=False
        )
        assert len(filtered) <= len(all_alignments)

    def test_load_with_redundant_filter(self):
        """Test that redundant filter is applied by default."""
        alignments_no_filter = load_alignments(BLAST6_FILE, filter_redundant=False)
        alignments_filtered = load_alignments(BLAST6_FILE, filter_redundant=True)
        # BLAST6 test file has redundant alignments
        assert len(alignments_filtered) <= len(alignments_no_filter)


class TestGetSequenceIdsFromAlignments:
    """Tests for get_sequence_ids_from_alignments function."""

    def test_extract_unique_ids(self):
        """Test extraction of unique sequence IDs."""
        alignments = [
            {'query_id': 'seq1', 'subject_id': 'seq2'},
            {'query_id': 'seq1', 'subject_id': 'seq3'},
            {'query_id': 'seq2', 'subject_id': 'seq3'},
        ]
        ids = get_sequence_ids_from_alignments(alignments)
        assert ids == {'seq1', 'seq2', 'seq3'}

    def test_empty_alignments(self):
        """Test with empty alignment list."""
        ids = get_sequence_ids_from_alignments([])
        assert ids == set()
