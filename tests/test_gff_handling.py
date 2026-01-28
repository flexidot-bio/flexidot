"""Tests for GFF file handling functions."""

from pathlib import Path


from flexidot.utils.alignments import load_alignments
from flexidot.utils.file_handling import read_gffs

# Test data paths
TEST_DATA_DIR = Path(__file__).parent / 'test-data'
EMPTY_GFF_FILE = TEST_DATA_DIR / 'empty.gff3'
EXAMPLE_GFF_FILE = TEST_DATA_DIR / 'example.gff3'


class TestReadGffs:
    """Tests for read_gffs function."""

    def test_read_empty_gff(self):
        """Test reading an empty GFF file with no annotation records."""
        # Should not raise an error even with empty GFF
        feat_dict = read_gffs(
            str(EMPTY_GFF_FILE),
            color_dict={'others': ('grey', 1, 0)},
            type_nuc=True,
            prefix='test',
            filetype='png',
        )
        # Should return empty dictionary
        assert feat_dict == {}

    def test_read_valid_gff(self):
        """Test reading a valid GFF file with annotations."""
        feat_dict = read_gffs(
            str(EXAMPLE_GFF_FILE),
            color_dict={
                'spacer1': ('blue', 1, 0),
                'repeat_region': ('red', 1, 0),
                'spacerzoom': ('green', 1, 0),
                'spacer2': ('yellow', 1, 0),
                'spacer3': ('purple', 1, 0),
                'others': ('grey', 1, 0),
            },
            type_nuc=True,
            prefix='test',
            filetype='png',
        )
        # Should return non-empty dictionary
        assert len(feat_dict) > 0
        assert 'Seq2' in feat_dict

    def test_read_multiple_gffs_with_empty(self):
        """Test reading multiple GFF files where one is empty."""
        feat_dict = read_gffs(
            [str(EXAMPLE_GFF_FILE), str(EMPTY_GFF_FILE)],
            color_dict={
                'spacer1': ('blue', 1, 0),
                'repeat_region': ('red', 1, 0),
                'spacerzoom': ('green', 1, 0),
                'spacer2': ('yellow', 1, 0),
                'spacer3': ('purple', 1, 0),
                'others': ('grey', 1, 0),
            },
            type_nuc=True,
            prefix='test',
            filetype='png',
        )
        # Should still work and return data from the non-empty file
        assert len(feat_dict) > 0
        assert 'Seq2' in feat_dict


class TestEmptyAlignments:
    """Tests for empty alignment file handling."""

    def test_empty_blast6_file(self, tmp_path):
        """Test reading an empty BLAST6 file."""
        empty_file = tmp_path / 'empty.blast6'
        empty_file.write_text('')

        alignments = load_alignments(str(empty_file), file_format='blast6')
        assert alignments == []

    def test_empty_paf_file(self, tmp_path):
        """Test reading an empty PAF file."""
        empty_file = tmp_path / 'empty.paf'
        empty_file.write_text('')

        alignments = load_alignments(str(empty_file), file_format='paf')
        assert alignments == []

    def test_blast6_file_with_only_comments(self, tmp_path):
        """Test reading a BLAST6 file with only comments."""
        comment_file = tmp_path / 'comments.blast6'
        comment_file.write_text('# This is a comment\n# Another comment\n')

        alignments = load_alignments(str(comment_file), file_format='blast6')
        assert alignments == []

    def test_paf_file_with_only_comments(self, tmp_path):
        """Test reading a PAF file with only comments."""
        comment_file = tmp_path / 'comments.paf'
        comment_file.write_text('# This is a comment\n# Another comment\n')

        alignments = load_alignments(str(comment_file), file_format='paf')
        assert alignments == []
