"""Tests for GFF file handling functions."""

from pathlib import Path


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
