from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

from flexidot._version import __version__


def parse_args():
    parser = ArgumentParser(
        prog='flexidot',
        description='FlexiDot: Flexible dotplot generation tool',
        formatter_class=ArgumentDefaultsHelpFormatter,
    )

    # Input/output arguments
    parser.add_argument(
        '-i',
        '--infiles',
        required=True,
        nargs='+',
        help='Input fasta files (fasta file name or space-separated file list.)',
    )
    parser.add_argument(
        '-o',
        '--output_prefix',
        default='flexidot_output',
        help='File prefix to be added to the generated filenames.',
    )
    parser.add_argument(
        '--outdir',
        default='.',
        help='Output directory. Default: current directory.',
    )

    # Collage arguments
    parser.add_argument(
        '-c',
        '--collage',
        action='store_true',
        default=False,
        help='Combine multiple dotplots in a collage.',
    )
    parser.add_argument(
        '--n_col',
        type=int,
        default=4,
        help='Number of columns per page (if collage is ON.',
    )
    parser.add_argument(
        '--n_row',
        type=int,
        default=5,
        help='Number of rows per page (if collage is ON).',
    )

    # File type
    parser.add_argument(
        '-f',
        '--filetype',
        choices=['png', 'pdf', 'svg', 'html'],
        default='png',
        help='Output file format: png, pdf, svg, html. '
        'html produces an interactive Plotly dotplot (zoom/pan/hover); '
        'collage layout, GFF shading, and LCS/custom-matrix shading are not supported in html output.',
    )

    # Sorting
    parser.add_argument(
        '-s',
        '--sort',
        action='store_true',
        default=False,
        help='Sort sequences alphabetically by name.',
    )

    # Calculation parameters
    parser.add_argument(
        '-k',
        '--wordsize',
        type=int,
        default=10,
        help='Wordsize (kmer length) for dotplot comparison.',
    )
    parser.add_argument(
        '-m',
        '--mode',
        action='append',
        choices=['0', '1', '2'],
        required=True,
        help='Mode of FlexiDot dotplotting. 0 = self, 1 = paired, 2 = poly (matrix with all-against-all dotplots). Call -m multiple times to run multiple modes.',
    )
    parser.add_argument(
        '-t',
        '--type_seq',
        choices=['aa', 'nuc'],
        default='nuc',
        help='Biological sequence type: aa (amino acid) or nuc (nucleotide).',
    )
    parser.add_argument(
        '-w',
        '--wobble_conversion',
        action='store_true',
        default=False,
        help='Ambiguity handling for relaxed matching. Note: This may make kmer matching slower.',
    )
    parser.add_argument(
        '-S',
        '--substitution_count',
        type=int,
        default=0,
        help='Number of substitutions allowed per window.',
    )
    parser.add_argument(
        '--max_n',
        type=float,
        default=10,
        help='Maximum percentage of Ns allowed in a kmer window. Applies only if --wobble_conversion is set, else kmers with Ns are skipped. Default: 10',
    )
    parser.add_argument(
        '-r',
        '--norev',
        action='store_true',
        default=False,
        help='Do not calculate reverse complementary matches (only for nucleotide sequences.)',
    )
    parser.add_argument(
        '-O',
        '--only_vs_first_seq',
        action='store_true',
        default=False,
        help='Limit pairwise comparisons to the 1st sequence only (if plotting mode=1 paired.)',
    )

    # Pre-calculated alignments
    parser.add_argument(
        '-a',
        '--alignment_file',
        default=None,
        help='Pre-calculated alignment file (BLAST6 or PAF format). When provided, alignments from this file will be plotted instead of performing k-mer matching.',
    )
    parser.add_argument(
        '--alignment_format',
        choices=['blast6', 'paf'],
        default=None,
        help="Format of the alignment file: 'blast6' or 'paf'. If not specified, format is auto-detected from file extension.",
    )
    parser.add_argument(
        '--min_identity',
        type=float,
        default=0.0,
        help='Minimum percent identity for filtering alignments (0-100). Default: 0.0',
    )
    parser.add_argument(
        '--min_length',
        type=int,
        default=0,
        help='Minimum alignment length for filtering. Default: 0',
    )

    # Graphic formatting
    parser.add_argument('-A', '--line_width', type=float, default=1, help='Line width')

    parser.add_argument('-B', '--line_col_for', default='black', help='Line color')

    parser.add_argument(
        '-C', '--line_col_rev', default='green', help='Reverse line color'
    )

    parser.add_argument(
        '-D',
        '--x_label_pos',
        choices=['top', 'bottom'],
        default='top',
        help="Position of the X-label. Default: 'top'",
    )

    parser.add_argument('-E', '--label_size', type=int, default=10, help='Font size')

    parser.add_argument(
        '-F',
        '--spacing',
        type=float,
        default=0.04,
        help='Spacing between dotplots (if plotting mode=2 polyplot).',
    )

    parser.add_argument(
        '-L',
        '--length_scaling',
        action='store_true',
        default=False,
        help='Scale plot size for pairwise comparison.',
    )

    parser.add_argument(
        '-M',
        '--mirror_y_axis',
        action='store_true',
        default=False,
        help='Flip y-axis (bottom-to-top or top-to-bottom)',
    )

    parser.add_argument('-P', '--plot_size', type=int, default=10, help='Plot size')

    parser.add_argument(
        '-R',
        '--representation',
        type=int,
        choices=[0, 1, 2],
        default=0,
        help='Region of plot to display. Only if plotting mode is 2: polyplot\n\
                                    0 = full [default]\n\
                                    1 = upper\n\
                                    2 = lower',
    )

    parser.add_argument(
        '-T',
        '--title_length',
        type=int,
        default=50,
        help='Limit title length for comparisons. Default: 50 characters',
    )

    # GFF shading
    parser.add_argument(
        '-g',
        '--gff',
        nargs='+',
        default=None,
        help='GFF3 files for markup in self-, paired, and poly dotplots. Provide a space-delimited list of GFF files.',
    )
    parser.add_argument(
        '-G',
        '--gff_color_config',
        default=None,
        help='Config file for custom GFF shading.',
    )

    # longest common subsequence (LCS) shading
    parser.add_argument(
        '-x',
        '--lcs_shading',
        action='store_true',
        default=False,
        help='Shade subdotplot based on longest common subsequence (LCS).',
    )

    parser.add_argument(
        '-X',
        '--lcs_shading_num',
        type=int,
        default=5,
        help='Number of shading intervals.',
    )

    parser.add_argument(
        '-y',
        '--lcs_shading_ref',
        type=int,
        choices=[0, 1, 2],
        default=0,
        help='Reference for LCS shading.\n\
                        0 = maximal LCS length [default]\n\
                        1 = maximally possible length (length of shorter sequence in pairwise comparison) \n\
                        2 = given interval sizes - DNA [default 100 bp] or proteins [default 10 aa] - see -Y',
    )

    parser.add_argument(
        '-Y',
        '--lcs_shading_interval_len',
        type=int,
        default=50,
        help='Length of intervals for LCS shading (only if --lcs_shading_ref=2) [default for nucleotides = 50; default for amino acids = 10]',
    )

    parser.add_argument(
        '-z',
        '--lcs_shading_ori',
        type=int,
        choices=[0, 1, 2],
        default=0,
        help='Shade subdotplots based on LCS\n\
                            0 = forward [default]\n\
                            1 = reverse, or\n\
                            2 = both strands (forward shading above diagonal, reverse shading on diagonal and below; if using --user_matrix_file, best LCS is used below diagonal)',
    )

    # Custom matrix shading
    parser.add_argument(
        '-u', '--user_matrix_file', help='Matrix file for shading above diagonal.'
    )
    parser.add_argument(
        '-U',
        '--user_matrix_print',
        action='store_true',
        default=False,
        help='Display matrix entries above diagonal.',
    )

    # Other options
    parser.add_argument(
        '--loglevel',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        help="Set logging level. Default: 'INFO'",
    )
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__),
    )
    parser.add_argument('--logfile', default=None, help='Name of log file')

    return parser.parse_args()
