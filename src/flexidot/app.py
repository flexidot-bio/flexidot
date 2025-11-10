import logging
import sys

from matplotlib import rc as mplrc
from matplotlib import rcParams
import pylab as P

from flexidot._version import __version__
from flexidot.plotting import pairdotplot, polydotplot, selfdotplot
from flexidot.utils.args import parse_args
from flexidot.utils.checks import check_kmer_length, print_summary
from flexidot.utils.file_handling import read_gff_color_config
from flexidot.utils.logs import init_logging

# Matplotlib settings

# Switch to non-interactive backend to avoid _tkinter.TclError on CentOs 7 servers see Github Issue #5
P.switch_backend('agg')

# Font settings
mplrc('pdf', fonttype=42, compression=0)

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = [
    'Helvetica',
    'Verdana',
    'Tahoma',
    'DejaVu Sans',
    'Droid Sans Mono',
    'Sans',
    'Liberation',
    'Ubuntu',
    'Arial',
]

__citation__ = (
    'Please remember to cite FlexiDot as follows:\n\n'
    'Kathrin M Seibt, Thomas Schmidt, Tony Heitkam,\n'
    'FlexiDot: highly customizable, ambiguity-aware dotplots for visual sequence analyses,\n'
    'Bioinformatics, Volume 34, Issue 20, October 2018, Pages 3575–3577,\n'
    'https://doi.org/10.1093/bioinformatics/bty395'
)

###############################
#        Function Call        #
###############################


def main():
    # Parse command line arguments
    args = parse_args()

    # Set up logging
    init_logging(loglevel=args.loglevel, logfile=args.logfile)

    # Print summary of arguments
    print_summary(args)

    # Log version and command line arguments if debug is enabled
    logging.debug('FlexiDot version: %s' % __version__)
    logging.debug('  ${0}\n\n'.format(' '.join(sys.argv)))
    logging.debug('Command line arguments: %s' % args)

    # Check valid kmer length
    check_kmer_length(args.wordsize)

    # Set up variables
    alphabetic_sorting = args.sort
    convert_wobbles = args.wobble_conversion
    filetype = args.filetype
    gff = args.gff
    gff_color_config_file = args.gff_color_config
    input_user_matrix_file = args.user_matrix_file
    label_size = args.label_size
    lcs_shading_interval_len = args.lcs_shading_interval_len
    lcs_shading_num = args.lcs_shading_num
    lcs_shading_ori = args.lcs_shading_ori
    lcs_shading_ref = args.lcs_shading_ref
    lcs_shading = args.lcs_shading
    length_scaling = args.length_scaling
    max_N_percentage = args.max_n
    mirror_y_axis = args.mirror_y_axis
    modes = args.mode
    multi = args.collage
    ncols = args.n_col
    nrows = args.n_row
    only_vs_first_seq = args.only_vs_first_seq
    plot_size = args.plot_size
    prefix = f'{args.outdir}/{args.output_prefix}'
    norevcomp = args.norev
    seq_list = args.infiles
    spacing = args.spacing
    substitution_count = args.substitution_count
    title_clip_pos = (
        'B'  # Note: This was processed out of title_length in previouis versions
    )
    title_length = args.title_length
    user_matrix_print = args.user_matrix_print
    wordsize = args.wordsize
    line_col_rev = args.line_col_rev
    line_col_for = args.line_col_for
    line_width = args.line_width

    # Set True if nucleotide sequence
    if args.type_seq == 'nuc':
        type_nuc = True
    elif args.type_seq == 'aa':
        type_nuc = False

    # Set x label position
    if args.x_label_pos == 'top':
        x_label_pos_top = True
    elif args.x_label_pos == 'bottom':
        x_label_pos_top = False

    # Read gff color config file if provided
    if args.gff:
        if gff_color_config_file:
            logging.info(
                f'Reading GFF color configuration file: {gff_color_config_file}'
            )
            gff_feat_colors = read_gff_color_config(gff_color_config_file)
    else:
        gff_feat_colors = {}
        if gff_color_config_file:
            logging.warning(
                f'Provide GFF annotation files to use configuration file: {gff_color_config_file}'
            )

    # If color is set to white, reverse complementary matches are skipped
    if norevcomp:  # if norev is set
        line_col_rev = 'white'  # reverse matches not calculated

    if not type_nuc and not norevcomp:
        logging.warning('Reverse complement deactivated for proteins.')
        line_col_rev = 'white'  # reverse matches not calculated
    elif not type_nuc:
        line_col_rev = 'white'

    # Log plotting modes
    mode_text = []
    mode_names = {'0': 'self', '1': 'paired', '2': 'poly'}
    for item in modes:
        mode_text.append(str(item) + ': ' + mode_names[item])

    logging.info(f'Requested plotting modes: {", ".join(mode_text)}\n\n{50 * "="}')

    # Create dotplots
    ##########################################

    # Init empty list for image file names
    list_of_png_names = []

    # self dotplots
    if '0' in modes:
        logging.info('Calling selfdotplot')
        list_of_png_names = selfdotplot(
            seq_list,
            wordsize,
            alphabetic_sorting=alphabetic_sorting,
            convert_wobbles=convert_wobbles,
            filetype=args.filetype,
            gff_color_dict=gff_feat_colors,
            gff_files=gff,
            label_size=label_size,
            line_col_rev=line_col_rev,
            line_col_for=line_col_for,
            max_N_percentage=max_N_percentage,
            mirror_y_axis=mirror_y_axis,
            multi=multi,
            ncols=ncols,
            nrows=nrows,
            plot_size=plot_size,
            prefix=prefix,
            substitution_count=substitution_count,
            title_clip_pos=title_clip_pos,
            title_length=title_length,
            type_nuc=type_nuc,
            line_width=line_width,
        )
        if list_of_png_names:
            joiner = ',\n\t\t\t'
            logging.info(
                f'\n-> Image file(s):\t{joiner.join(list_of_png_names)}\n\n{50 * "="}'
            )
        else:
            logging.warning(f'No image files were created!\n\n{50 * "="}\n')

    # paired dotplots
    if '1' in modes:
        if multi:
            logging.info('Calling pairdotplot with collage')
            list_of_png_names = pairdotplot(
                seq_list,
                wordsize,
                alphabetic_sorting=alphabetic_sorting,
                convert_wobbles=convert_wobbles,
                filetype=filetype,
                label_size=label_size,
                length_scaling=length_scaling,
                line_col_rev=line_col_rev,
                line_col_for=line_col_for,
                max_N_percentage=max_N_percentage,
                mirror_y_axis=mirror_y_axis,
                multi=multi,
                ncols=ncols,
                nrows=nrows,
                only_vs_first_seq=only_vs_first_seq,
                plot_size=plot_size,
                prefix=prefix,
                substitution_count=substitution_count,
                title_clip_pos=title_clip_pos,
                title_length=title_length,
                type_nuc=type_nuc,
                x_label_pos_top=x_label_pos_top,
                line_width=line_width,
            )
            # t1 = time_track(t1)
        else:
            if not length_scaling:
                logging.info(
                    'Pairwise dotplot with individual output files scaled by sequence length automatically.'
                )

            logging.info('Calling pairdotplot')
            list_of_png_names = pairdotplot(
                seq_list,
                wordsize,
                alphabetic_sorting=alphabetic_sorting,
                convert_wobbles=convert_wobbles,
                filetype=filetype,
                label_size=label_size,
                length_scaling=True,
                line_col_rev=line_col_rev,
                line_col_for=line_col_for,
                max_N_percentage=max_N_percentage,
                mirror_y_axis=mirror_y_axis,
                multi=multi,
                ncols=ncols,
                nrows=nrows,
                only_vs_first_seq=only_vs_first_seq,
                plot_size=plot_size,
                prefix=prefix,
                substitution_count=substitution_count,
                title_clip_pos=title_clip_pos,
                title_length=title_length,
                type_nuc=type_nuc,
                line_width=line_width,
            )
        if list_of_png_names:
            joiner = ',\n\t\t\t'
            logging.info(
                f'\n-> Image file(s):\t{joiner.join(list_of_png_names)}\n\n{50 * "="}'
            )
        else:
            logging.warning(f'No image files were created!\n\n{50 * "="}\n')

    # all-against-all dotplot
    if '2' in modes:
        logging.info('Calling polydotplot')
        list_of_png_names = polydotplot(
            seq_list,
            wordsize,
            alphabetic_sorting=alphabetic_sorting,
            convert_wobbles=convert_wobbles,
            filetype=filetype,
            gff_color_dict=gff_feat_colors,
            gff_files=gff,
            input_user_matrix_file=input_user_matrix_file,
            label_size=label_size,
            line_col_rev=line_col_rev,
            line_col_for=line_col_for,
            lcs_shading_interval_len=lcs_shading_interval_len,
            lcs_shading_num=lcs_shading_num,
            lcs_shading_ori=lcs_shading_ori,
            lcs_shading_ref=lcs_shading_ref,
            lcs_shading=lcs_shading,
            max_N_percentage=max_N_percentage,
            mirror_y_axis=mirror_y_axis,
            plot_size=plot_size,
            prefix=prefix,
            representation=args.representation,
            spacing=spacing,
            substitution_count=substitution_count,
            title_clip_pos=title_clip_pos,
            title_length=title_length,
            type_nuc=type_nuc,
            user_matrix_print=user_matrix_print,
            line_width=line_width,
        )

        if list_of_png_names:
            joiner = ',\n\t\t\t'
            logging.info(
                f'\n-> Image file(s):\t{joiner.join(list_of_png_names)}\n\n{50 * "="}'
            )
        else:
            logging.warning(f'No image files were created!\n\n{50 * "="}\n')

    logging.info(f'\nFinished! Thank you for using FlexiDot.\n\n{__citation__}')


######################
# FlexiDot Execution #
######################
if __name__ == '__main__':
    main()
