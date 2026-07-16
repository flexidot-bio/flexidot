###############################
#     Dot Plot Functions      #
###############################

import logging

import matplotlib.collections as cllct
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import numpy as np
import pylab as P

from flexidot.utils.alignments import alignments_to_coordinates
from flexidot.utils.file_handling import legend_figure, read_gffs, read_matrix, read_seq
from flexidot.utils.matching import find_match_pos_diag, find_match_pos_regex
from flexidot.utils.plotly_export import (
    save_pairdotplot_html,
    save_polydotplot_html,
    save_selfdotplot_html,
)
from flexidot.utils.utils import (
    calc_fig_ratio,
    create_color_list,
    sanitize_filename,
    shorten_name,
    unicode_name,
)


def selfdotplot(
    input_fasta,
    wordsize,
    alignments=None,
    alphabetic_sorting=False,
    convert_wobbles=False,
    filetype='png',
    gff_color_dict=None,
    gff_files=None,
    label_size=10,
    line_col_for='#000000',  # defalut black
    line_col_rev='#009243',  # default green
    line_width=1,
    max_N_percentage=10,
    mirror_y_axis=False,
    multi=True,
    ncols=4,
    nrows=5,
    plot_size=10,
    prefix=None,
    substitution_count=0,
    title_clip_pos='B',
    title_length=float('Inf'),
    type_nuc=True,
    teviewer_integration=False,
):
    """
    self-against-self dotplot
    partially from biopython cookbook
    """

    # Initialize default gff_color_dict if None
    if gff_color_dict is None:
        gff_color_dict = {'others': ('grey', 1, 0)}

    # read sequences
    if gff_files is None:
        gff_files = []
    seq_dict, sequences = read_seq(input_fasta)
    if seq_dict == {}:
        logging.warning('Failed to load sequences.')
        return []

    if type_nuc:
        aa_bp_unit = 'bp'
    else:
        aa_bp_unit = 'aa'

    if alphabetic_sorting:
        sequences = sorted(sequences)

    # check if at least one input sequence
    if len(sequences) == 0:
        text = '\n%s\n\nCreating %s selfdotplot images\n%s\n\n=>' % (
            50 * '=',
            len(sequences),
            28 * '-',
        )
        text += ' No sequences provided for selfdotplot!\n\nTerminating polydotplot!'
        logging.info(text)
        return
    elif len(sequences) == 1 and multi:
        text = '\n\nCreating collage output for single selfdotplot!'
        text += "\nRecommendation: Change to individual mode by using '--collage_output n'!\n\n"
        logging.info(text)

    if multi and (ncols == 0 or nrows == 0):
        ncols = max(ncols, 1)
        nrows = max(nrows, 1)
        text = (
            '\n\nSelfdotplot Collage: Invalid collage - correcting number of rows and columns:\n\tncols=%d, nrows=%d\n'
            % (ncols, nrows)
        )
        logging.info(text)

    if multi and ncols > len(sequences):
        ncols = len(sequences)
        nrows = 1
        text = (
            '\n\nSelfdotplot Collage: Few sequences - correcting number of rows and columns:\n\tncols=%d, nrows=%d\n'
            % (ncols, nrows)
        )
        logging.info(text)
    elif multi and ncols * (nrows - 1) > len(sequences):
        nrows = ((len(sequences) - 1) // ncols) + 1
        text = (
            '\n\nSelfdotplot Collage: Few sequences - correcting number of rows:\n\tncols=%d, nrows=%d\n'
            % (ncols, nrows)
        )
        logging.info(text)

    if multi and not (nrows == 1 and ncols == 1) and plot_size <= label_size / 2:
        label_size = plot_size * 3 // 2
        text = 'Reducing label size for better visualization to %d\n' % label_size
        logging.info(text)

    # Initialize gff_files if None
    if gff_files is None:
        gff_files = []

    # read gff annotation data if provided for shading
    if gff_files:
        text = '\n%s\n\nReading %s GFF annotation files\n%s\n\n=> %s\n' % (
            50 * '=',
            len(gff_files),
            28 * '-',
            ', '.join(gff_files),
        )
        logging.info(text)
        if prefix:
            legend_prefix = prefix + '-Selfdotplot'
        else:
            legend_prefix = 'Selfdotplot'
        feat_dict = read_gffs(
            gff_files,
            color_dict=gff_color_dict,
            type_nuc=type_nuc,
            prefix=legend_prefix,
            filetype=filetype,
        )

    log_txt = '\n%s\n\nCreating %s selfdotplot images\n%s\n\n=>' % (
        50 * '=',
        len(sequences),
        28 * '-',
    )

    # preparations for file name
    name_graph = 'Selfdotplots'
    if prefix:
        if not prefix[-1] == '-':
            prefix = prefix + '-'
    else:
        prefix = ''
    suffix = ''
    if convert_wobbles:
        suffix += '_wobbles'
    if substitution_count != 0:
        suffix += '_S%d' % substitution_count
    if multi:
        suffix += '_collage'

    # Create wordsize suffix only when not using pre-calculated alignments
    wordsize_suffix = '' if alignments is not None else '_wordsize%i' % wordsize

    # calculate fig ratios
    if not multi:
        ncols = 1
        nrows = 1
    figsize_x, figsize_y = calc_fig_ratio(ncols, nrows, plot_size)

    if filetype != 'html':
        P.cla()  # clear any prior graph
        if multi:
            P.figure(figsize=(figsize_x, figsize_y))
    if multi:
        page_counter = 1
    list_of_png_names = []

    counter = 0
    for seq_name in sequences:
        log_txt += '\n- ' + seq_name

        counter += 1
        if not multi:
            P.cla()  # clear any prior graph

        # read sequence
        seq_record = seq_dict[seq_name]
        name_seq = seq_record.id
        seq_one = seq_record.seq.upper()
        length_seq = len(seq_one)

        # get positions of matches
        if alignments is not None:
            # Use pre-calculated alignments
            x_lists, y_lists, x_lists_rc, y_lists_rc = alignments_to_coordinates(
                alignments, seq_name, seq_name, length_seq
            )
        elif substitution_count != 0:
            # print "RE"
            x_lists, y_lists, x_lists_rc, y_lists_rc = find_match_pos_regex(
                seq_one,
                seq_one,
                wordsize,
                substitution_count=substitution_count,
                convert_wobbles=convert_wobbles,
                max_N_percentage=max_N_percentage,
                type_nuc=type_nuc,
            )
        else:
            # print "DIAG",
            x_lists, y_lists, x_lists_rc, y_lists_rc = find_match_pos_diag(
                seq_one,
                seq_one,
                wordsize,
                convert_wobbles=convert_wobbles,
                max_N_percentage=max_N_percentage,
                type_nuc=type_nuc,
            )

        # interactive HTML plotting with plotly (no collage)
        if filetype == 'html':
            if multi and counter == 1:
                logging.info(
                    'Collage layout is not supported for HTML output - generating one interactive file per sequence.'
                )
            html_suffix = ''
            if convert_wobbles:
                html_suffix += '_wobbles'
            if substitution_count != 0:
                html_suffix += '_S%d' % substitution_count
            fig_name = '%s%s-%d_%s%s%s.%s' % (
                prefix,
                name_graph,
                counter,
                sanitize_filename(
                    shorten_name(
                        name_seq, max_len=title_length, title_clip_pos=title_clip_pos
                    )
                ),
                wordsize_suffix,
                html_suffix,
                filetype,
            )
            save_selfdotplot_html(
                name_seq=name_seq,
                length_seq=length_seq,
                x_lists=x_lists,
                y_lists=y_lists,
                x_lists_rc=x_lists_rc,
                y_lists_rc=y_lists_rc,
                fig_name=fig_name,
                aa_bp_unit=aa_bp_unit,
                line_col_for=line_col_for,
                line_col_rev=line_col_rev,
                line_width=line_width,
                mirror_y_axis=mirror_y_axis,
                title_length=title_length,
                title_clip_pos=title_clip_pos,
                label_size=label_size,
                seq=seq_one,
                gff_features=feat_dict.get(seq_name) if gff_files else None,
                gff_color_dict=gff_color_dict,
                teviewer_integration=teviewer_integration,
            )
            list_of_png_names.append(fig_name)
            continue

        # plotting with matplotlib
        #################################

        # combined plotting
        if multi:
            # plotting subplot with matplotlib
            ax = P.subplot(nrows, ncols, counter)  # rows, columns, plotnumber

            # shade annotated regions
            if gff_files:
                if seq_name in list(feat_dict.keys()):
                    features = feat_dict[seq_name]
                    for item in features:
                        feat_type, start, stop = item
                        feat_color, strength, zoom = gff_color_dict[feat_type.lower()]
                        start = max(0, start - zoom - 0.5)
                        stop = min(length_seq + 1, stop + zoom + 0.5)
                        width = stop - start
                        ax.add_patch(
                            patches.Rectangle(
                                (start, start),  # (x,y)
                                width,
                                width,  # width, height
                                edgecolor=None,
                                linewidth=line_width + zoom,
                                fill=True,
                                facecolor=feat_color,
                                alpha=strength,
                            )
                        )

            # collect lines
            lines = []
            color_list = []
            for x_lines, y_lines, col in [
                (x_lists_rc, y_lists_rc, line_col_rev),
                (x_lists, y_lists, line_col_for),
            ]:
                # If color is not white, add lines to plot
                if col != 'white':
                    for ldx in range(len(x_lines)):
                        lines.append(
                            [
                                (x_lines[ldx][0], y_lines[ldx][0]),
                                (x_lines[ldx][-1], y_lines[ldx][-1]),
                            ]
                        )
                        color_list.append(col)
            color_list = np.array(color_list)

            # draw lines
            lc = cllct.LineCollection(lines, colors=color_list, linewidths=line_width)
            ax.add_collection(lc)

            # format axes
            # print P.xticks()[0], P.yticks()[0]
            P.axis('scaled')  # make images quadratic
            P.xlim(0, length_seq + 1)
            if mirror_y_axis:
                P.ylim(0, length_seq + 1)  # rotate y axis (point upwards)
            else:
                P.ylim(length_seq + 1, 0)  # rotate y axis (point downwards)
            P.xlabel('[%s]' % aa_bp_unit, fontsize=label_size)
            P.ylabel('[%s]' % aa_bp_unit, fontsize=label_size)
            P.tick_params(axis='both', which='major', labelsize=label_size * 0.9)

            # # use same tick labels for x and y axis
            # tick_locs, tick_labels = P.yticks()
            # P.xticks(tick_locs)
            # P.xlim(0, length_seq+1)

            P.title(
                unicode_name(
                    shorten_name(
                        name_seq, max_len=title_length, title_clip_pos=title_clip_pos
                    )
                ),
                fontsize=label_size,
                fontweight='bold',
            )
            # P.title(unicode_name(name_seq), fontsize=label_size*1.3, fontweight='bold')

            # save figure and reinitiate if page is full
            if counter == ncols * nrows:
                # finalize layout - margins & spacing between plots
                try:
                    P.tight_layout(h_pad=0.02, w_pad=0.02)
                except Exception as e:
                    logging.info(
                        'Attention - pylab.tight_layout failed! Please check sequence names and layout settings! Error: %s'
                        % str(e)
                    )
                P.subplots_adjust(
                    hspace=0.5, wspace=0.5
                )  # space between rows - def 0.4

                # name and create output files (names derived from SEQNAME)
                fig_name = '%s%s%s%s-%.3d.%s' % (
                    prefix,
                    name_graph,
                    wordsize_suffix,
                    suffix,
                    page_counter,
                    filetype,
                )
                P.savefig(fig_name, bbox_inches='tight')
                P.close()
                P.cla()

                list_of_png_names.append(fig_name)

                counter = 0
                page_counter += 1

                P.figure(figsize=(figsize_x, figsize_y))

        # plotting separate figure files
        else:  # not multi
            P.figure(figsize=(plot_size, plot_size))  # figure size needs to be a square
            ax = P.subplot(1, 1, 1)  # rows, columns, plotnumber

            # shade annotated regions
            if gff_files:
                if seq_name in list(feat_dict.keys()):
                    features = feat_dict[seq_name]
                    for item in features:
                        feat_type, start, stop = item
                        feat_color, strength, zoom = gff_color_dict[feat_type.lower()]
                        start = max(0, start - zoom - 0.5)
                        stop = min(length_seq + 1, stop + zoom + 0.5)
                        width = stop - start
                        ax.add_patch(
                            patches.Rectangle(
                                (start, start),  # (x,y)
                                width,
                                width,  # width, height
                                edgecolor=None,
                                linewidth=line_width + zoom,
                                fill=True,
                                facecolor=feat_color,
                                alpha=strength,
                            )
                        )

            # collect lines
            lines = []
            color_list = []
            for x_lines, y_lines, col in [
                (x_lists_rc, y_lists_rc, line_col_rev),
                (x_lists, y_lists, line_col_for),
            ]:
                # If color is not white, add lines to plot
                if col != 'white':
                    for ldx in range(len(x_lines)):
                        lines.append(
                            [
                                (x_lines[ldx][0], y_lines[ldx][0]),
                                (x_lines[ldx][-1], y_lines[ldx][-1]),
                            ]
                        )
                        color_list.append(col)

            color_list = np.array(color_list)

            # draw lines
            lc = cllct.LineCollection(lines, colors=color_list, linewidths=line_width)
            ax.add_collection(lc)

            # format axes
            P.axis('scaled')  # make images quadratic
            P.xlim(0, length_seq + 1)
            if mirror_y_axis:
                P.ylim(0, length_seq + 1)  # rotate y axis (point upwards)
            else:
                P.ylim(length_seq + 1, 0)  # rotate y axis (point downwards)
            P.xlabel('[%s]' % aa_bp_unit, fontsize=label_size)
            P.ylabel('[%s]' % aa_bp_unit, fontsize=label_size)
            P.tick_params(axis='both', which='major', labelsize=label_size * 0.9)

            # # use same tick labels for x and y axis
            # tick_locs, tick_labels = P.yticks()
            # P.xticks(tick_locs)
            # P.xlim(0, length_seq+1)

            P.title(
                unicode_name(
                    shorten_name(
                        name_seq, max_len=title_length, title_clip_pos=title_clip_pos
                    )
                ),
                fontsize=label_size * 1.3,
                fontweight='bold',
            )

            # name and create output files (names derived from SEQNAME)
            fig_name = '%s%s-%d_%s%s%s.%s' % (
                prefix,
                name_graph,
                counter,
                sanitize_filename(
                    shorten_name(
                        name_seq, max_len=title_length, title_clip_pos=title_clip_pos
                    )
                ),
                wordsize_suffix,
                suffix,
                filetype,
            )
            P.savefig(fig_name, bbox_inches='tight')

            P.close()
            P.cla()  # clear any prior graph

            list_of_png_names.append(fig_name)

    if filetype != 'html' and multi and counter >= 1:
        # finalize layout - margins & spacing between plots
        try:
            P.tight_layout(h_pad=0.02, w_pad=0.02)
        except Exception:
            logging.info(
                'Attention - pylab.tight_layout failed! Please check sequence names and layout settings!'
            )
        P.subplots_adjust(hspace=0.5, wspace=0.5)  # space between rows - def 0.4

        # name and create output files (names derived from SEQNAME)
        fig_name = '%s%s%s%s-%.3d.%s' % (
            prefix,
            name_graph,
            wordsize_suffix,
            suffix,
            page_counter,
            filetype,
        )
        P.savefig(fig_name, bbox_inches='tight')
        P.close()
        P.cla()  # clear any prior graph

        list_of_png_names.append(fig_name)

    log_txt += '\n\nDrawing selfdotplots done.\n'
    logging.info(log_txt)

    return list_of_png_names


def pairdotplot(
    input_fasta,
    wordsize,
    alignments=None,
    alphabetic_sorting=False,
    convert_wobbles=False,
    filetype='png',
    gff_color_dict=None,
    gff_files=None,
    label_size=10,
    length_scaling=True,
    line_col_for='#000000',  # defalut black
    line_col_rev='#009243',  # default green
    line_width=1,
    max_N_percentage=10,
    mirror_y_axis=False,
    multi=True,
    ncols=4,
    nrows=5,
    only_vs_first_seq=False,
    plot_size=10,
    prefix=None,
    scale_delim_col='red',
    substitution_count=0,
    title_clip_pos='B',
    title_length=float('Inf'),
    type_nuc=True,
    x_label_pos_top=True,
    teviewer_integration=False,
):
    """
    pairwise dotplot (all-against-all)

    GFF annotations are shown as full-height/full-width shaded bands: a
    feature on the x-axis sequence spans the full plot height, a feature on
    the y-axis sequence spans the full plot width.
    """

    # Initialize default gff_color_dict if None
    if gff_color_dict is None:
        gff_color_dict = {'others': ('grey', 1, 0)}

    # Initialize default gff_files if None
    if gff_files is None:
        gff_files = []

    # read sequences
    seq_dict, sequences = read_seq(input_fasta)
    if seq_dict == {}:
        logging.warning('Failed to load sequences.')
        return []

    if type_nuc:
        aa_bp_unit = 'bp'
    else:
        aa_bp_unit = 'aa'

    if alphabetic_sorting:
        sequences = sorted(sequences)

    # check if at least two input sequences
    if len(sequences) < 2:
        text = '\n%s\n\nCreating %d paired dotplot image \n%s\n\n=>' % (
            50 * '=',
            len(sequences) * (len(sequences) - 1) / 2,
            36 * '-',
        )
        text += ' Please provide at least two sequences for pairdotplot!\n\nTerminating paired dotplot!'
        logging.info(text)
        return
    elif len(sequences) == 2 and multi:
        text = '\n\nCreating collage output for single pairdotplot!'
        text += "\nRecommendation: Change to individual mode by using '--collage_output n'!\n\n"
        logging.info(text)

    if multi and (ncols == 0 or nrows == 0):
        ncols = max(ncols, 1)
        nrows = max(nrows, 1)
        text = (
            '\n\nPairdotplot Collage: Invalid collage settings - correcting number of rows and columns:\n\tncols=%d, nrows=%d\n'
            % (ncols, nrows)
        )
        logging.info(text)

    if multi and ncols > len(sequences) * (len(sequences) - 1):
        ncols = len(sequences)
        nrows = 1
        text = (
            '\n\nPairdotplot Collage: Few sequences - correcting number of rows and columns:\n\tncols=%d, nrows=%d\n'
            % (ncols, nrows)
        )
        logging.info(text)
    elif multi and ncols * (nrows - 1) > len(sequences) * (len(sequences) - 1):
        nrows = ((len(sequences) - 1) // ncols) + 1
        text = (
            '\n\nPairdotplot Collage: Few sequences - correcting number of rows:\n\tncols=%d, nrows=%d\n'
            % (ncols, nrows)
        )
        logging.info(text)

    if not only_vs_first_seq:
        text = '\n%s\n\nCreating %d paired dotplot image for\n%s\n\n=>' % (
            50 * '=',
            len(sequences) * (len(sequences) - 1) / 2,
            36 * '-',
        )
        text += ',\n'.join(sequences) + '\n'
    else:
        text = (
            "\n%s\n\nCreating %d paired dotplot images against 1st sequence '%s':\n%s\n\n=>"
            % (50 * '=', len(sequences) - 1, sequences[0], 36 * '-')
        )
        text += ',\n'.join(sequences[1:]) + '\n'
    logging.info(text)

    if multi and not (nrows == 1 and ncols == 1) and plot_size <= label_size / 2:
        label_size = plot_size * 3 // 2
        text = 'Reducing label size for better visualization to %d\n' % label_size
        logging.info(text)

    # read gff annotation data if provided for shading
    if gff_files:
        text = '\n%s\n\nReading %s GFF annotation files\n%s\n\n=> %s\n' % (
            50 * '=',
            len(gff_files),
            28 * '-',
            ', '.join(gff_files),
        )
        logging.info(text)
        if prefix:
            legend_prefix = prefix + '-Pairdotplot'
        else:
            legend_prefix = 'Pairdotplot'
        feat_dict = read_gffs(
            gff_files,
            color_dict=gff_color_dict,
            type_nuc=type_nuc,
            prefix=legend_prefix,
            filetype=filetype,
        )

    y_label_rotation = 'vertical'
    # for cartesian coordinate system with mirrored y-axis: plot x labels below plot
    if mirror_y_axis:
        x_label_pos_top = False

    # preparations for file name
    name_graph = 'Pairdotplot'
    if prefix:
        if not prefix[-1] == '-':
            prefix = prefix + '-'
    else:
        prefix = ''
    suffix = ''
    if convert_wobbles:
        suffix += '_wobbles'
    if substitution_count != 0:
        suffix += '_S%d' % substitution_count
    if length_scaling:
        suffix += '_scaled'
    if multi:
        suffix += '_collage'

    # Create wordsize suffix only when not using pre-calculated alignments
    wordsize_suffix = '' if alignments is not None else '_wordsize%i' % wordsize

    # calculate fig ratios
    if not multi:
        ncols = 1
        nrows = 1
    figsize_x, figsize_y = calc_fig_ratio(ncols, nrows, plot_size)

    if filetype != 'html':
        P.cla()  # clear any prior graph
        if multi:
            P.figure(figsize=(figsize_x, figsize_y))
    list_of_png_names = []
    if multi:
        page_counter = 1

    # prepare LCS data file
    lcs_data_file = open(
        '%sPairdotplot%s_lcs_data_file%s.txt'
        % (
            prefix,
            wordsize_suffix,
            suffix.replace('_scaled', '').replace('_collage', ''),
        ),
        'w',
    )
    lcs_data_file.write(
        '\t'.join(
            [
                '#title1',
                'title2',
                'len_seq1',
                'len_seq2',
                'len_lcs_for',
                '%_min_seq_len',
                'len_lcs_rev',
                '%_min_seq_len',
            ]
        )
        + '\n'
    )

    counter, seq_counter = 0, 0
    log_txt = '\nDrawing pairwise dotplots'

    seq_text = ''
    for idx in range(len(sequences) - 1):
        logging.debug('\n%d\t%s vs.' % ((seq_counter + 1), sequences[idx]))
        seq_text += '\n%d\t%s vs.' % ((seq_counter + 1), sequences[idx])

        rec_two = seq_dict[sequences[idx]]
        name_two = rec_two.id
        seq_two = rec_two.seq
        len_two = len(seq_two)

        for jdx in range(idx + 1, len(sequences)):
            rec_one = seq_dict[sequences[jdx]]
            name_one = rec_one.id
            seq_one = rec_one.seq
            len_one = len(seq_one)

            counter += 1
            seq_counter += 1

            logging.debug(sequences[jdx])
            seq_text += ' ' + sequences[jdx]

            if not seq_counter % 25:
                log_txt += ' ' + str(seq_counter)

            # get positions of matches
            if alignments is not None:
                # Use pre-calculated alignments
                x1, y1, x2, y2 = alignments_to_coordinates(
                    alignments, name_one, name_two, len_two
                )
                # Calculate LCS from alignment lengths (use longest alignment)
                # Each element in x1/x2 is a numpy array of [start, end]
                lcs_for = max(
                    (abs(int(x[1]) - int(x[0])) for x in x1 if len(x) > 1), default=0
                )
                lcs_rev = max(
                    (abs(int(x[1]) - int(x[0])) for x in x2 if len(x) > 1), default=0
                )
            elif substitution_count != 0:
                # print "RE"
                x1, y1, x2, y2, lcs_for, lcs_rev = find_match_pos_regex(
                    seq_one,
                    seq_two,
                    wordsize,
                    substitution_count=substitution_count,
                    convert_wobbles=convert_wobbles,
                    max_N_percentage=max_N_percentage,
                    report_lcs=True,
                    type_nuc=type_nuc,
                )
            else:
                # print "DIAG"
                x1, y1, x2, y2, lcs_for, lcs_rev = find_match_pos_diag(
                    seq_one,
                    seq_two,
                    wordsize,
                    convert_wobbles=convert_wobbles,
                    max_N_percentage=max_N_percentage,
                    report_lcs=True,
                    type_nuc=type_nuc,
                )

            # write LCS data file
            lcs_data_file.write(
                '\t'.join(
                    [
                        name_one,
                        name_two,
                        str(len_one),
                        str(len_two),
                        str(lcs_for),
                        str(round((lcs_for * 100.0 / min(len_one, len_two)), 3)),
                        str(lcs_rev),
                        str(round((lcs_rev * 100.0 / min(len_one, len_two)), 3)),
                    ]
                )
                + '\n'
            )

            # interactive HTML plotting with plotly (no collage)
            if filetype == 'html':
                if multi and counter == 1:
                    logging.info(
                        'Collage layout is not supported for HTML output - generating one interactive file per pair.'
                    )
                html_suffix = ''
                if convert_wobbles:
                    html_suffix += '_wobbles'
                if substitution_count != 0:
                    html_suffix += '_S%d' % substitution_count
                fig_name = '%s%s-%d%s%s.%s' % (
                    prefix,
                    name_graph,
                    counter,
                    wordsize_suffix,
                    html_suffix,
                    filetype,
                )
                save_pairdotplot_html(
                    name_one=name_one,
                    name_two=name_two,
                    len_one=len_one,
                    len_two=len_two,
                    x1=x1,
                    y1=y1,
                    x2=x2,
                    y2=y2,
                    fig_name=fig_name,
                    aa_bp_unit=aa_bp_unit,
                    line_col_for=line_col_for,
                    line_col_rev=line_col_rev,
                    line_width=line_width,
                    mirror_y_axis=mirror_y_axis,
                    title_length=title_length,
                    title_clip_pos=title_clip_pos,
                    label_size=label_size,
                    x_label_pos_top=x_label_pos_top,
                    seq_one=seq_one,
                    seq_two=seq_two,
                    gff_features_one=feat_dict.get(name_one) if gff_files else None,
                    gff_features_two=feat_dict.get(name_two) if gff_files else None,
                    gff_color_dict=gff_color_dict,
                    teviewer_integration=teviewer_integration,
                )
                list_of_png_names.append(fig_name)
                continue

            # Plotting with matplotlib
            #################################

            # combined plotting
            if multi:
                # plotting subplot with matplotlib
                ax = P.subplot(nrows, ncols, counter)  # rows, columns, plotnumber

            else:
                # calculate figure size for separate figures
                if len_one >= len_two:
                    pass
                    # sizing = (plot_size, max(2, (plot_size) * len_two * 1.0 / len_one))
                else:
                    pass
                    # sizing = (max(2, (plot_size) * len_one * 1.0 / len_two), plot_size)
                P.figure(figsize=(plot_size, plot_size))

                ax = P.subplot(1, 1, 1)

            # shade annotated regions if gff file(s) provided
            # x-axis sequence (name_one) -> full-height vertical band
            # y-axis sequence (name_two) -> full-width horizontal band
            if gff_files:
                for seq_name, features_len, is_x_axis in [
                    (name_one, len_two, True),
                    (name_two, len_one, False),
                ]:
                    if seq_name not in list(feat_dict.keys()):
                        continue
                    for item in feat_dict[seq_name]:
                        feat_type, start, stop = item
                        feat_color, strength, zoom = gff_color_dict[feat_type.lower()]
                        start = max(0, start - zoom - 0.5)
                        stop = stop + zoom + 0.5
                        width = stop - start
                        if is_x_axis:
                            xy = (start, 0)
                            rect_w, rect_h = width, features_len + 1
                        else:
                            xy = (0, start)
                            rect_w, rect_h = features_len + 1, width
                        ax.add_patch(
                            patches.Rectangle(
                                xy,
                                rect_w,
                                rect_h,
                                edgecolor=None,
                                linewidth=line_width + zoom,
                                fill=True,
                                facecolor=feat_color,
                                alpha=strength,
                            )
                        )

            # collect lines
            lines = []
            color_list = []
            for x_lines, y_lines, col in [
                (x2, y2, line_col_rev),
                (x1, y1, line_col_for),
            ]:
                # If color is not white, add lines to plot
                if col != 'white':
                    for ldx in range(len(x_lines)):
                        lines.append(
                            [
                                (x_lines[ldx][0], y_lines[ldx][0]),
                                (x_lines[ldx][-1], y_lines[ldx][-1]),
                            ]
                        )
                        color_list.append(col)
            color_list = np.array(color_list)

            # draw lines
            lc = cllct.LineCollection(lines, colors=color_list, linewidths=line_width)
            ax.add_collection(lc)

            # format axes
            P.xlabel(
                unicode_name(
                    shorten_name(
                        name_one, max_len=title_length, title_clip_pos=title_clip_pos
                    )
                )
                + ' [%s]' % aa_bp_unit,
                fontsize=label_size,
                fontweight='bold',
                labelpad=4,
            )
            P.ylabel(
                unicode_name(
                    shorten_name(
                        name_two, max_len=title_length, title_clip_pos=title_clip_pos
                    )
                )
                + ' [%s]' % aa_bp_unit,
                fontsize=label_size,
                fontweight='bold',
                labelpad=4,
            )
            P.tick_params(axis='both', which='major', labelsize=label_size * 0.9)

            # P.axis('scaled') # make images scaled by size ### optional update ###
            if not multi:
                if length_scaling:
                    ax.set_aspect(aspect='equal', adjustable='box', anchor='NW')
                P.xlim(0, len_one + 1)
                # xlimit = [0, len_one+1]
                if mirror_y_axis:
                    P.ylim(0, len_two + 1)  # rotate y axis (point upwards)
                else:
                    P.ylim(len_two + 1, 0)  # rotate y axis (point downwards)
            elif not length_scaling:
                P.xlim(0, len_one + 1)
                # xlimit = [0, len_one+1]
                if mirror_y_axis:
                    P.ylim(0, len_two + 1)  # rotate y axis (point upwards)
                else:
                    P.ylim(len_two + 1, 0)  # rotate y axis (point downwards)
            else:
                max_len = max(len_one, len_two)
                P.xlim(0, max_len + 1)
                # xlimit = [0, max_len+1]
                if mirror_y_axis:
                    P.ylim(0, max_len + 1)  # rotate y axis (point upwards)
                else:
                    P.ylim(max_len + 1, 0)  # rotate y axis (point downwards)

                # plot line deliminating shorter sequence
                if max_len != len_one:
                    ax.plot(
                        (len_one + 1, len_one + 1),
                        (0, len_two),
                        marker='',
                        linestyle='--',
                        color=scale_delim_col,
                        markerfacecolor='r',
                    )
                if max_len != len_two:
                    ax.plot(
                        (0, len_one),
                        (len_two + 1, len_two + 1),
                        marker='',
                        linestyle='--',
                        color=scale_delim_col,
                        markerfacecolor='r',
                    )

            # # use same tick labels for x and y axis
            # if P.xlim() == P.ylim():
            #     tick_locs, tick_labels = P.yticks()
            #     P.xticks(tick_locs)
            #     P.xlim(xlimit)

            # evtl. switch x axis position
            if x_label_pos_top:
                ax.xaxis.tick_top()
                ax.xaxis.set_label_position('top')
            P.setp(ax.get_xticklabels(), fontsize=label_size * 0.9)
            P.setp(ax.get_yticklabels(), fontsize=label_size * 0.9)

            # save figure and reinitiate if page is full
            if multi and counter == ncols * nrows:
                # finalize layout - margins & spacing between plots
                try:
                    P.tight_layout(h_pad=0.02, w_pad=0.02)
                except KeyError:
                    logging.warning(
                        'Attention - pylab.tight_layout failed! Please check sequence names and layout settings!'
                    )
                if x_label_pos_top:
                    P.subplots_adjust(
                        hspace=0.5, wspace=0.5, top=0.95
                    )  # space between rows - def 0.4
                else:
                    P.subplots_adjust(
                        hspace=0.5, wspace=0.5, bottom=0.05
                    )  # space between rows - def 0.4

                # name and create output files (names derived from SEQNAME)
                fig_name = '%s%s%s%s-%.3d.%s' % (
                    prefix,
                    name_graph,
                    wordsize_suffix,
                    suffix,
                    page_counter,
                    filetype,
                )
                P.savefig(fig_name, bbox_inches='tight')
                P.close()
                P.cla()

                list_of_png_names.append(fig_name)

                counter = 0
                page_counter += 1

                P.figure(figsize=(figsize_x, figsize_y))

            # plotting separate figure files
            elif not multi:
                # finalize layout - margins & spacing between plots
                try:
                    P.tight_layout(h_pad=0.02, w_pad=0.02)
                except Exception:
                    logging.warning(
                        'Attention - pylab.tight_layout failed! Please check sequence names and layout settings!'
                    )
                if y_label_rotation == 'horizontal':
                    if x_label_pos_top:
                        P.subplots_adjust(
                            hspace=0.02, wspace=0.02, left=0.13, top=0.95
                        )  # space between rows - def 0.4
                    else:
                        P.subplots_adjust(
                            hspace=0.02, wspace=0.02, left=0.13, bottom=0.05
                        )  # space between rows - def 0.4
                else:
                    P.subplots_adjust(
                        hspace=0.02, wspace=0.02
                    )  # space between rows - def 0.4

                # name and create output files
                fig_name = '%s%s-%d%s%s.%s' % (
                    prefix,
                    name_graph,
                    counter,
                    wordsize_suffix,
                    suffix,
                    filetype,
                )
                P.savefig(fig_name)
                P.close()
                P.cla()

                list_of_png_names.append(fig_name)
                P.figure()  # clear any prior graph

        if only_vs_first_seq:
            break

    # save figure
    if filetype != 'html' and multi and counter >= 1:
        # finalize layout - margins & spacing between plots
        try:
            P.tight_layout(h_pad=0.02, w_pad=0.02)
        except Exception:
            logging.warning(
                'Attention - pylab.tight_layout failed! Please check sequence names and layout settings!'
            )
        if x_label_pos_top:
            P.subplots_adjust(
                hspace=0.5, wspace=0.5, top=0.95
            )  # space between rows - def 0.4
        else:
            P.subplots_adjust(
                hspace=0.5, wspace=0.5, bottom=0.05
            )  # space between rows - def 0.4

        # name and create output files (names derived from SEQNAME)
        fig_name = '%s%s%s%s-%.3d.%s' % (
            prefix,
            name_graph,
            wordsize_suffix,
            suffix,
            page_counter,
            filetype,
        )
        P.savefig(fig_name, bbox_inches='tight')
        P.close()
        P.cla()

        list_of_png_names.append(fig_name)

    log_txt += '\n%d done' % seq_counter
    logging.info(log_txt)

    logging.debug(seq_text)

    return list_of_png_names


def polydotplot(
    input_fasta,
    wordsize=10,
    alignments=None,
    gff_files=None,
    alphabetic_sorting=False,
    convert_wobbles=False,
    filetype='png',
    gff_color_dict=None,
    input_user_matrix_file='',
    label_size=10,
    lcs_shading_interval_len=100,
    lcs_shading_num=5,
    lcs_shading_ori=0,
    lcs_shading_ref=0,
    lcs_shading=True,
    line_col_for='#000000',  # defalut black
    line_col_rev='#009243',  # default green
    line_width=1,
    max_N_percentage=10,
    mirror_y_axis=False,
    plot_size=10,
    prefix=None,
    representation=0,
    rotate_labels=False,
    spacing=0.04,
    substitution_count=0,
    title_clip_pos='B',
    title_length=float('Inf'),
    type_nuc=True,
    user_matrix_print=True,
    x_label_pos_top=True,
    teviewer_integration=False,
):
    """
    all-against-all dotplot
    derived from dotplot function

    lcs_shading_refs:
        0 color relative to maximum lcs observed in dataset [default]
        1 color by coverage of shorter sequence (e.g. lcs = 70% of seq1)
    lcs_shading_ori
        0 forward only
        1 reverse only
        2 both orientations (in opposite plot)
    """

    # Initialize default gff_color_dict if None
    if gff_color_dict is None:
        gff_color_dict = {'others': ('grey', 1, 0)}

    # Initialize default gff_files if None
    if gff_files is None:
        gff_files = []

    # read sequences
    seq_dict, sequences = read_seq(input_fasta)
    if seq_dict == {}:
        logging.warning('Failed to load sequences.')
        return []

    if type_nuc:
        aa_bp_unit = 'bp'
    else:
        aa_bp_unit = 'aa'

    if alphabetic_sorting:
        sequences = sorted(sequences)

    if len(sequences) == 0:
        text = '\n%s\n\nCreating %dx%d polydotplot image\n%s\n\n=>' % (
            50 * '=',
            len(sequences),
            len(sequences),
            30 * '-',
        )
        text += ' No sequences provided for polydotplot!\n\nTerminating polydotplot!'
        logging.info(text)
        return
    elif len(sequences) == 1:
        text = '\n\nCreating polydotplot for single sequence!'
        text += "\nRecommendation: Use selfdotplot via '--plotting_mode 0'!\n\n"
        logging.info(text)

    text = '\n%s\n\nCreating %dx%d polydotplot image\n%s\n\n=>' % (
        50 * '=',
        len(sequences),
        len(sequences),
        30 * '-',
    )
    text += ' ' + ' '.join(sequences) + '\n'
    logging.info(text)

    # read gff annotation data if provided for shading
    if gff_files is not None and gff_files != []:
        text = '\n%s\n\nReading %s GFF annotation files\n%s\n\n=> %s\n' % (
            50 * '=',
            len(gff_files),
            28 * '-',
            ', '.join(gff_files),
        )
        logging.info(text)
        if prefix is not None and prefix != '':
            legend_prefix = prefix + '-Polydotplot'
        else:
            legend_prefix = 'Polydotplot'
        feat_dict = read_gffs(
            gff_files,
            color_dict=gff_color_dict,
            type_nuc=type_nuc,
            prefix=legend_prefix,
            filetype=filetype,
        )

    if lcs_shading and not type_nuc:
        if lcs_shading_ori != 0:
            lcs_shading_ori = 0
            text = 'Protein shading does not support reverse complementary matching!\n'
            logging.info(text)

    # read custom shading matrix & match names of sequences to fasta
    if input_user_matrix_file != '' and input_user_matrix_file is not None:
        logging.info('Reading user matrix file: %s' % input_user_matrix_file)
        # lcs_shading_ori = 2
        custom_dict = read_matrix(input_user_matrix_file)
        if custom_dict != {}:
            custom_shading = True
            custom_similarity_dict = {}
            invalid_entries = []
            custom_max = 0
            custom_min = float('Inf')
            for key in list(custom_dict.keys()):
                number_key = []

                # convert number into float
                try:
                    value = float(custom_dict[key])
                    if '.' not in custom_dict[key]:
                        value = int(custom_dict[key])
                    custom_max = max(custom_max, value)
                    custom_min = min(custom_min, value)
                except Exception:
                    value = custom_dict[key]
                    if value == '':
                        value = None
                    invalid_entries.append(key)
                # match matrix names with sequence names
                for item in key:
                    if item in sequences:
                        number_key.append(sequences.index(item))
                    else:
                        number_key.append(-1)
                # dictionary with tuple of sorted sequence indices as key and number as value
                custom_similarity_dict[tuple(sorted(number_key))] = value
            if len(invalid_entries) != 0:
                text = (
                    'No valid number in custom similarity matrix for %d entries: \n\t'
                    % (len(invalid_entries))
                )
                for key in invalid_entries:
                    text += str(key) + ' - ' + str(custom_dict[key]) + '; '
                logging.info(text[:-2] + '\n')

        text = 'Custom user matrix given: min %.2f, max %.2f\n' % (
            custom_min,
            custom_max,
        )

        # artificially rounding intervals if likely identity/divergence percentages
        if 0 <= custom_min < 1 and 0 < custom_max <= 1:
            rounding_factor = 5
            multi_factor = 100
            text += (
                ' > artificially rounding custom shading intervals: old (%.2f, %.2f) - '
                % (custom_min, custom_max)
            )
            custom_min = max(
                0,
                (multi_factor * custom_min // rounding_factor)
                * (1.0 * rounding_factor / multi_factor),
            )
            custom_max = min(
                (multi_factor * custom_max // rounding_factor)
                * (1.0 * rounding_factor / multi_factor),
                1,
            )
            text += 'new (%.2f, >%2f)\n' % (custom_min, custom_max)

        elif 0 <= custom_min < 100 and 0 < custom_max <= 100:
            rounding_factor = 5
            text += (
                ' > artificially rounding custom shading intervals: old (%.2f, %.2f) - '
                % (custom_min, custom_max)
            )
            custom_min = max(0, (custom_min // rounding_factor) * rounding_factor)
            custom_max = min((custom_max // rounding_factor) * rounding_factor, 100)
            text += 'new (%d, >%d)\n' % (custom_min, custom_max)

        logging.info(text)

    else:
        custom_shading = False

    name_graph = 'Polydotplot'
    suffix = ''
    if convert_wobbles:
        suffix += '_wobbles'
    if substitution_count != 0:
        suffix += '_S%d' % substitution_count
    if custom_shading:
        suffix += '_matrix'
    if lcs_shading:
        suffix += '_%dshades_ref%d_ori%s' % (
            lcs_shading_num + 1,
            lcs_shading_ref,
            lcs_shading_ori,
        )
        if 'ref2' in suffix and type_nuc:
            suffix = suffix.replace('ref2', '%dbp' % lcs_shading_interval_len)
        elif 'ref2' in suffix:
            suffix = suffix.replace('ref2', '%daa' % lcs_shading_interval_len)

    # Create wordsize suffix only when not using pre-calculated alignments
    wordsize_suffix = '' if alignments is not None else '_wordsize%i' % wordsize

    # name and create output files (names derived from SEQNAME)
    if prefix:
        prefix = str(prefix) + '-'
    else:
        prefix = ''

    # preparations for background shading
    if lcs_shading or custom_shading:
        # create color range white to grey
        colors = create_color_list(lcs_shading_num + 1, color_map='Greys')
        colors_2 = create_color_list(lcs_shading_num + 1, color_map='OrRd')

        if custom_shading:
            text = 'Custom Matrix Colors: ' + ', '.join(colors_2)

    # write lcs lengths to file
    lcs_data_file = open(
        '%sPolydotplot%s_lcs_data_file%s.txt'
        % (
            prefix,
            wordsize_suffix,
            suffix.replace('_scaled', '').replace('_collage', ''),
        ),
        'w',
    )
    lcs_data_file.write(
        '\t'.join(
            [
                '#title1',
                'title2',
                'len_seq1',
                'len_seq2',
                'len_lcs_for',
                '%_min_seq_len',
                'len_lcs_rev',
                '%_min_seq_len',
            ]
        )
        + '\n'
    )

    # compare sequences pairwise - save lcs and line information in dictionary for plotting
    data_dict = {}  # keys = tuple(idx, jdx), value = x1, y1, x2, y2 (line positions)
    lcs_dict = {}  # keys = tuple(idx, jdx), value = length of lcs: lcs_len or (lcs_for, lcs_rev)
    for_lcs_set = set()  # keep lengths to calculate max (excluding self comparisons)
    rev_lcs_set = set()  # keep lengths to calculate max (all)

    text = '\nTotal plot count:   %d' % (len(sequences) * (len(sequences)))
    text += '\nTotal calculations: %d' % (len(sequences) * (len(sequences) + 1) / 2)
    logging.info(text)

    logging.info(
        '\nCalculating shared regions and lengths of longest_common_substring...'
    )
    log_txt = '\nCalculating shared regions and lengths of longest_common_substring...'
    # determine  matches and length of lcs by comparing all sequence pairs

    seq_text = ''
    counter = 0
    for idx in range(len(sequences)):
        logging.debug('\n%d\t%s vs.' % ((counter + 1), sequences[idx]))
        seq_text += '\n%d\t%s vs.' % ((counter + 1), sequences[idx])
        rec_two = seq_dict[sequences[idx]]
        name_two = rec_two.id
        seq_two = rec_two.seq
        len_two = len(seq_two)

        for jdx in range(idx, len(sequences)):
            rec_one = seq_dict[sequences[jdx]]
            name_one = rec_one.id
            seq_one = rec_one.seq
            len_one = len(seq_one)

            counter += 1
            logging.debug(sequences[jdx])
            seq_text += ' ' + sequences[jdx]

            if len(sequences) < 5:
                log_txt += '\n\t%s (%d %s), %s (%d %s)' % (
                    name_one,
                    len_one,
                    aa_bp_unit,
                    name_two,
                    len_two,
                    aa_bp_unit,
                )
            else:
                if not counter % 25:
                    print(counter)
                    log_txt += str(counter)

            # Get positions of matches &  length of longest common substring based on match lengths
            if alignments is not None:
                # Use pre-calculated alignments
                x1, y1, x2, y2 = alignments_to_coordinates(
                    alignments, name_one, name_two, len_two
                )
                # Calculate LCS from alignment lengths (use longest alignment)
                lcs_for = max(
                    (abs(int(x[1]) - int(x[0])) for x in x1 if len(x) > 1), default=0
                )
                lcs_rev = max(
                    (abs(int(x[1]) - int(x[0])) for x in x2 if len(x) > 1), default=0
                )
            elif substitution_count != 0:
                # print "RE"
                x1, y1, x2, y2, lcs_for, lcs_rev = find_match_pos_regex(
                    seq_one,
                    seq_two,
                    wordsize,
                    substitution_count=substitution_count,
                    convert_wobbles=convert_wobbles,
                    max_N_percentage=max_N_percentage,
                    report_lcs=True,
                    type_nuc=type_nuc,
                )
            else:
                # print "DIAG"
                x1, y1, x2, y2, lcs_for, lcs_rev = find_match_pos_diag(
                    seq_one,
                    seq_two,
                    wordsize,
                    convert_wobbles=convert_wobbles,
                    max_N_percentage=max_N_percentage,
                    report_lcs=True,
                    type_nuc=type_nuc,
                )
            data_dict[(idx, jdx)] = x1[:], y1[:], x2[:], y2[:]
            lcs_dict[idx, jdx] = lcs_for, lcs_rev

            if idx != jdx:
                for_lcs_set.add(lcs_for)
            rev_lcs_set.add(lcs_rev)

            lcs_data_file.write(
                '\t'.join(
                    [
                        name_one,
                        name_two,
                        str(len_one),
                        str(len_two),
                        str(lcs_for),
                        str(round((lcs_for * 100.0 / min(len_one, len_two)), 3)),
                        str(lcs_rev),
                        str(round((lcs_rev * 100.0 / min(len_one, len_two)), 3)),
                    ]
                )
                + '\n'
            )

    log_txt += '\n' + str(len(sequences) * (len(sequences) + 1) / 2) + ' done'

    logging.info(log_txt)

    logging.debug('\n\nlcs_dict\n' + str(lcs_dict))
    if custom_shading:
        logging.debug('\ncustom_dict\n' + str(custom_dict))
        logging.debug('\ncustom_similarity_dict\n\n' + str(custom_similarity_dict))

    logging.info(seq_text + '\n')

    # interactive HTML plotting with plotly (no GFF/LCS/custom-matrix shading)
    if filetype == 'html':
        if gff_files:
            logging.info(
                'GFF shading is not supported for HTML output - skipping annotations.'
            )
        if lcs_shading or custom_shading:
            logging.info(
                'LCS/custom-matrix shading is not supported for HTML output - skipping shading.'
            )
        html_suffix = ''
        if convert_wobbles:
            html_suffix += '_wobbles'
        if substitution_count != 0:
            html_suffix += '_S%d' % substitution_count
        fig_name = '%s%s%s%s.%s' % (
            prefix,
            name_graph,
            wordsize_suffix,
            html_suffix,
            filetype,
        )
        save_polydotplot_html(
            sequences=sequences,
            seq_dict=seq_dict,
            data_dict=data_dict,
            fig_name=fig_name,
            aa_bp_unit=aa_bp_unit,
            line_col_for=line_col_for,
            line_col_rev=line_col_rev,
            line_width=line_width,
            mirror_y_axis=mirror_y_axis,
            title_length=title_length,
            title_clip_pos=title_clip_pos,
            label_size=label_size,
            plot_size=plot_size,
            teviewer_integration=teviewer_integration,
        )
        return [fig_name]

    if lcs_shading_ref == 2:
        color_bins = []
        text = '\nLCS lengh bins: '
        for idx in range(lcs_shading_num):
            color_bins.append(lcs_shading_interval_len * (idx + 1))
            text += ' ' + str(lcs_shading_interval_len * (idx + 1))
        logging.info(text)

    # Calculate maximum lcs length
    if lcs_shading_ori == 0:  # forward only
        if len(for_lcs_set) != 0:
            max_lcs = max(for_lcs_set)
        else:
            max_lcs = None
    elif lcs_shading_ori == 1:  # reverse complement only
        if len(rev_lcs_set) != 0:
            max_lcs = max(rev_lcs_set)
        else:
            max_lcs = None
    else:  # both orientations
        if len(rev_lcs_set) != 0 and len(for_lcs_set) != 0:
            max_lcs = max(max(rev_lcs_set), max(for_lcs_set))
        elif len(rev_lcs_set) != 0:
            max_lcs = max(rev_lcs_set)
        elif len(for_lcs_set) != 0:
            max_lcs = max(for_lcs_set)
        else:
            max_lcs = None

    if max_lcs:
        text = 'Maximum LCS: %d %s' % (max_lcs, aa_bp_unit)
        logging.info(text)
    if custom_shading:
        text = 'Maximum custom value: %d\n' % custom_max
        logging.info(text)

    # count sequences
    ncols = len(sequences)
    nrows = len(sequences)

    # get sequence lengths to scale plot widths and heights accordingly
    size_ratios = []
    for item in sequences:
        size_ratios.append(len(seq_dict[item].seq))

    P.cla()  # clear any prior graph
    # use GridSpec to resize plots according to sequence length
    if mirror_y_axis:
        height_ratios = size_ratios[::-1]
    else:
        height_ratios = size_ratios[:]
    gs = gridspec.GridSpec(
        nrows, ncols, width_ratios=size_ratios, height_ratios=height_ratios
    )
    P.figure(figsize=(plot_size, plot_size))

    # for cartesian coordinate system with mirrored y-axis: plot x labels below plot
    if mirror_y_axis and representation == 1:
        x_label_pos_top = True
    elif mirror_y_axis or representation == 2:
        x_label_pos_top = False

    # print y labels on the right, if upper right triangle is displayed
    if (representation == 1 and not mirror_y_axis) or (
        representation == 2 and mirror_y_axis
    ):
        y_label_pos = 0  # last column
    else:  # left y label
        y_label_pos = 1  # first column

    # determine label orientations
    if len(sequences) > 5 or rotate_labels:
        x_label_rotation = 45
        y_label_rotation = 'horizontal'
        if x_label_pos_top:
            xhalign = 'left'
            xvalign = 'bottom'
        else:
            xhalign = 'right'
            xvalign = 'top'
        yhalign = 'right'
    else:
        x_label_rotation = 'horizontal'
        y_label_rotation = 'vertical'
        xvalign = 'center'
        xhalign = 'center'
        yhalign = 'center'
    yvalign = 'center'

    # check combination of shading parameters for triangular output
    if (
        representation != 0 and lcs_shading and custom_shading
    ):  # both directions in triangle
        logging.info(
            '\nAttention: For triangular output custom-shading and LCS shading cannot be combined!\n'
        )
    elif (
        representation != 0 and lcs_shading and lcs_shading_ori == 2
    ):  # both directions in triangle
        logging.info(
            '\nAttention: For triangular output LCS shading for both orientations is combined to max of both orientations!\n'
        )

    log_txt = '\nDrawing polydotplot...'

    # draw subplots
    if lcs_shading and custom_shading:
        lcs_text = (
            '\n'
            + '\t'.join(
                [
                    '#Seq1',
                    'Seq2',
                    'LCS for [%s]' % aa_bp_unit,
                    'LCS for [%s]' % aa_bp_unit,
                    'Custom matrix value',
                    'Matrix color index',
                    'LCS color index',
                ]
            )
            + '\n'
        )
    elif lcs_shading:
        lcs_text = (
            '\n'
            + '\t'.join(
                [
                    '#Seq1',
                    'Seq2',
                    'LCS for [%s]' % aa_bp_unit,
                    'LCS for [%s]' % aa_bp_unit,
                    'LCS color index for',
                    'LCS color index rev',
                ]
            )
            + '\n'
        )
    elif custom_shading:
        lcs_text = (
            '\n'
            + '\t'.join(
                [
                    '#Seq1',
                    'Seq2',
                    'Custom matrix value',
                    'Color index for',
                    'Color index rev',
                ]
            )
            + '\n'
        )

    seq_text = ''

    counter, seq_counter = 0, 0
    for idx in range(len(sequences)):
        logging.debug('\n%d\t%s vs.' % ((seq_counter + 1), sequences[idx]))
        seq_text += '\n%d\t%s vs.' % ((seq_counter + 1), sequences[idx])

        rec_two = seq_dict[sequences[idx]]
        len_two = len(rec_two.seq)
        name_two = rec_two.id

        for jdx in range(idx, len(sequences)):
            rec_one = seq_dict[sequences[jdx]]
            len_one = len(rec_one.seq)
            name_one = rec_one.id

            counter += 1
            seq_counter += 1

            logging.debug(sequences[jdx])
            seq_text += ' ' + sequences[jdx]

            if not seq_counter % 25:
                # print(seq_counter)
                log_txt += str(seq_counter)

            # optional shade background according to length of LCS and/or user matrix
            #########################################################################

            # get interval based on LCS
            background_colors = [None, None]
            if lcs_shading and (
                lcs_shading_ref == 1 or lcs_shading_ref == 2 or max_lcs is not None
            ):  # self plot max_lcs_for == None
                lcs_len = lcs_dict[(idx, jdx)]
                l1 = lcs_len[0]  # forward
                l2 = lcs_len[1]  # reverse complement

                lcs_shading_bool = True

                # calculate shading acc. to chosen option
                if lcs_shading_ref == 1:  # percentage of shorter sequence
                    color_idx0 = min(
                        len(colors) - 1, l1 * lcs_shading_num // min(len_one, len_two)
                    )
                    color_idx1 = min(
                        len(colors) - 1, l2 * lcs_shading_num // min(len_one, len_two)
                    )
                elif lcs_shading_ref == 2:  # by given interval size
                    color_idx0 = min(len(colors) - 1, l1 // lcs_shading_interval_len)
                    color_idx1 = min(len(colors) - 1, l2 // lcs_shading_interval_len)
                    if color_idx0 >= len(colors):
                        color_idx0 = len(colors)
                    if color_idx1 >= len(colors):
                        color_idx1 = len(colors)
                else:  # percentage of maximum lcs length
                    color_idx0 = min(len(colors) - 1, l1 * lcs_shading_num // max_lcs)
                    color_idx1 = min(len(colors) - 1, l2 * lcs_shading_num // max_lcs)
            else:
                lcs_shading_bool = False

            # get interval based on custom matrix
            if custom_shading:
                # matrix value
                try:
                    custom_value = custom_similarity_dict[(idx, jdx)]
                except Exception:
                    custom_value = ''

                # bottom left triangle = LCS forward/reverse or best of both
                if lcs_shading_bool:
                    if lcs_shading_ori == 0:  # forward
                        color_idx1 = color_idx0
                    elif lcs_shading_ori == 2:  # both directions
                        color_idx1 = max(color_idx0, color_idx1)

                # top right triangle = custom value (not colored if text matrix provided)
                if type(custom_value) is int or type(custom_value) is float:
                    color_idx0 = int(
                        (custom_value - custom_min)
                        * lcs_shading_num
                        // (custom_max - custom_min)
                    )
                # no color if string is proviced
                else:
                    color_idx0 = 0

            # use best LCS of both orientations for coloring triangle with two-ori-LCS
            if (
                representation != 0 and lcs_shading_ori == 2
            ):  # both directions in triangle
                color_idx0, color_idx1 = (
                    max(color_idx0, color_idx1),
                    max(color_idx0, color_idx1),
                )

            # set colors dependent on lcs dependent on orientation
            if lcs_shading_bool and not custom_shading:
                if idx != jdx:
                    if lcs_shading_ori == 0:
                        color_idx1 = color_idx0
                    elif lcs_shading_ori == 1:
                        color_idx0 = color_idx1
                    background_colors[0] = colors[color_idx0]
                    background_colors[1] = colors[color_idx1]
                # for selfcomparison, only color reverse complement
                elif lcs_shading_ori != 0 and not custom_shading:
                    background_colors[0] = colors[color_idx1]
            # set different colors for shading by LCS + user matrix
            elif lcs_shading_bool and custom_shading:
                # print colors, background_colors, color_idx0, color_idx1
                background_colors[0] = colors_2[color_idx0]
                background_colors[1] = colors[color_idx1]
            # set grey color range for user matrix if no LCS shading
            elif custom_shading:
                background_colors[0] = colors[color_idx0]
                background_colors[1] = colors[color_idx0]

            if custom_shading and lcs_shading_bool:
                lcs_text += (
                    '\t'.join(
                        [
                            name_one,
                            name_two,
                            str(lcs_len[0]),
                            str(lcs_len[1]),
                            str(custom_value),
                            str(color_idx0),
                            str(color_idx1),
                        ]
                    )
                    + '\n'
                )
            elif lcs_shading_bool:
                lcs_text += (
                    '\t'.join(
                        [
                            name_one,
                            name_two,
                            str(lcs_len[0]),
                            str(lcs_len[1]),
                            str(color_idx0),
                            str(color_idx1),
                        ]
                    )
                    + '\n'
                )
            elif custom_shading:
                lcs_text += (
                    '\t'.join(
                        [
                            name_one,
                            name_two,
                            str(custom_value),
                            str(color_idx0),
                            str(color_idx1),
                        ]
                    )
                    + '\n'
                )

            # calculate figure position in polyplot
            # diagonal (self-dotplots)
            if idx == jdx:
                if mirror_y_axis:
                    seq_num = sequences.index(name_one) + 1
                    counter1 = seq_num + len(sequences) * (len(sequences) - seq_num)
                    counter = counter + (counter - 1) // (nrows)
                else:
                    # skip positions below diagonal
                    counter1 = counter + (counter - 1) // (nrows)  # + row_pos
                    counter = counter1
                counters = [counter1]

            # draw both graphs at once (due to symmetry)
            else:
                if mirror_y_axis:
                    col_pos = sequences.index(name_two) + 1
                    row_pos = len(sequences) - (sequences.index(name_one) + 1)
                    counter1 = row_pos * ncols + col_pos
                    counter2 = (ncols - col_pos) * ncols + ncols - row_pos
                else:
                    counter1 = counter
                    col_pos = (counter - 1) % ncols
                    row_pos = (counter - 1) // (nrows)
                    counter2 = col_pos * ncols + row_pos + 1
                counters = [counter1, counter2]  # lower, upper

            if len(counters) == 2:
                seq_counter += 1
                if not seq_counter % 25:
                    # print(seq_counter)
                    log_txt += str(seq_counter)

            x_lists, y_lists, x_lists_rc, y_lists_rc = data_dict[(idx, jdx)]

            # plot diagram(s)
            for kdx in range(len(counters)):
                if (
                    representation == 0
                    or len(counters) == 1
                    or (representation == 1 and kdx == 0)
                    or (representation == 2 and kdx == 1)
                ):
                    fig_pos = counters[kdx]
                    # plotting subplot with matplotlib
                    ax = P.subplot(gs[fig_pos - 1])  # rows, columns, plotnumber

                    # shade annotated regions if gff file(s) provided
                    if idx == jdx and gff_files:
                        if name_one in list(feat_dict.keys()):
                            features = feat_dict[name_one]
                            if len_two != len_one:
                                logging.info(
                                    'Polydot GFF shading for diagonal fields - nequal length error!'
                                )
                                return
                            for item in features:
                                feat_type, start, stop = item
                                feat_color, strength, zoom = gff_color_dict[
                                    feat_type.lower()
                                ]
                                start = max(0, start - zoom - 0.5)
                                stop = min(len_one + 1, stop + zoom + 0.5)
                                width = stop - start
                                ax.add_patch(
                                    patches.Rectangle(
                                        (start, start),  # (x,y)
                                        width,
                                        width,  # width, height
                                        edgecolor=None,
                                        linewidth=line_width + zoom,
                                        fill=True,
                                        facecolor=feat_color,
                                        alpha=strength,
                                    )
                                )

                    # if custom matrix value printed into upper matrix triangle, skip data plotting
                    # text print in top triangle
                    if user_matrix_print and custom_shading and kdx == 0 and idx != jdx:
                        data_plotting = False
                    # dotplot in bottom triangle
                    else:
                        data_plotting = True

                    # mirror plot, if plotting below diagonal
                    if kdx == 0:
                        l1, l2 = len_one, len_two
                        n1, n2 = name_one, name_two
                        x1, y1 = x_lists, y_lists
                        x2, y2 = x_lists_rc, y_lists_rc
                    else:
                        l2, l1 = len_one, len_two
                        n2, n1 = name_one, name_two
                        x1, y1 = y_lists, x_lists
                        x2, y2 = y_lists_rc, x_lists_rc

                    if mirror_y_axis:
                        x1, y1, x2, y2 = y1, x1, y2, x2
                        n1, n2 = n2, n1

                    if data_plotting:
                        # collect lines
                        lines = []
                        color_list = []
                        for x_lines, y_lines, col in [
                            (x2, y2, line_col_rev),
                            (x1, y1, line_col_for),
                        ]:
                            # If color is not white, add lines to plot
                            if col != 'white':
                                for ldx in range(len(x_lines)):
                                    lines.append(
                                        [
                                            (x_lines[ldx][0], y_lines[ldx][0]),
                                            (x_lines[ldx][-1], y_lines[ldx][-1]),
                                        ]
                                    )
                                    color_list.append(col)
                        color_list = np.array(color_list)

                        # draw lines
                        lc = cllct.LineCollection(
                            lines, colors=color_list, linewidths=line_width
                        )
                        ax.add_collection(lc)

                    # plot value provided by customer instead of dotplot
                    else:
                        alignment = {
                            'horizontalalignment': 'center',
                            'verticalalignment': 'center',
                        }
                        # P.text(0.5, 0.5, custom_value, size='medium', transform=ax.transAxes, **alignment)
                        P.text(
                            0.5,
                            0.5,
                            custom_value,
                            size=label_size * 1.5,
                            transform=ax.transAxes,
                            **alignment,
                        )
                        # P.text(0.5, 0.5, custom_value, size=label_size*1.5, transform=ax.transAxes,
                        # horizontalalignment='center', verticalalignment='center', color="black")

                    if custom_shading:
                        # omit diagonal
                        if idx == jdx:
                            ax.set_facecolor('white')
                        # use white background for text fields (top right triangle only [kdx 0])
                        elif (
                            type(custom_value) is not int
                            and type(custom_value) is not float
                            and kdx == 0
                        ):
                            ax.set_facecolor('white')
                        else:
                            ax.set_facecolor(background_colors[kdx])
                    # set background color if lcs shading
                    elif lcs_shading_bool and background_colors[kdx]:
                        ax.set_facecolor(background_colors[kdx])

                    # set axis limits
                    # P.xlim(0, l1+1)
                    if mirror_y_axis:
                        P.xlim(0, l2 + 1)
                        P.ylim(0, l1 + 1)  # rotate y axis (point upwards)
                    else:
                        P.xlim(0, l1 + 1)
                        P.ylim(l2 + 1, 0)  # rotate y axis (point downwards)

                    ## axis labelling
                    ##################

                    # determine axis positions
                    if x_label_pos_top:
                        ax.xaxis.tick_top()
                        ax.xaxis.set_label_position('top')
                        x_label_bool = fig_pos <= ncols
                        x_tick_bool = fig_pos > ncols * (ncols - 1)
                    else:
                        x_label_bool = fig_pos > ncols * (ncols - 1)
                        x_tick_bool = fig_pos <= ncols

                    # settings for y labels on right side
                    if y_label_pos == 0:  # right label
                        ax.yaxis.tick_right()
                        ax.yaxis.set_label_position('right')
                        label_dist = 30
                    else:
                        label_dist = 8

                    # x axis labels dependent on plot position/number
                    if x_label_bool:  # x title and labels on top or bottom
                        P.xlabel(
                            unicode_name(
                                shorten_name(
                                    n1,
                                    max_len=title_length,
                                    title_clip_pos=title_clip_pos,
                                )
                            ),
                            fontsize=label_size,
                            rotation=x_label_rotation,
                            verticalalignment=xvalign,
                            horizontalalignment=xhalign,
                            fontweight='bold',
                            labelpad=8,
                        )  # axis naming
                        if x_label_rotation not in ['horizontal', 'vertical']:
                            P.setp(
                                ax.get_xticklabels(),
                                fontsize=label_size * 0.9,
                                rotation='vertical',
                            )
                        else:
                            P.setp(
                                ax.get_xticklabels(),
                                fontsize=label_size * 0.9,
                                rotation=x_label_rotation,
                            )
                    elif x_tick_bool and x_label_pos_top:  # x ticks on bottom row
                        ax.xaxis.tick_bottom()  # ticks without labels on bottom
                        P.setp(
                            ax.get_xticklabels(),
                            fontsize=label_size,
                            rotation=x_label_rotation,
                            visible=False,
                        )
                    elif x_tick_bool:  # x ticks on top row
                        ax.xaxis.tick_top()  # # ticks without labels on top
                        P.setp(
                            ax.get_xticklabels(),
                            fontsize=label_size,
                            rotation=x_label_rotation,
                            visible=False,
                        )  # inner diagrams without labelling
                    elif idx == jdx and representation != 0:
                        if not mirror_y_axis and representation == 1:  # upper
                            ax.xaxis.tick_bottom()
                        elif mirror_y_axis and representation == 2:  # lower
                            ax.xaxis.tick_top()
                        elif mirror_y_axis and representation == 1:  # upper
                            ax.xaxis.tick_bottom()
                        elif not mirror_y_axis and representation == 2:  # lower
                            ax.xaxis.tick_top()
                        P.setp(
                            ax.get_xticklabels(), visible=False
                        )  # inner diagrams without labelling
                    else:  # no x ticks on internal rows
                        ax.axes.get_xaxis().set_visible(False)

                    # y axis labels dependent on plot position/number
                    if fig_pos % ncols == y_label_pos or (
                        ncols == 1 and nrows == 1
                    ):  # y title and labels in 1st column
                        P.ylabel(
                            unicode_name(
                                shorten_name(
                                    n2,
                                    max_len=title_length,
                                    title_clip_pos=title_clip_pos,
                                )
                            ),
                            fontsize=label_size,
                            rotation=y_label_rotation,
                            verticalalignment=yvalign,
                            horizontalalignment=yhalign,
                            fontweight='bold',
                            labelpad=label_dist,
                        )
                        P.setp(
                            ax.get_yticklabels(), fontsize=label_size * 0.9
                        )  # axis naming
                    elif fig_pos % ncols == 0:  # y ticks in last column
                        ax.yaxis.tick_right()
                        P.setp(
                            ax.get_yticklabels(), visible=False
                        )  # inner diagrams without labelling
                    elif idx == jdx and representation != 0:
                        if not mirror_y_axis and representation == 1:  # upper
                            ax.yaxis.tick_left()
                        elif mirror_y_axis and representation == 2:  # lower
                            ax.yaxis.tick_left()
                        elif mirror_y_axis and representation == 1:  # upper
                            ax.yaxis.tick_right()
                        elif not mirror_y_axis and representation == 2:  # lower
                            ax.yaxis.tick_right()
                        P.setp(
                            ax.get_yticklabels(), visible=False
                        )  # inner diagrams without labelling
                    else:
                        ax.axes.get_yaxis().set_visible(False)

    log_txt += '\n%d done' % seq_counter
    logging.info(log_txt)

    try:
        logging.debug(lcs_text)
    except Exception:
        pass

    # finalize layout - margins & spacing between plots
    P.tick_params(axis='both', which='major', labelsize=label_size * 0.9)
    try:
        P.tight_layout(h_pad=0.02, w_pad=0.02)
    except Exception:
        logging.info(
            'Attention - pylab.tight_layout failed! Please check sequence names and layout settings!'
        )
    # gs.tight_layout(fig, h_pad=.02, w_pad=.02) # less overlapping tick labels, but also disturbingly large spacing
    if y_label_rotation == 'horizontal':
        if x_label_pos_top:
            P.subplots_adjust(
                hspace=spacing, wspace=spacing, left=0.13, top=0.87
            )  # space between rows - def 0.4
        else:
            P.subplots_adjust(
                hspace=spacing, wspace=spacing, left=0.13, bottom=0.13
            )  # space between rows - def 0.4
    else:
        P.subplots_adjust(
            hspace=spacing, wspace=spacing
        )  # space between rows - def 0.4

    # save figure and close instance
    fig_name = '%s%s%s%s.%s' % (prefix, name_graph, wordsize_suffix, suffix, filetype)
    P.savefig(fig_name)
    P.close()
    P.cla()

    # create figure color legend
    if lcs_shading:
        if lcs_shading_ref == 1:  # percentage of shorter sequence
            legend_file_name = legend_figure(
                colors, lcs_shading_num, unit='%', filetype=filetype, prefix=prefix
            )
        elif lcs_shading_ref == 2:  # interval sizes
            legend_file_name = legend_figure(
                colors,
                lcs_shading_num,
                unit=aa_bp_unit,
                filetype=filetype,
                prefix=prefix,
                bins=color_bins,
            )
        else:  # relative of maximum lcs
            legend_file_name = legend_figure(
                colors,
                lcs_shading_num,
                unit=aa_bp_unit,
                filetype=filetype,
                prefix=prefix,
                max_len=max_lcs,
            )

    if custom_shading:
        custom_prefix = 'custom-matrix-' + prefix
        legend_file_name_custom = legend_figure(
            colors_2,
            lcs_shading_num,
            unit='%',
            filetype=filetype,
            prefix=custom_prefix,
            max_len=custom_max,
            min_len=custom_min,
        )

    if lcs_shading and custom_shading:
        return [fig_name, legend_file_name, legend_file_name_custom]
    elif lcs_shading:
        return [fig_name, legend_file_name]
    elif custom_shading:
        return [fig_name, legend_file_name_custom]
    else:
        return [fig_name]
