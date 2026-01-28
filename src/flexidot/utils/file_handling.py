import logging
import os
from typing import Dict, List, Tuple, Union

from Bio import SeqIO
import matplotlib.colors as mcolors
import pylab as P


# Default color palette for auto-generated GFF colors
_DEFAULT_GFF_COLORS = [
    '#e41a1c',  # red
    '#377eb8',  # blue
    '#4daf4a',  # green
    '#984ea3',  # purple
    '#ff7f00',  # orange
    '#ffff33',  # yellow
    '#a65628',  # brown
    '#f781bf',  # pink
    '#999999',  # grey
    '#66c2a5',  # teal
    '#fc8d62',  # coral
    '#8da0cb',  # light blue
    '#e78ac3',  # light pink
    '#a6d854',  # lime
    '#ffd92f',  # gold
    '#e5c494',  # tan
    '#b3b3b3',  # silver
]


def extract_gff_annotation_types(
    gff_files: Union[str, List[str]],
) -> Dict[str, Tuple[str, float, int]]:
    """
    Extract unique annotation types from GFF3 files and generate auto-colors.

    When no GFF color config file is provided, this function extracts all unique
    annotation types (column 3) from the GFF files and assigns a unique color
    to each type.

    Parameters
    ----------
    gff_files : str or list of str
        Path(s) to GFF3 file(s).

    Returns
    -------
    dict
        Dictionary with annotation type as key (lowercase) and tuple
        (color, alpha, zoom) as value.

    Examples
    --------
    >>> colors = extract_gff_annotation_types(['annotations.gff3'])
    >>> print(colors)
    {'micro': ('#e41a1c', 0.7, 0), 'ssatar1': ('#377eb8', 0.7, 0), ...}
    """
    if isinstance(gff_files, str):
        gff_files = [gff_files]

    # Extract unique annotation types from all GFF files
    annotation_types = set()
    for gff_file in gff_files:
        if not os.path.exists(gff_file):
            logging.warning(f'GFF file not found: {gff_file}')
            continue

        with open(gff_file, 'r') as f:
            for line in f:
                if not line.startswith('#') and line.strip():
                    fields = line.strip().split('\t')
                    if len(fields) >= 3:
                        feat_type = fields[2].lower()
                        annotation_types.add(feat_type)

    # Generate color dictionary with unique colors for each annotation type
    gff_feat_colors = {}
    sorted_types = sorted(annotation_types)
    for idx, feat_type in enumerate(sorted_types):
        color = _DEFAULT_GFF_COLORS[idx % len(_DEFAULT_GFF_COLORS)]
        gff_feat_colors[feat_type] = (color, 0.7, 0)

    # Always include 'others' for unknown annotation types
    if 'others' not in gff_feat_colors:
        gff_feat_colors['others'] = ('grey', 0.5, 0)

    # Log the auto-generated colors
    if gff_feat_colors:
        logging.warning(
            'No GFF color config file provided. Auto-generating colors for annotation types.\n'
            'Consider using --gff_color_config for custom colors.\n'
            'Auto-generated colors:'
        )
        for feat_type, (color, alpha, zoom) in sorted(gff_feat_colors.items()):
            logging.info(f'  {feat_type}: {color} (alpha={alpha}, zoom={zoom})')

    return gff_feat_colors


def read_seq(input_fasta, degap=False):
    """
    Read fasta sequences from (all) file(s)
    """

    # Check if file provided
    if input_fasta == [] or input_fasta == '':
        raise ValueError('No file names provided: %s' % input_fasta)

    # Initialize variables
    concat_created = False

    # Combine sequence files, if required
    if type(input_fasta) is list:
        # Concatenate fasta files
        if len(input_fasta) > 1:
            logging.info(
                'Concatenating sequences from multiple files: {}'.format(input_fasta)
            )
            input_fasta_combi = concatenate_files(input_fasta)
            concat_created = True
        else:
            input_fasta_combi = input_fasta[0]
    else:
        # Single fasta file
        input_fasta_combi = input_fasta

    # read sequences
    logging.info(f'Reading sequences from {input_fasta_combi}')
    try:
        # SeqIO.index returns a dictionary with the sequence id as key and the sequence as value
        seq_dict = SeqIO.index(input_fasta_combi, 'fasta')
    except ValueError as e:
        logging.error(
            f'ValueError: {e} - please check input files, e.g. for duplicate names!'
        )
        return {}, []
    except FileNotFoundError as e:
        logging.error(f'FileNotFoundError: {e} - please check if the file exists!')
        return {}, []
    except Exception as e:
        logging.error(f'Unexpected error: {e} - please check input files!')
        return {}, []

    for seq_id in seq_dict:
        if '-' in seq_dict[seq_id].seq:
            logging.warning('Gaps detected in sequence: %s' % seq_id)
            return read_seq(degap_fasta(input_fasta), degap=True)

    # Get sequence names for sorting
    sequences = []
    for item in SeqIO.parse(input_fasta_combi, 'fasta'):
        sequences.append(item.id)

    # If degap=True remove input file after processing
    if degap and type(input_fasta) is list:
        logging.info('Removing temp degapped input fasta files: {}'.format(input_fasta))
        for item in input_fasta:
            os.remove(item)
    elif degap:
        logging.info('Removing temp degapped input fasta file: {}'.format(input_fasta))
        os.remove(input_fasta)

    # If concatenation was required, remove combined file
    if concat_created:
        logging.info('Removing concatenated tempfile: {}'.format(input_fasta_combi))
        os.remove(input_fasta_combi)

    return seq_dict, sequences


def read_gff_color_config(gff_color_config_file=''):
    """
    define coloring options for gff-based color shading of self-dotplots
    """

    # default aestetics for annotation shading (e.g. if no user config file is provided)
    # dictionary with feature_type as key and tuple(color, transparency, zoom) as value
    gff_feat_colors = {
        'orf': ('#b41a31', 0.2, 0),
        'orf_rev': ('#ff773b', 0.3, 0),
        'gene': ('#b41a31', 0.2, 0),
        'cds': ('darkorange', 0.2, 0),
        'exon': ('orange', 0.2, 0),
        'intron': ('lightgrey', 0.2, 0),
        'utr': ('lightblue', 0.2, 0),
        'repeat_region': ('green', 0.3, 0),
        'repeat': ('green', 0.3, 0),
        'tandem_repeat': ('red', 0.3, 0),
        'transposable_element': ('blue', 0.3, 0),
        'ltr_retrotransposon': ('#cccccc', 0.5, 0),
        'ltr-retro': ('#cccccc', 0.5, 0),
        'long_terminal_repeat': ('#2dd0f0', 0.75, 2),
        'ltr': ('#2dd0f0', 0.75, 2),
        'pbs': ('purple', 0.75, 2),
        'ppt': ('#17805a', 0.5, 2),
        'target_site_duplication': ('red', 0.75, 2),
        'misc_feature': ('grey', 0.3, 0),
        'misc_feat': ('grey', 0.3, 0),
        'misc': ('grey', 0.3, 0),
        'others': ('grey', 0.5, 0),
    }
    if not gff_color_config_file or not os.path.exists(str(gff_color_config_file)):
        logging.info('No custom GFF color configuration found. Using defaults.')
        return gff_feat_colors

    logging.info('Updating GFF color configuration with custom specifications.')

    # Read custom gff_color_config_file
    in_file = open(gff_color_config_file, 'r')
    overwritten = set()
    for line in in_file:
        if len(line.strip().split('\t')) >= 4 and not line.startswith('#'):
            data = line.strip().split('\t')
            feat = data[0].lower()
            color = data[1].lower()

            # Check, if settings are valid
            if not mcolors.is_color_like(color):
                color = 'grey'
                text = 'Invalid color specified for %s: %s - default grey' % (
                    data[0],
                    data[1],
                )
                logging.info(text)
            try:
                alpha = float(data[2])
            except ValueError:
                alpha = 0.75
                text = 'Invalid alpha specified for %s: %s - default 0.75' % (
                    data[0],
                    data[2],
                )
                logging.info(text)
            try:
                zoom = float(data[3])
            except Exception:
                zoom = 0
                text = 'Invalid zoom specified for %s: %s - default 0' % (
                    data[0],
                    data[3],
                )
                logging.info(text)

            # Track changes of predefined settings
            if feat in list(gff_feat_colors.keys()):
                overwritten.add(data[0].lower())

            gff_feat_colors[feat] = (color, alpha, zoom)

    in_file.close()

    # Default coloring for unknown annotations
    if 'others' not in list(gff_feat_colors.keys()):
        gff_feat_colors['others'] = ('grey', 0.5, 0)

    # Print configuration
    text = '\n\nGFF color specification:\n%s\n' % (60 * '.')
    for item in sorted(gff_feat_colors.keys()):
        text += '%-30s\t%-10s\t%-5s\t%s\n' % (
            item,
            str(gff_feat_colors[item][0]),
            str(gff_feat_colors[item][1]),
            str(gff_feat_colors[item][2]),
        )
    logging.debug(text)

    # print overwritting feature type specifications
    if len(overwritten) != 0:
        text = '%d feature type specifications overwritten:' % len(overwritten)
        text += '\n\t' + ', '.join(overwritten) + '\n'
        logging.info(text)

    text = 'GFF color specification updated acc. to %s\n\t%s\n\n' % (
        gff_color_config_file,
        ', '.join(gff_feat_colors),
    )
    logging.info(text)

    return gff_feat_colors


def read_gffs(
    input_gff_files,
    color_dict=None,
    type_nuc=True,
    prefix='',
    filetype='png',
):
    """
    Create feature dictionary from input_gff
    sequence name as key and (feature type, start, stop) as value
    """
    if color_dict is None:
        color_dict = {'others': ('grey', 1, 0)}

    if type(input_gff_files) is not list:
        input_gff_files = [input_gff_files]

    # Create dictionary with seq_name as key and (type, start and stop) as value
    unknown_feats = set()
    used_feats = set()
    feat_dict = {}
    for input_gff in input_gff_files:
        text = '...reading ' + input_gff
        logging.info(text)

        in_file = open(input_gff, 'r')
        for line in in_file:
            if not line.startswith('#') and line.strip() != '':
                data = line.strip().split('\t')
                feat_type = data[2].lower()
                if data[6] == '-':
                    feat_type += '_rev'
                if feat_type.lower() not in list(color_dict.keys()):
                    if feat_type.lower().replace('_rev', '') in list(color_dict.keys()):
                        feat_type = feat_type.replace('_rev', '')
                    else:
                        unknown_feats.add(feat_type)
                        feat_type = 'others'
                used_feats.add(feat_type)
                if data[0] not in list(feat_dict.keys()):
                    feat_dict[data[0]] = [
                        (feat_type, int(data[3]), int(data[4]))
                    ]  # feature type, start, stop
                else:
                    feat_dict[data[0]].append(
                        (feat_type, int(data[3]), int(data[4]))
                    )  # feature type, start, stop

        text = '\nAnnotations for: %s\n' % ', '.join(list(feat_dict.keys())[:10])
        if len(list(feat_dict.keys())) > 10:
            text = text[:-1] + ', ...\n'
        logging.info(text)

        in_file.close()

    # print feature types without specific shading settings
    if len(unknown_feats) != 0:
        text = 'Missing shading specification for %d feature type(s):\n\t%s\n' % (
            len(unknown_feats),
            ', '.join(sorted(unknown_feats)),
        )
        logging.info(text)

    # create color legend
    colors, alphas = [], []
    for item in sorted(used_feats):
        colors.append(color_dict[item][0])
        alphas.append(color_dict[item][1])
    legend_figure(
        colors=colors,
        lcs_shading_num=len(used_feats),
        type_nuc=type_nuc,
        bins=sorted(used_feats),
        alphas=alphas,
        gff_legend=True,
        prefix=prefix,
        filetype=filetype,
    )

    # print settings
    text = 'GFF Feature Types: %s\nGFF Colors:        %s' % (
        ', '.join(sorted(used_feats)),
        ', '.join(sorted(colors)),
    )
    logging.info(text)

    return feat_dict


def read_matrix(matrix_file_name, delim='\t', symmetric=True, recursion=False):
    input_file = open(matrix_file_name, 'r')

    # read sequence names from first column
    names = []
    for line in input_file:
        if not line.startswith('#') and not line.startswith(delim) and delim in line:
            names.append(line.strip().split(delim)[0])
    logging.info(
        "Delimiter '%s': %d names - %s\n" % (delim, len(names), ', '.join(names))
    )

    # check if names were found - otherwise try another delimiter
    if names == [] and not recursion:
        if delim == '\t':
            new_delim = ','
        else:
            new_delim = '\t'
        logging.info(
            "\nMatrix file not containing data delimited by '%s' - trying to read matrix with delimiter '%s'"
            % (delim.replace('\t', '\\t'), new_delim)
        )
        info_dict = read_matrix(
            matrix_file_name,
            delim=new_delim,
            symmetric=symmetric,
            recursion=True,
        )
        return info_dict
    elif names == []:
        logging.info('Empty matrix file with alternative delimiter!')
        return info_dict
    input_file.close()

    input_file = open(matrix_file_name, 'r')
    # read matrix entries as values in dictionary with tuple(names) as key
    info_dict = {}
    contradictory_entries = []
    for line in input_file:
        if not line.startswith('#') and not line.startswith(delim) and delim in line:
            data = line.strip().split(delim)
            for idx in range(len(data[1:])):
                # print tuple(sorted([data[0], names[idx]])), data[idx+1]
                if symmetric:
                    key = tuple(sorted([names[idx], data[0]]))
                else:
                    key = tuple(names[idx], data[0])
                if key in list(info_dict.keys()):
                    if (
                        symmetric
                        and info_dict[key] != data[idx + 1]
                        and data[idx + 1] not in ['', '-']
                        and info_dict[key] not in ['', '-']
                    ):
                        contradictory_entries.append(key)
                info_dict[key] = data[idx + 1]
    input_file.close()

    if len(contradictory_entries) != 0:
        try:
            logging.info(
                '\nContradictory entries in matrix file %s:\n\t%s'
                % (matrix_file_name, ', '.join(contradictory_entries))
            )
        except Exception:
            log_txt = '\nContradictory entries in matrix file %s:\n\t' % (
                matrix_file_name
            )
            for item in contradictory_entries:
                log_txt += str(item).replace("'", '') + ', '
            log_txt = log_txt[:-2]
            logging.info(log_txt)
        logging.info('Using value from bottom left triangle!')

    logging.debug(f'Matrix information: {", ".join(names)}')

    return info_dict


def concatenate_files(file_list, combi_filename='temp_combined.fasta'):
    """
    concatenate content of all files in file_list into a combined file named combi_filename
    """
    out_file = open(combi_filename, 'w')
    text = ''
    for item in file_list:
        text += item + ' '
        # read in_file linewise and write to out_file
        in_file = open(item, 'r')
        for line in in_file:
            out_file.write(line.strip() + '\n')
        in_file.close()
    out_file.close()

    logging.debug(text)

    return combi_filename


def degap_fasta(input_fasta):
    """
    Remove gaps from fasta - Write new degapped fasta.
    """
    logging.info('Removing gaps from fasta files: %s' % input_fasta)
    # degap all sequence files
    output_fastas = []
    if type(input_fasta) is not list:
        input_fasta = list(input_fasta)
    for input_fas in input_fasta:
        output_fas = input_fas[: input_fas.rfind('.')] + '_degapped.fas'
        in_file = open(input_fas, 'r')
        out_file = open(output_fas, 'w')
        for line in in_file:
            if line.startswith('>'):
                out_file.write(line.strip() + '\n')
            else:
                out_file.write(line.strip().replace('-', '') + '\n')
        out_file.close()
        in_file.close()
        output_fastas.append(output_fas)
    return output_fastas


def legend_figure(
    colors,
    lcs_shading_num,
    type_nuc=True,
    unit='%',
    filetype='png',
    max_len=None,
    min_len=0,
    bins=None,
    alphas=None,
    gff_legend=False,
    prefix='',
):
    """
    create figure color legend
    """
    max_legend_length_row = 8
    max_legend_length_col = 4

    # Handle default value for bins
    if bins is None:
        bins = []

    # define output file
    if filetype not in ['png', 'pdf', 'svg']:
        text = 'Provide valid file type - png, pdf, or svg'
        logging.info(text)
        filetype = 'png'

    # check if length of information fit
    if not gff_legend and (
        (bins != [] and len(colors) != lcs_shading_num + 1)
        or (bins != [] and len(colors) != len(bins) + 1)
    ):
        if bins != [] and len(colors) != lcs_shading_num + 1:
            text = (
                '**Attention**\nlcs_shading_num (%d) does not match number of colors (%d)!\n'
                % (lcs_shading_num, len(bins))
            )
        elif bins != [] and len(colors) != len(bins) + 1:
            text = (
                '**Attention**\nnumber of LCS length bins (%d) does not match number of colors (%d)!\n'
                % (len(colors), len(bins))
            )
        logging.info(text)
    elif gff_legend and len(bins) != len(colors):
        text = (
            '**Attention**\nnumber of GFF Feature Types (%d) does not match number of colors (%d)!\n'
            % (len(colors), len(bins))
        )
        logging.info(text)

    # set alpha values to opaque if none are provided
    if alphas is None:
        alphas = []
        alphas = [1] * len(colors)

    # legend data points
    data_points = list(range(len(colors)))
    if not gff_legend:
        # specify intervals, if max_len provided
        if max_len is not None:
            _multi_factor = 100  # one digit
            if max_len <= 1:
                _multi_factor = 1000  # two digits
            # len_interval_size = (max_len-min_len) * multi_factor *1. // lcs_shading_num * (1./ multi_factor)
            len_interval_size = (max_len - min_len) * 1.0 / lcs_shading_num
            len_pos = [float('%.2f' % (min_len))]
            # calculate interval positions
            len_pos.extend(
                [
                    float('%.2f' % (len_pos[-1] + len_interval_size))
                    for _ in range(lcs_shading_num)
                ]
            )

            if prefix.startswith('custom-matrix') and (
                0 <= max_len <= 100 and 0 <= min_len <= 100
            ):
                unit = '%'
            elif prefix.startswith('custom-matrix'):
                unit = ''

            text = (
                '\n%d Legend intervals from %.2f to %.2f: \n\t%s - number: %d, step: %.2f, unit: %s\n'
                % (
                    lcs_shading_num + 1,
                    min_len,
                    max_len,
                    str(len_pos),
                    len(len_pos),
                    len_interval_size,
                    unit,
                )
            )
            logging.info(text)
            pos = len_pos
            interval_size = len_interval_size
        # generate legend labels acc. to standard interval notation
        else:
            # use default max_len = 100 and min_len = 0
            len_interval_size = 100.0 / lcs_shading_num
            pos = [float('%.2f' % (0))]
            # calculate interval positions
            pos.extend(
                [
                    float('%.2f' % (pos[-1] + len_interval_size))
                    for _ in range(lcs_shading_num)
                ]
            )

            # interval_size = 100 // lcs_shading_num
            # pos = range(interval_size, 101+interval_size, interval_size)

        # remove unneccessary zeros in decimal places (i.e. if x.x00 in all entries)
        while True:
            last_digit_all_zero = True
            no_delim = False
            for idx in range(len(pos)):
                # only process if fraction with decimal places
                if '.' not in str(pos[idx]):
                    no_delim = True
                    break
                # only process when all entries end in zero
                elif str(pos[idx])[-1] != '0':
                    last_digit_all_zero = False
                    break
            if not last_digit_all_zero or no_delim:
                break
            # remove last decimal place (== 0) from all entries
            else:
                temp_pos = pos[:]
                for idx in range(len(pos)):
                    if not str(pos[idx])[-2] == '.':
                        pos[idx] = float(str(pos[idx])[:-1])
                    else:
                        pos[idx] = int(str(pos[idx])[:-2])
                logging.info('Shortening legend entries: %s - %s' % (temp_pos, pos))

        # eliminate fractions if unit == bp/aa
        if unit in ['aa', 'bp']:
            for idx in range(len(pos)):
                temp_pos = pos[:]
                rounded_unit = False
                if '.' in str(pos[idx]):
                    rounded_unit = True
                    # round values up to next integer (keep integer, if not a fraction)
                    pos[idx] = int(pos[idx] / 1) + int(pos[idx] % 1 > 0)
                    if idx == len(pos) - 1 and pos[idx] == 101:
                        pos[idx] = 100
            if rounded_unit:
                logging.info(
                    "Fractions not permitted for unit '%s': %s -> %s"
                    % (unit, temp_pos, pos)
                )

        if bins != []:  # labels provided
            legend_labels = bins[:]
            legend_labels.append('max')
            legend_labels_lengths = []
            for item in bins:
                legend_labels_lengths.append(
                    '[%d %s, %d %s)' % (item - min(bins), unit, item, unit)
                )
            if len(bins) == len(colors) - 1:
                legend_labels_lengths.append(
                    '[%d %s, %s]' % (max(bins), unit, '\u221e')
                )  # infinite

        else:
            legend_labels = []
            legend_labels_lengths = []
            for idx in range(len(pos)):
                num = pos[idx]
                try:
                    legend_labels.append('[%d%%, %d%%)' % (num - interval_size, num))
                except Exception:
                    legend_labels.append(
                        '[%d%%, %d%%)' % (num - len_interval_size, num)
                    )
                if max_len is not None:
                    num = len_pos[idx]
                    # as int or float
                    if num == int(num) and int(len_interval_size) == len_interval_size:
                        legend_labels_lengths.append(
                            '[%d %s, %d %s)'
                            % (num, unit, num + len_interval_size, unit)
                        )
                    else:
                        legend_labels_lengths.append(
                            '[%.2f %s, %.2f %s)'
                            % (num, unit, num + len_interval_size, unit)
                        )
            legend_labels[-1] = '100' + unit
            if max_len is not None:
                if num == int(num) and int(len_interval_size) == len_interval_size:
                    legend_labels_lengths[-1] = '[%d %s, \u221e]' % (max_len, unit)
                else:
                    legend_labels_lengths[-1] = '[%.2f %s, \u221e]' % (max_len, unit)

    # set labels and choose file name
    if gff_legend:
        label_text = bins[:]
        edge_col = None
        legend_file_name = 'GFF_Shading_Legend_n%d.' % lcs_shading_num + filetype
    elif max_len is not None:
        label_text = legend_labels_lengths[:]
        edge_col = 'black'
        legend_file_name = (
            'Polydotplot_LCS_Shading_Legend_max%d%s_n%d.'
            % (max_len, unit, lcs_shading_num + 1)
            + filetype
        )
    elif bins != []:
        label_text = legend_labels_lengths[:]
        edge_col = 'black'
        legend_file_name = (
            'Polydotplot_LCS_Shading_Legend_%d%s_n%d.'
            % (bins[0], unit, lcs_shading_num + 1)
            + filetype
        )
    else:
        label_text = legend_labels[:]
        edge_col = 'black'
        legend_file_name = (
            'Polydotplot_LCS_Shading_Legend_%%len_n%d.' % (lcs_shading_num + 1)
            + filetype
        )

    if prefix is not None and prefix != '':
        if not prefix.endswith('-'):
            prefix = prefix + '-'
        legend_type = 'LCS'
        if prefix.startswith('custom-matrix'):
            prefix = prefix.replace('custom-matrix', '')[1:]
            legend_type = 'CustomMatrix'
        legend_file_name = prefix + legend_file_name.replace('LCS', legend_type)

    # plot legend figure
    fig, ax = P.subplots(3, 1, figsize=(len(colors) * 2, len(colors) * 2))
    for idx in range(len(colors)):
        ax[0].bar(
            data_points[idx] + 1,
            data_points[idx] + 1,
            color=colors[idx],
            label=label_text[idx],
            alpha=alphas[idx],
            edgecolor=edge_col,
        )
        ax[1].bar(
            data_points[idx] + 1,
            0,
            color=colors[idx],
            label=label_text[idx],
            alpha=alphas[idx],
            edgecolor=edge_col,
        )
        ax[2].bar(
            data_points[idx] + 1,
            0,
            color=colors[idx],
            label=label_text[idx],
            alpha=alphas[idx],
            edgecolor=edge_col,
        )
    ax[1].set_ylim(0, 1)
    ax[2].set_ylim(0, 1)
    ax[1].legend(
        ncol=((len(colors) - 1) // max_legend_length_row) + 1, framealpha=1
    )  # vertical legend
    col_num = len(colors)
    if len(colors) > max_legend_length_col:
        remainder = 0
        if len(colors) % max_legend_length_col != 0:
            remainder = 1
        row_num = len(colors) // max_legend_length_col + remainder
        remainder = 0
        if len(colors) % row_num != 0:
            remainder = 1
        col_num = len(colors) // row_num + remainder
    ax[2].legend(ncol=col_num, framealpha=1)  # horizontal legend

    P.savefig(legend_file_name)

    return legend_file_name
