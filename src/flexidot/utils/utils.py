import logging
import re
import time

from colormap import rgb2hex
from colour import Color
from matplotlib import colormaps


# TODO: Remove internal logging message. Return delta time instead of printing it.
# Check if "now" value is used elsewhere.
def time_track(starting_time, show=True):
    """
    calculate time passed since last time measurement
    """
    now = time.time()
    delta = now - starting_time
    if show:
        logging.info(f'{delta} seconds')
    return delta


def calc_fig_ratio(ncols, nrows, plot_size):
    """
    Calculate size ratio for given number of columns (ncols) and rows (nrows)
    with plot_size as maximum width and length
    """
    ratio = ncols * 1.0 / nrows
    logging.debug(' '.join([str(ncols), str(nrows), str(ratio)]))
    if ncols >= nrows:
        figsize_x = plot_size
        figsize_y = plot_size / ratio
    else:
        figsize_x = plot_size * ratio
        figsize_y = plot_size
    return figsize_x, figsize_y


def shorten_name(seq_name, max_len=20, title_clip_pos='B'):  # , delim="_"):
    """
    shorten sequence names (for diagram titles)
    """

    if len(seq_name) <= max_len:
        return seq_name

    # take last characters
    if title_clip_pos == 'E':
        name = seq_name[len(seq_name) - max_len :]

    # take first characters
    else:
        name = seq_name[:max_len]

    """# keep first and last part if multiple parts separated by delimiter (e.g. species_prefix + sequence_id)
    if delim in seq_name:
        if seq_name.count(delim) >= 2:
            name = "%s..." % delim.join(seq_name.split(delim)[:1]) + seq_name.split(delim)[-1] # .replace("_000", "-")
        else:
            name = seq_name[:((max_len-2)//2)] + "..." + seq_name[((max_len-2)//2):]

        if len(name) > max_len:
            name = name[:((max_len-2)//2)] + "..." + name[((max_len-2)//2):]
    else:
        name = seq_name[:((max_len-2)//2)] + "..." + seq_name[((max_len-2)//2):]
    """

    return name


def sanitize_filename(name):
    """
    replace characters that are unsafe in filenames (e.g. path separators
    found in sequence names such as "TE#Class/Family") with underscores
    """
    return re.sub(r'[\\/]', '_', name)


def unicode_name(name):
    """
    replace non-ascii characters in string (e.g. for use in matplotlib)
    """
    unicode_string = eval('u"%s"' % name)
    return ''.join(char for char in unicode_string if ord(char) < 128)
    # return unicodedata.normalize("NFKD", unicode_string).encode("ascii", "ignore")


def create_color_list(
    number: int, color_map: str = 'Greys', max_grey: str = '#595959'
) -> list:
    """
    Create color list with given number of entries
    grey by default, matplotlib color_map can be provided
    """
    if color_map not in list(colormaps):
        logging.warning(
            f'Invalid color_map {color_map} provided! - Examples: {list(colormaps)}.'
        )
        logging.warning('See https://matplotlib.org/users/colormaps.html\n')

    try:
        # create pylab colormap
        cmap = eval('P.cm.' + color_map)

        # get descrete color list from pylab
        cmaplist = [cmap(i) for i in range(cmap.N)]  # extract colors from map

        # determine positions for number of colors required
        steps = round((len(cmaplist) - 1) / (number))
        numbers = list(range(0, len(cmaplist), int(steps)))

        # extract RGB color and convert to hex code
        colors = []
        for idx in numbers[:-1]:
            rgba_color = cmaplist[idx]
            rgb_color = rgba_color[:3]
            col = rgb2hex(
                int(rgb_color[0] * 255),
                int(rgb_color[1] * 255),
                int(rgb_color[2] * 255),
            )
            colors.append(col)

    # Default to Greys color scheme if an error occurs
    except Exception as e:
        logging.warning(f'An error occurred: {e}')
        logging.warning('Using grey color scheme instead.')

        colors = list(Color('#FFFFFF').range_to(Color(max_grey), number))  # grey
        for idx in range(len(colors)):
            colors[idx] = str(colors[idx]).replace('Color ', '')
            if '#' in colors[idx] and len(colors[idx]) != 7:
                # print colors[idx]
                colors[idx] = colors[idx] + colors[idx][-(7 - len(colors[idx])) :]

    logging.info('%d Colors: %s' % (len(colors), ', '.join(colors)))

    if len(colors) < number:
        logging.info(
            '\nError in color range definition! %d colors missing\n'
            % (number - len(colors))
        )

    return colors
