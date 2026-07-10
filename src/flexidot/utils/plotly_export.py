###################################
#   Interactive HTML Dot Plots   #
###################################
"""
Interactive HTML dotplot rendering using Plotly.

Covers the core dotplot view (match lines with zoom/pan/hover/click) plus
GFF annotation shading for the self and paired plotting modes. LCS shading,
custom matrix shading, and collage layouts are not represented in HTML
output; GFF shading is also not yet supported for the poly mode.
"""

import logging

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from flexidot.utils.utils import shorten_name, unicode_name

# solid box border around the plot area, similar to matplotlib's default axes spines
_AXIS_BOX_STYLE = dict(showline=True, linewidth=1, linecolor='black', mirror=True)

_DIV_ID = 'flexidot-dotplot'


def _extract_seq(seq, start, end):
    """Return the sequence substring spanning the (1-based, inclusive) start/end."""
    if seq is None:
        return ''
    lo, hi = sorted((int(start), int(end)))
    return str(seq[lo - 1 : hi])


def _line_traces(
    x_lists,
    y_lists,
    x_lists_rc,
    y_lists_rc,
    line_col_for,
    line_col_rev,
    line_width,
    aa_bp_unit='bp',
    seq_x=None,
    seq_y=None,
):
    """Build Plotly line traces for forward and reverse-complement matches.

    All segments of one color are combined into a single trace (separated by
    None) so plots with many matches stay responsive. Each segment carries
    customdata (x start/end, y start/end, length, aligned sequence text) used
    for the hover tooltip and the click-to-reveal sequence box.
    """
    traces = []
    for x_lines, y_lines, col, name in [
        (x_lists_rc, y_lists_rc, line_col_rev, 'reverse complement'),
        (x_lists, y_lists, line_col_for, 'forward'),
    ]:
        if col == 'white' or len(x_lines) == 0:
            continue
        xs, ys, customdata = [], [], []
        for ldx in range(len(x_lines)):
            x0, x1 = x_lines[ldx][0], x_lines[ldx][-1]
            y0, y1 = y_lines[ldx][0], y_lines[ldx][-1]
            length = abs(int(x1) - int(x0)) + 1
            seq_text = 'X %d-%d: %s\nY %d-%d: %s' % (
                x0,
                x1,
                _extract_seq(seq_x, x0, x1) or 'n/a',
                y0,
                y1,
                _extract_seq(seq_y, y0, y1) or 'n/a',
            )
            point = [x0, x1, y0, y1, length, seq_text]
            xs += [x0, x1, None]
            ys += [y0, y1, None]
            customdata += [point, point, point]
        traces.append(
            go.Scattergl(
                x=xs,
                y=ys,
                mode='lines',
                line=dict(color=col, width=line_width),
                name=name,
                customdata=customdata,
                hovertemplate=(
                    '<b>%s match</b><br>' % name.title()
                    + 'X: %{customdata[0]}-%{customdata[1]}<br>'
                    + 'Y: %{customdata[2]}-%{customdata[3]}<br>'
                    + 'Length: %{customdata[4]} '
                    + aa_bp_unit
                    + '<extra></extra>'
                ),
                showlegend=False,
            )
        )
    return traces


def _click_seq_postscript(div_id=_DIV_ID):
    """JS injected after the plot: shows the clicked match's aligned sequence in a box below it."""
    return """
var plotDiv = document.getElementById('%(div_id)s');
var seqBox = document.createElement('div');
seqBox.id = '%(div_id)s-seqbox';
seqBox.style.cssText = 'font-family: monospace; font-size: 13px; white-space: pre-wrap; word-break: break-all; padding: 10px 14px; margin-top: 10px; border: 1px solid #ccc; border-radius: 4px; background: #f7f7f7; min-height: 1.4em;';
seqBox.innerText = 'Click a match line to show its aligned sequence.';
plotDiv.parentNode.insertBefore(seqBox, plotDiv.nextSibling);
plotDiv.on('plotly_click', function(evt) {
    var pt = evt.points[0];
    var cd = pt.customdata;
    if (!cd) { return; }
    seqBox.innerText = cd[5];
});
""" % {'div_id': div_id}


def _self_gff_shapes(features, gff_color_dict, length_seq):
    """Diagonal square shading bands for self-dotplot GFF features."""
    shapes = []
    for feat_type, start, stop in features:
        feat_color, strength, zoom = gff_color_dict[feat_type.lower()]
        start_c = max(0, start - zoom - 0.5)
        stop_c = min(length_seq + 1, stop + zoom + 0.5)
        shapes.append(
            dict(
                type='rect',
                x0=start_c,
                x1=stop_c,
                y0=start_c,
                y1=stop_c,
                fillcolor=feat_color,
                opacity=strength,
                line=dict(width=0),
                layer='below',
            )
        )
    return shapes


def _pair_gff_shapes(features_one, features_two, gff_color_dict, len_one, len_two):
    """Full-height/full-width shading bands for paired dotplot GFF features.

    A feature on the x-axis sequence (name_one) spans the full plot height;
    a feature on the y-axis sequence (name_two) spans the full plot width.
    """
    shapes = []
    for feat_type, start, stop in features_one:
        feat_color, strength, zoom = gff_color_dict[feat_type.lower()]
        start_c = max(0, start - zoom - 0.5)
        stop_c = stop + zoom + 0.5
        shapes.append(
            dict(
                type='rect',
                x0=start_c,
                x1=stop_c,
                y0=0,
                y1=len_two + 1,
                fillcolor=feat_color,
                opacity=strength,
                line=dict(width=0),
                layer='below',
            )
        )
    for feat_type, start, stop in features_two:
        feat_color, strength, zoom = gff_color_dict[feat_type.lower()]
        start_c = max(0, start - zoom - 0.5)
        stop_c = stop + zoom + 0.5
        shapes.append(
            dict(
                type='rect',
                x0=0,
                x1=len_one + 1,
                y0=start_c,
                y1=stop_c,
                fillcolor=feat_color,
                opacity=strength,
                line=dict(width=0),
                layer='below',
            )
        )
    return shapes


def save_selfdotplot_html(
    name_seq,
    length_seq,
    x_lists,
    y_lists,
    x_lists_rc,
    y_lists_rc,
    fig_name,
    aa_bp_unit,
    line_col_for,
    line_col_rev,
    line_width,
    mirror_y_axis,
    title_length,
    title_clip_pos,
    label_size,
    seq=None,
    gff_features=None,
    gff_color_dict=None,
):
    """Render an interactive self-dotplot to an HTML file."""
    traces = _line_traces(
        x_lists,
        y_lists,
        x_lists_rc,
        y_lists_rc,
        line_col_for,
        line_col_rev,
        line_width,
        aa_bp_unit=aa_bp_unit,
        seq_x=seq,
        seq_y=seq,
    )
    fig = go.Figure(data=traces)

    if gff_features:
        fig.update_layout(
            shapes=_self_gff_shapes(gff_features, gff_color_dict, length_seq)
        )

    title = unicode_name(
        shorten_name(name_seq, max_len=title_length, title_clip_pos=title_clip_pos)
    )
    fig.update_layout(
        title=title,
        xaxis_title='[%s]' % aa_bp_unit,
        yaxis_title='[%s]' % aa_bp_unit,
        template='plotly_white',
        font=dict(size=label_size),
    )
    fig.update_xaxes(range=[0, length_seq + 1], constrain='domain', **_AXIS_BOX_STYLE)
    if mirror_y_axis:
        fig.update_yaxes(
            range=[0, length_seq + 1],
            constrain='domain',
            scaleanchor='x',
            scaleratio=1,
            **_AXIS_BOX_STYLE,
        )
    else:
        fig.update_yaxes(
            range=[length_seq + 1, 0],
            constrain='domain',
            scaleanchor='x',
            scaleratio=1,
            **_AXIS_BOX_STYLE,
        )

    fig.write_html(
        fig_name,
        include_plotlyjs=True,
        div_id=_DIV_ID,
        post_script=_click_seq_postscript(),
    )
    logging.info('Interactive HTML selfdotplot written to %s' % fig_name)
    return fig_name


def save_pairdotplot_html(
    name_one,
    name_two,
    len_one,
    len_two,
    x1,
    y1,
    x2,
    y2,
    fig_name,
    aa_bp_unit,
    line_col_for,
    line_col_rev,
    line_width,
    mirror_y_axis,
    title_length,
    title_clip_pos,
    label_size,
    x_label_pos_top,
    seq_one=None,
    seq_two=None,
    gff_features_one=None,
    gff_features_two=None,
    gff_color_dict=None,
):
    """Render an interactive paired dotplot to an HTML file.

    Each axis always spans its own sequence's length (matching matplotlib's
    non-collage pairdotplot behavior) - matches never get stretched past the
    end of the shorter sequence.
    """
    traces = _line_traces(
        x1,
        y1,
        x2,
        y2,
        line_col_for,
        line_col_rev,
        line_width,
        aa_bp_unit=aa_bp_unit,
        seq_x=seq_one,
        seq_y=seq_two,
    )
    fig = go.Figure(data=traces)

    if gff_features_one or gff_features_two:
        fig.update_layout(
            shapes=_pair_gff_shapes(
                gff_features_one or [], gff_features_two or [], gff_color_dict, len_one, len_two
            )
        )

    x_title = (
        unicode_name(shorten_name(name_one, max_len=title_length, title_clip_pos=title_clip_pos))
        + ' [%s]' % aa_bp_unit
    )
    y_title = (
        unicode_name(shorten_name(name_two, max_len=title_length, title_clip_pos=title_clip_pos))
        + ' [%s]' % aa_bp_unit
    )
    fig.update_layout(
        title='%s vs. %s' % (name_one, name_two),
        xaxis_title=x_title,
        yaxis_title=y_title,
        template='plotly_white',
        font=dict(size=label_size),
    )

    x_range = [0, len_one + 1]
    y_range = [0, len_two + 1] if mirror_y_axis else [len_two + 1, 0]

    fig.update_xaxes(
        range=x_range,
        constrain='domain',
        side='top' if x_label_pos_top else 'bottom',
        **_AXIS_BOX_STYLE,
    )
    fig.update_yaxes(
        range=y_range, constrain='domain', scaleanchor='x', scaleratio=1, **_AXIS_BOX_STYLE
    )

    fig.write_html(
        fig_name,
        include_plotlyjs=True,
        div_id=_DIV_ID,
        post_script=_click_seq_postscript(),
    )
    logging.info('Interactive HTML pairdotplot written to %s' % fig_name)
    return fig_name


def save_polydotplot_html(
    sequences,
    seq_dict,
    data_dict,
    fig_name,
    aa_bp_unit,
    line_col_for,
    line_col_rev,
    line_width,
    mirror_y_axis,
    title_length,
    title_clip_pos,
    label_size,
    plot_size,
):
    """Render an interactive all-against-all polydotplot grid to a single HTML file."""
    n = len(sequences)
    names = [
        unicode_name(
            shorten_name(seq_dict[s].id, max_len=title_length, title_clip_pos=title_clip_pos)
        )
        for s in sequences
    ]
    lengths = [len(seq_dict[s].seq) for s in sequences]
    seqs = [seq_dict[s].seq for s in sequences]

    fig = make_subplots(
        rows=n,
        cols=n,
        horizontal_spacing=min(0.02, 0.3 / n),
        vertical_spacing=min(0.02, 0.3 / n),
    )

    for idx in range(n):
        for jdx in range(idx, n):
            x_lists, y_lists, x_lists_rc, y_lists_rc = data_dict[(idx, jdx)]

            row, col = idx + 1, jdx + 1
            for trace in _line_traces(
                x_lists,
                y_lists,
                x_lists_rc,
                y_lists_rc,
                line_col_for,
                line_col_rev,
                line_width,
                aa_bp_unit=aa_bp_unit,
                seq_x=seqs[jdx],
                seq_y=seqs[idx],
            ):
                fig.add_trace(trace, row=row, col=col)
            fig.update_xaxes(range=[0, lengths[jdx] + 1], row=row, col=col, **_AXIS_BOX_STYLE)
            if mirror_y_axis:
                fig.update_yaxes(range=[0, lengths[idx] + 1], row=row, col=col, **_AXIS_BOX_STYLE)
            else:
                fig.update_yaxes(range=[lengths[idx] + 1, 0], row=row, col=col, **_AXIS_BOX_STYLE)

            if idx != jdx:
                # mirrored cell (swap x/y so the grid is symmetric)
                row_m, col_m = jdx + 1, idx + 1
                for trace in _line_traces(
                    y_lists,
                    x_lists,
                    y_lists_rc,
                    x_lists_rc,
                    line_col_for,
                    line_col_rev,
                    line_width,
                    aa_bp_unit=aa_bp_unit,
                    seq_x=seqs[idx],
                    seq_y=seqs[jdx],
                ):
                    fig.add_trace(trace, row=row_m, col=col_m)
                fig.update_xaxes(
                    range=[0, lengths[idx] + 1], row=row_m, col=col_m, **_AXIS_BOX_STYLE
                )
                if mirror_y_axis:
                    fig.update_yaxes(
                        range=[0, lengths[jdx] + 1], row=row_m, col=col_m, **_AXIS_BOX_STYLE
                    )
                else:
                    fig.update_yaxes(
                        range=[lengths[jdx] + 1, 0], row=row_m, col=col_m, **_AXIS_BOX_STYLE
                    )

    # sequence name labels along the bottom row and left column
    for jdx in range(n):
        fig.update_xaxes(title_text=names[jdx], row=n, col=jdx + 1)
    for idx in range(n):
        fig.update_yaxes(title_text=names[idx], row=idx + 1, col=1)

    fig.update_layout(
        template='plotly_white',
        font=dict(size=label_size),
        showlegend=False,
        height=max(600, plot_size * 100),
        width=max(600, plot_size * 100),
        title='Polydotplot',
    )

    fig.write_html(
        fig_name,
        include_plotlyjs=True,
        div_id=_DIV_ID,
        post_script=_click_seq_postscript(),
    )
    logging.info('Interactive HTML polydotplot written to %s' % fig_name)
    return fig_name
