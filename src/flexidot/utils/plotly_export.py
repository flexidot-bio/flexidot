###################################
#   Interactive HTML Dot Plots   #
###################################
"""
Interactive HTML dotplot rendering using Plotly.

Covers the core dotplot view (match lines with zoom/pan/hover/click) plus
GFF annotation shading for the self and paired plotting modes. LCS shading,
custom matrix shading, and collage layouts are not represented in HTML
output; GFF shading is also not yet supported for the poly mode.

Embedder contract: clicking a match line shows its position/sequence in an info
box directly below the plot (Go-to buttons above the sequence text). Clicking
one of those buttons (not the match line itself) sets ``window.__flexidotLastClick``
in the page:
    {x_seq, x_start, x_end, y_seq, y_start, y_end, length, unit, ts}
(``y_seq``/``y_start``/``y_end`` are null for a single-axis "go to" request.)
External tools (e.g. TEviewer) can poll this global (mirroring how they already
poll other Plotly-embedded state) to match ``x_seq``/``y_seq`` against their own
open records and jump to that region - only on the explicit button click, not on
every match click. The Go-to buttons (and this global) only appear when Flexidot
is run with ``--teviewer-integration``; standalone/browser use just shows the
click position/sequence with no host-app affordances.
"""

import json
import logging

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from flexidot.utils.utils import shorten_name, unicode_name

# solid box border around the plot area, similar to matplotlib's default axes spines
_AXIS_BOX_STYLE = {
    'showline': True,
    'linewidth': 1,
    'linecolor': 'black',
    'mirror': True,
}

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
    name_x='',
    name_y='',
):
    """Build Plotly line traces for forward and reverse-complement matches.

    All segments of one color are combined into a single trace (separated by
    None) so plots with many matches stay responsive. Each segment carries
    customdata (x/y sequence name, x/y start/end, length, aligned sequence
    text, orientation) used for the hover tooltip and the click-to-reveal
    sequence box / ``window.__flexidotLastClick`` (see module docstring).
    """
    traces = []
    for x_lines, y_lines, col, orientation in [
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
            point = [
                x0,
                x1,
                y0,
                y1,
                length,
                name_x,
                name_y,
                _extract_seq(seq_x, x0, x1) or 'n/a',
                _extract_seq(seq_y, y0, y1) or 'n/a',
                orientation,
            ]
            xs += [x0, x1, None]
            ys += [y0, y1, None]
            customdata += [point, point, point]
        traces.append(
            go.Scatter(
                x=xs,
                y=ys,
                mode='lines',
                line={'color': col, 'width': line_width},
                name=orientation,
                customdata=customdata,
                hovertemplate=(
                    '<b>%s match</b><br>' % orientation.title()
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


def _click_seq_postscript(unit='bp', div_id=_DIV_ID, show_goto_buttons=False):
    """JS injected after the plot: shows the clicked match's position/sequence in an
    info box directly below the plot. When ``show_goto_buttons`` is True (only set
    when Flexidot is launched with ``--teviewer-integration``), a "Go to <seq>
    <start>-<end>" button per axis is also shown (above the sequence text) that
    publishes ``window.__flexidotLastClick`` for the host app (e.g. TEviewer) to
    poll. Standalone/browser use never shows these buttons or touches that global."""
    if show_goto_buttons:
        make_button_fn = """
    function makeGoToButton(label, axis, seq, start, end) {
        var btn = document.createElement('button');
        btn.textContent = label;
        btn.style.cssText = 'margin: 0 8px 8px 0; font-family: monospace; font-size: 12px; padding: 2px 8px; cursor: pointer;';
        btn.addEventListener('click', function() {
            var payload = {
                x_seq: null, x_start: null, x_end: null,
                y_seq: null, y_start: null, y_end: null,
                length: Math.abs(end - start) + 1, unit: %s,
                ts: Date.now()
            };
            if (axis === 'x') {
                payload.x_seq = seq;
                payload.x_start = start;
                payload.x_end = end;
            } else {
                payload.y_seq = seq;
                payload.y_start = start;
                payload.y_end = end;
            }
            window.__flexidotLastClick = payload;
            btn.textContent = 'Sent \\u2192 ' + label;
            setTimeout(function() { btn.textContent = label; }, 1200);
        });
        return btn;
    }
    """ % json.dumps(unit)
        buttons_js = """
        var buttons = document.createElement('div');
        buttons.appendChild(
            makeGoToButton('Go to ' + cd[5] + ' ' + cd[0] + '-' + cd[1], 'x', cd[5], cd[0], cd[1])
        );
        buttons.appendChild(
            makeGoToButton('Go to ' + cd[6] + ' ' + cd[2] + '-' + cd[3], 'y', cd[6], cd[2], cd[3])
        );
        infoBox.appendChild(buttons);
        """
    else:
        make_button_fn = ''
        buttons_js = ''

    template = (
        """
(function() {
    var plotDiv = document.getElementById('%(div_id)s');
    if (!plotDiv) { return; }
    var infoBox = document.createElement('div');
    infoBox.id = '%(div_id)s-infobox';
    infoBox.style.cssText = 'font-family: monospace; font-size: 13px; padding: 10px 14px; margin-top: 6px; border: 1px solid #ccc; border-radius: 4px; background: #f7f7f7;';
    infoBox.innerText = 'Click a match line to show its position and sequence.';
    plotDiv.parentNode.insertBefore(infoBox, plotDiv.nextSibling);
"""
        + make_button_fn
        + """
    plotDiv.on('plotly_click', function(evt) {
        var pt = evt.points[0];
        var cd = pt.customdata;
        if (!cd) { return; }
        infoBox.innerHTML = '';
"""
        + buttons_js
        + """
        var text = document.createElement('div');
        text.style.cssText = 'white-space: pre-wrap; word-break: break-all;';
        text.textContent = 'X ' + cd[5] + ' ' + cd[0] + '-' + cd[1] + ': ' + cd[7] +
            '\\nY ' + cd[6] + ' ' + cd[2] + '-' + cd[3] + ': ' + cd[8];
        infoBox.appendChild(text);
    });
})();
"""
    )
    return template % {'div_id': div_id}


def _self_gff_shapes(features, gff_color_dict, length_seq):
    """Diagonal square shading bands for self-dotplot GFF features."""
    shapes = []
    for feat_type, start, stop in features:
        feat_color, strength, zoom = gff_color_dict[feat_type.lower()]
        start_c = max(0, start - zoom - 0.5)
        stop_c = min(length_seq + 1, stop + zoom + 0.5)
        shapes.append(
            {
                'type': 'rect',
                'x0': start_c,
                'x1': stop_c,
                'y0': start_c,
                'y1': stop_c,
                'fillcolor': feat_color,
                'opacity': strength,
                'line': {'width': 0},
                'layer': 'below',
            }
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
            {
                'type': 'rect',
                'x0': start_c,
                'x1': stop_c,
                'y0': 0,
                'y1': len_two + 1,
                'fillcolor': feat_color,
                'opacity': strength,
                'line': {'width': 0},
                'layer': 'below',
            }
        )
    for feat_type, start, stop in features_two:
        feat_color, strength, zoom = gff_color_dict[feat_type.lower()]
        start_c = max(0, start - zoom - 0.5)
        stop_c = stop + zoom + 0.5
        shapes.append(
            {
                'type': 'rect',
                'x0': 0,
                'x1': len_one + 1,
                'y0': start_c,
                'y1': stop_c,
                'fillcolor': feat_color,
                'opacity': strength,
                'line': {'width': 0},
                'layer': 'below',
            }
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
    teviewer_integration=False,
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
        name_x=name_seq,
        name_y=name_seq,
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
        font={'size': label_size},
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
        post_script=_click_seq_postscript(
            unit=aa_bp_unit, show_goto_buttons=teviewer_integration
        ),
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
    teviewer_integration=False,
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
        name_x=name_one,
        name_y=name_two,
    )
    fig = go.Figure(data=traces)

    if gff_features_one or gff_features_two:
        fig.update_layout(
            shapes=_pair_gff_shapes(
                gff_features_one or [],
                gff_features_two or [],
                gff_color_dict,
                len_one,
                len_two,
            )
        )

    x_title = (
        unicode_name(
            shorten_name(name_one, max_len=title_length, title_clip_pos=title_clip_pos)
        )
        + ' [%s]' % aa_bp_unit
    )
    y_title = (
        unicode_name(
            shorten_name(name_two, max_len=title_length, title_clip_pos=title_clip_pos)
        )
        + ' [%s]' % aa_bp_unit
    )
    fig.update_layout(
        title='%s vs. %s' % (name_one, name_two),
        xaxis_title=x_title,
        yaxis_title=y_title,
        template='plotly_white',
        font={'size': label_size},
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
        range=y_range,
        constrain='domain',
        scaleanchor='x',
        scaleratio=1,
        **_AXIS_BOX_STYLE,
    )

    fig.write_html(
        fig_name,
        include_plotlyjs=True,
        div_id=_DIV_ID,
        post_script=_click_seq_postscript(
            unit=aa_bp_unit, show_goto_buttons=teviewer_integration
        ),
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
    teviewer_integration=False,
):
    """Render an interactive all-against-all polydotplot grid to a single HTML file."""
    n = len(sequences)
    names = [
        unicode_name(
            shorten_name(
                seq_dict[s].id, max_len=title_length, title_clip_pos=title_clip_pos
            )
        )
        for s in sequences
    ]
    lengths = [len(seq_dict[s].seq) for s in sequences]
    seqs = [seq_dict[s].seq for s in sequences]
    record_ids = [seq_dict[s].id for s in sequences]

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
                name_x=record_ids[jdx],
                name_y=record_ids[idx],
            ):
                fig.add_trace(trace, row=row, col=col)
            fig.update_xaxes(
                range=[0, lengths[jdx] + 1], row=row, col=col, **_AXIS_BOX_STYLE
            )
            if mirror_y_axis:
                fig.update_yaxes(
                    range=[0, lengths[idx] + 1], row=row, col=col, **_AXIS_BOX_STYLE
                )
            else:
                fig.update_yaxes(
                    range=[lengths[idx] + 1, 0], row=row, col=col, **_AXIS_BOX_STYLE
                )

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
                    name_x=record_ids[idx],
                    name_y=record_ids[jdx],
                ):
                    fig.add_trace(trace, row=row_m, col=col_m)
                fig.update_xaxes(
                    range=[0, lengths[idx] + 1], row=row_m, col=col_m, **_AXIS_BOX_STYLE
                )
                if mirror_y_axis:
                    fig.update_yaxes(
                        range=[0, lengths[jdx] + 1],
                        row=row_m,
                        col=col_m,
                        **_AXIS_BOX_STYLE,
                    )
                else:
                    fig.update_yaxes(
                        range=[lengths[jdx] + 1, 0],
                        row=row_m,
                        col=col_m,
                        **_AXIS_BOX_STYLE,
                    )

    # sequence name labels along the bottom row and left column
    for jdx in range(n):
        fig.update_xaxes(title_text=names[jdx], row=n, col=jdx + 1)
    for idx in range(n):
        fig.update_yaxes(title_text=names[idx], row=idx + 1, col=1)

    fig.update_layout(
        template='plotly_white',
        font={'size': label_size},
        showlegend=False,
        height=max(600, plot_size * 100),
        width=max(600, plot_size * 100),
        title='Polydotplot',
    )

    fig.write_html(
        fig_name,
        include_plotlyjs=True,
        div_id=_DIV_ID,
        post_script=_click_seq_postscript(
            unit=aa_bp_unit, show_goto_buttons=teviewer_integration
        ),
    )
    logging.info('Interactive HTML polydotplot written to %s' % fig_name)
    return fig_name
