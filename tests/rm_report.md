# FlexiDot HTML Output

`-f html` swaps the rendering backend from matplotlib (static PNG/PDF/SVG) to
**Plotly**, producing a self-contained, interactive `.html` file per plot
instead of a static image. It's implemented in
`src/flexidot/utils/plotly_export.py` (566 lines) and hooked into the three
plotting modes in `plotting.py`.

## Where it plugs in

The k-mer matching stage (`find_match_pos_diag` / chunked matcher) runs
exactly as normal and produces the same `x1, y1, x2, y2` match-coordinate
arrays regardless of output format. Each of the three plot functions in
`plotting.py` (`selfdotplot`, `pairdotplot`, `polydotplot`) has an early
`if filetype == 'html': ... continue` branch (plotting.py:232, 815, 1532)
that hands those coordinates to Plotly instead of falling through to the
matplotlib code below it.

## Core building block: `_line_traces`

Every mode funnels through `_line_traces` (plotly_export.py:48). For each
match segment it takes only the **first and last point**
(`x_lines[ldx][0]`, `x_lines[ldx][-1]`) — segments are contiguous diagonal
runs, so a straight line between the two endpoints is visually identical to
plotting every point. Forward and reverse-complement matches become two
separate `go.Scatter` traces (one color each), with all of that color's
segments concatenated into a single trace separated by `None` gaps — this
keeps a plot with thousands of matches from creating thousands of separate
Plotly traces, which is what actually keeps large plots responsive in the
browser.

Each segment also carries `customdata`: both sequence names, start/end on
each axis, match length, the actual aligned subsequence text (via
`_extract_seq`), and orientation — this is what powers the hover tooltip and
the click info box.

## Three entry points, one shared structure

- **`save_selfdotplot_html`** — one HTML file per sequence.
- **`save_pairdotplot_html`** — one HTML file per sequence pair; each axis is
  scaled to its own sequence's length (never stretched to match the other).
- **`save_polydotplot_html`** — a single HTML file containing the full
  all-vs-all grid, built with `plotly.subplots.make_subplots` (one subplot
  per pair, mirrored across the diagonal by swapping x/y for the lower
  triangle).

All three call `fig.write_html(..., include_plotlyjs=True, post_script=...)`
— `include_plotlyjs=True` bundles the entire Plotly.js library inline (why
each file is ~3.5MB and fully self-contained/offline-viewable), and
`post_script` injects custom JS.

## Interactivity (`_click_seq_postscript`)

Beyond Plotly's built-in zoom/pan/hover, a small injected JS snippet adds a
click handler: clicking a match line shows an info box below the plot with
the exact coordinates and sequence text on both axes. With
`--teviewer-integration`, it additionally renders "Go to X / Go to Y"
buttons that publish `window.__flexidotLastClick` (a documented global) so
an external tool like TEviewer can poll it and jump to that region — this
only appears when that flag is passed; standalone use never touches that
global.

## GFF shading

Supported for self (`_self_gff_shapes`, diagonal square bands) and paired
(`_pair_gff_shapes`, full-height/full-width bands) plots via Plotly
`shapes`, but **not yet for polydotplot**.

## Known limitations (explicit, logged to the user)

- **No collage layout** — HTML always emits one file per sequence/pair, even
  if `-c` was passed (logged once: "Collage layout is not supported for
  HTML output").
- **No LCS shading or custom-matrix shading** — matplotlib-only features.
- These are documented directly in the module docstring and in the `-f` CLI
  help text, not silent gaps.
