"""
plotly_tools.py
===============
Interactive Plotly figures for ecfplat (Streamlit applets).

Convention
----------
- These functions return a `plotly.graph_objects.Figure` ready for
  `st.plotly_chart(fig, on_select="rerun", selection_mode="points")`.
- The matplotlib equivalents live in `graphic_tools.py`; this module is the
  interactive counterpart.
- Colour palette mirrors `graphic_tools.py` (steelblue points, red highlight).
"""

import math

import numpy as np
import plotly.graph_objects as go

from ecqf_tools import abc_to_tau
from palette import row_colors, text_color_for

_BLUE = "steelblue"
_RED  = "red"
_FD_FILL = "rgba(178,178,210,0.25)"   # faint fill for the FD interior
_FD_LINE = "black"

# Corners of the fundamental domain where the unit arc meets Re = ±1/2:
#   e^{iπ/3}   = ( 1/2, √3/2)
#   e^{i2π/3}  = (-1/2, √3/2)
_ARC_Y = np.sqrt(3) / 2.0


def _qf_tau(qf, invert: bool = False):
    """τ = (Re, Im) for a quadratic form, forced into the upper half-plane.

    `abc_to_tau` returns one root of A x² + B x + C; for a (reduced) form the
    two roots are complex conjugates, so we take |Im| to land in the UHP.

    When ``invert`` is True we apply S: τ ↦ −1/τ, which on the form (a, b, c)
    is just (a, b, c) ↦ (c, −b, a). This sends the cusp (the high-up trivial
    class) to the origin and packs every reduced form into the unit disc, so
    the picture has a fixed size regardless of the discriminant.
    """
    a, b, c = qf
    if invert:
        a, b, c = c, -b, a
    x, y = abc_to_tau((a, b, c))
    return float(x), abs(float(y))


def _arc_points(t0: float, t1: float, n: int = 80):
    """Polyline on the unit circle for angles in [t0, t1]."""
    thetas = np.linspace(t0, t1, n)
    return np.cos(thetas), np.sin(thetas)


def _circle_arc(cx: float, cy: float, r: float, t0: float, t1: float, n: int = 60):
    """Polyline on the circle of radius r about (cx, cy) for angles [t0, t1]."""
    thetas = np.linspace(t0, t1, n)
    return cx + r * np.cos(thetas), cy + r * np.sin(thetas)


def _fd_shapes(y_top: float):
    """Layout shapes for the FD boundary: faint interior fill + black border.

    Drawn as layout shapes (not traces) so the figure has a single trace
    (the τ points), keeping every selection's curve number == 0.
    """
    ax, ay = _arc_points(np.pi / 3.0, 2.0 * np.pi / 3.0)

    # Filled interior: down the left wall, across the arc, up the right wall.
    fill_path = f"M {-0.5},{y_top} L {-0.5},{_ARC_Y} "
    fill_path += " ".join(f"L {x:.5f},{y:.5f}" for x, y in zip(ax, ay))
    fill_path += f" L {0.5},{y_top} Z"

    # Arc-only path for the visible bottom border.
    arc_path = "M " + " L ".join(f"{x:.5f},{y:.5f}" for x, y in zip(ax, ay))

    return [
        dict(type="path", path=fill_path, fillcolor=_FD_FILL,
             line=dict(width=0), layer="below"),
        dict(type="line", x0=-0.5, y0=_ARC_Y, x1=-0.5, y1=y_top,
             line=dict(color=_FD_LINE, width=1)),
        dict(type="line", x0=0.5, y0=_ARC_Y, x1=0.5, y1=y_top,
             line=dict(color=_FD_LINE, width=1)),
        dict(type="path", path=arc_path, line=dict(color=_FD_LINE, width=1)),
    ]


def _path_from_xy(xs, ys, move=True):
    """SVG path string from coordinate arrays."""
    cmds = [("M" if (move and i == 0) else "L") + f" {x:.5f},{y:.5f}"
            for i, (x, y) in enumerate(zip(xs, ys))]
    return " ".join(cmds)


def _fd_inverted_shapes():
    """Boundary shapes for the S-image (τ ↦ −1/τ) of the fundamental domain.

    That image is the ideal curvilinear triangle with vertices e^{iπ/3},
    e^{i2π/3} and 0, bounded by three unit-radius arcs:
      - top:   |τ| = 1   between e^{iπ/3} and e^{i2π/3}
      - right: circle about (+1, 0)  (image of the wall Re τ = +1/2)
      - left:  circle about (−1, 0)  (image of the wall Re τ = −1/2)
    Every reduced form lands inside this region, so the picture is bounded.
    """
    tx, ty = _arc_points(2.0 * np.pi / 3.0, np.pi / 3.0)        # e^{i2π/3} → e^{iπ/3}
    rx, ry = _circle_arc(1.0, 0.0, 1.0, 2.0 * np.pi / 3.0, np.pi)   # e^{iπ/3} → 0
    lx, ly = _circle_arc(-1.0, 0.0, 1.0, 0.0, np.pi / 3.0)         # 0 → e^{i2π/3}

    # Filled interior: 0 → (left arc) → top arc → (right arc) → 0.
    fill_path = (_path_from_xy(lx, ly)
                 + " " + _path_from_xy(tx, ty, move=False)
                 + " " + _path_from_xy(rx, ry, move=False) + " Z")

    return [
        dict(type="path", path=fill_path, fillcolor=_FD_FILL,
             line=dict(width=0), layer="below"),
        dict(type="path", path=_path_from_xy(tx, ty),
             line=dict(color=_FD_LINE, width=1)),
        dict(type="path", path=_path_from_xy(rx, ry),
             line=dict(color=_FD_LINE, width=1)),
        dict(type="path", path=_path_from_xy(lx, ly),
             line=dict(color=_FD_LINE, width=1)),
    ]


def fd_points_figure(df, selected_row=None, view: str = "standard",
                     colors=None, y_cap: float = 4.0):
    """Interactive fundamental-domain picture, one point per dataframe row.

    Parameters
    ----------
    df : pandas.DataFrame
        Must have columns `qf_coefs`, `ec_invs`. Points are plotted in df-row
        order, so the Plotly point index equals the df row index.
    selected_row : int or None
        Row to highlight (enlarged, with a red ring; keeps its own colour).
    colors : list[str] or None
        Per-row marker colours (df-row order). Defaults to `row_colors(n)` so
        the picture matches the isogeny graph's per-class colours.
    view : "standard" | "inverted"
        "standard" draws the usual fundamental domain (Im τ capped at `y_cap`;
        the trivial class sits high near the cusp and may need autoscale to
        see). "inverted" applies S: τ ↦ −1/τ, packing every class into the
        ideal triangle inside the unit disc — bounded and a fixed size, so the
        trivial class is always visible near the origin.
    y_cap : float
        Standard view only: initial cap on the visible Im τ. Points above the
        cap still exist — use the Plotly modebar's autoscale or pan to reach
        them, or switch to the inverted view.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    invert = (view == "inverted")
    xs, ys, texts = [], [], []
    for i, row in enumerate(df.itertuples(index=False)):
        qf = row.qf_coefs
        x, y = _qf_tau(qf, invert=invert)
        xs.append(x)
        ys.append(y)
        # Hover always shows the original form so the curve is identifiable
        # regardless of which view moved the point.
        texts.append(
            f"row {i}<br>qf = {tuple(qf)}<br>ec_invs = {row.ec_invs}"
        )

    n = len(xs)
    marker_colors = list(colors) if colors is not None else row_colors(n)
    sizes       = [12] * n
    line_colors = ["black"] * n
    line_widths = [0.7] * n
    if selected_row is not None and 0 <= selected_row < n:
        # Keep the row's own colour; mark selection with size + a red ring.
        sizes[selected_row]       = 20
        line_colors[selected_row] = _RED
        line_widths[selected_row] = 3.0

    fig = go.Figure(
        go.Scatter(
            x=xs, y=ys,
            mode="markers",
            marker=dict(color=marker_colors, size=sizes,
                        line=dict(color=line_colors, width=line_widths)),
            text=texts,
            hovertemplate="%{text}<extra></extra>",
            customdata=list(range(n)),
            # Keep every point fully opaque even while one is "selected" in the
            # Plotly sense — our ring/size is the only highlight we want.
            selected=dict(marker=dict(opacity=1.0)),
            unselected=dict(marker=dict(opacity=1.0)),
        )
    )

    if invert:
        shapes = _fd_inverted_shapes()
        x_range, y_range = [-0.62, 0.62], [-0.02, 1.06]
    else:
        max_y = max(ys) if ys else 1.0
        # Draw the domain walls all the way past the tallest point so that
        # panning/zooming up always stays inside the fundamental domain.
        shape_top = max(max_y + 0.6, _ARC_Y + 0.4)
        shapes = _fd_shapes(shape_top)
        # The *initial* view is still capped so the cluster is legible; the
        # high cusp points (e.g. the trivial class) are reached by panning up
        # or autoscaling — or via the unit-disc view.
        view_top = max(min(max_y + 0.4, y_cap), _ARC_Y + 0.4)
        x_range, y_range = [-1.15, 1.15], [0, view_top]

    fig.update_layout(
        shapes=shapes,
        margin=dict(l=10, r=10, t=10, b=10),
        showlegend=False,
        plot_bgcolor="white",
        height=460,
    )
    fig.update_xaxes(
        range=x_range,
        zeroline=False, showgrid=False,
        title_text="Re τ",
    )
    fig.update_yaxes(
        range=y_range,
        scaleanchor="x", scaleratio=1,
        zeroline=False, showgrid=False,
        title_text="Im τ",
    )
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Isogeny graph — structure-aware layout (volcanoes / cycles)
# ─────────────────────────────────────────────────────────────────────────────

_H_EDGE = "steelblue"   # horizontal (rim / cycle) isogenies
_V_EDGE = "orange"      # vertical (descending) isogenies

_NODE_SPACING = 1.0     # target distance between adjacent nodes on a cycle/rim
_LEVEL_GAP    = 1.15    # radial gap between volcano levels
_COMP_GAP     = 1.4     # gap between packed components


def _ngon_radius(n: int) -> float:
    """Radius of a regular n-gon whose adjacent vertices are _NODE_SPACING apart."""
    if n <= 1:
        return 0.0
    if n == 2:
        return _NODE_SPACING / 2.0
    return _NODE_SPACING / (2.0 * math.sin(math.pi / n))


def _ngon_positions(cycle, center=(0.0, 0.0)):
    """Positions of a cycle as a regular n-gon (start at top, go clockwise)."""
    n = len(cycle)
    cx, cy = center
    if n == 1:
        return {cycle[0]: (cx, cy)}, 0.0
    R = _ngon_radius(n)
    pos = {}
    for i, qf in enumerate(cycle):
        ang = math.pi / 2 - 2 * math.pi * i / n
        pos[qf] = (cx + R * math.cos(ang), cy + R * math.sin(ang))
    return pos, R


def _vl(c: int, l: int) -> int:
    """l-adic valuation."""
    v = 0
    while c % l == 0:
        c //= l
        v += 1
    return v


def _chi_str(a, p):
    """Characteristic polynomial χ(x) = x² − a·x + p as a tidy unicode string."""
    s = "x²"
    b = -a  # coefficient of x in x² − a x + p
    if b == 1:
        s += " + x"
    elif b == -1:
        s += " − x"
    elif b > 0:
        s += f" + {b}x"
    elif b < 0:
        s += f" − {-b}x"
    return s + f" + {p}"


def _dedupe_cycle(cyc):
    """Order-preserving de-duplication.

    `isog_cycle` returns `[v, v]` for a vertex whose only ℓ-isogenies are
    self-isogenies (a self-loop). Collapse that to a single node so it is laid
    out as one vertex; the self-loop is drawn separately from the multiplicity
    data.
    """
    return list(dict.fromkeys(cyc))


def _horiz_edge_mult(isoclass, l):
    """Horizontal isogeny multiplicities from the (repeat-preserving) data.

    Returns (pair_mult, self_mult):
      pair_mult: {frozenset({v1, v2}): count}  for v1 ≠ v2
      self_mult: {v: count}                    self-isogenies (self-loops)
    The counts are exactly what the adjacency matrix records.
    """
    from collections import Counter
    horiz = isoclass.get_isog_neighbors_horz(l)
    pair_mult, self_mult = {}, {}
    for v in isoclass.qfs_ordered:
        for u, c in Counter(horiz[v]).items():
            if u == v:
                self_mult[v] = c
            else:
                key = frozenset((u, v))
                pair_mult[key] = max(pair_mult.get(key, 0), c)
    return pair_mult, self_mult


def _cycle_components(isoclass, l):
    """(endo_cond, ordered_cycle) for each component when l ∤ cond."""
    cond_dict = isoclass.endo_cond_dict
    seen, comps = set(), []
    for qf in isoclass.qfs_ordered:
        if qf in seen:
            continue
        cyc = _dedupe_cycle(isoclass.isog_cycle(qf, l))
        seen |= set(cyc)
        comps.append((cond_dict[qf], cyc))
    comps.sort(key=lambda c: (c[0], -len(c[1])))
    return comps


def _volcano_components(isoclass, l):
    """One dict per volcano (connected component) when l | cond.

    Each: {positions, h_edges, v_edges, radius} centred at the origin.
    Surface (maximal at l) forms the rim cycle; descending trees follow the
    ascending-isogeny data, one level per power of l in the conductor.
    """
    cond_dict = isoclass.endo_cond_dict
    horiz = isoclass.get_isog_neighbors_horz(l)
    asc   = isoclass.get_isog_neighbors_asc(l)        # child -> parent
    children = {}
    for child, parent in asc.items():
        children.setdefault(parent, []).append(child)

    surface = [qf for qf in isoclass.qfs_ordered if _vl(cond_dict[qf], l) == 0]
    seen, rims = set(), []
    for qf in surface:
        if qf in seen:
            continue
        cyc = _dedupe_cycle(isoclass.isog_cycle(qf, l))
        seen |= set(cyc)
        rims.append(cyc)

    comps = []
    for rim in rims:
        m = len(rim)
        pos, h_edges, v_edges = {}, [], []
        R0 = _ngon_radius(m)
        for i, qf in enumerate(rim):
            ang = math.pi / 2 - (2 * math.pi * i / m if m > 1 else 0.0)
            pos[qf] = (R0 * math.cos(ang), R0 * math.sin(ang))
        if m >= 3:
            h_edges = [(rim[i], rim[(i + 1) % m]) for i in range(m)]
        elif m == 2:
            h_edges = [(rim[0], rim[1])]

        max_level = [0]

        def _descend(parent, base, half, level):
            kids = children.get(parent, [])
            n = len(kids)
            for j, kid in enumerate(kids):
                lo = base - half + (2 * half) * (j / n)
                hi = base - half + (2 * half) * ((j + 1) / n)
                ang = 0.5 * (lo + hi)
                r = R0 + level * _LEVEL_GAP
                pos[kid] = (r * math.cos(ang), r * math.sin(ang))
                v_edges.append((parent, kid))
                max_level[0] = max(max_level[0], level)
                _descend(kid, ang, 0.5 * (hi - lo), level + 1)

        for i, qf in enumerate(rim):
            base = math.pi / 2 - (2 * math.pi * i / m if m > 1 else 0.0)
            half = math.pi / m if m > 1 else math.pi
            _descend(qf, base, half, 1)

        comps.append({
            "positions": pos, "h_edges": h_edges, "v_edges": v_edges,
            "radius": R0 + max_level[0] * _LEVEL_GAP,
        })
    return comps


def _pack_grid(comp_pos_radii, ncols=None):
    """Lay component centres on a grid; return merged positions."""
    n = len(comp_pos_radii)
    if ncols is None:
        ncols = max(1, math.ceil(math.sqrt(n)))
    cell = 2 * max((r for _, r in comp_pos_radii), default=1.0) + _COMP_GAP
    merged = {}
    for idx, (pos, _r) in enumerate(comp_pos_radii):
        cx, cy = (idx % ncols) * cell, -(idx // ncols) * cell
        for qf, (x, y) in pos.items():
            merged[qf] = (cx + x, cy + y)
    return merged


def _pack_grouped(groups, ncols=6):
    """Stack groups vertically, each as a wrapped grid; return (positions, annotations).

    groups: list of (label, [(positions, radius), ...]).
    annotations: list of (x, y, label) placed at the left of each band.
    """
    maxR = max((r for _, comps in groups for _, r in comps), default=1.0)
    cell = 2 * maxR + _COMP_GAP
    merged, annotations = {}, []
    y = 0.0
    for label, comps in groups:
        rows = max(1, math.ceil(len(comps) / ncols))
        annotations.append((-0.5 * cell, y + maxR, label))
        for idx, (pos, _r) in enumerate(comps):
            cx = (idx % ncols) * cell
            cy = y - (idx // ncols) * cell
            for qf, (x, yy) in pos.items():
                merged[qf] = (cx + x, cy + yy)
        y -= rows * cell + 0.6 * _COMP_GAP
    return merged, annotations


def _edge_shapes(pos, edges, color):
    return [dict(type="line", x0=pos[a][0], y0=pos[a][1],
                 x1=pos[b][0], y1=pos[b][1],
                 line=dict(color=color, width=1.5), layer="below")
            for a, b in edges if a in pos and b in pos]


def isogeny_graph_figure(isoclass, l, colors, qf_to_row, qf_to_hover,
                         selected_row=None):
    """Clickable Plotly l-isogeny graph laid out by its volcano / cycle structure.

    - l | cond(Z[π])  → disjoint volcanoes (rim cycle + descending trees).
    - l ∤ cond         → disjoint regular-n-gon cycles, grouped by endo ring.

    Nodes are labelled by their df-row number and carry `customdata` = that row,
    so a click selects the class via the same mechanism as the FD picture.

    `colors` is the per-class palette (df-row order); `qf_to_row` maps each form
    to its df-row index; `qf_to_hover` maps each form to its hover string.
    """
    is_volcano = (isoclass.cond % l == 0)
    v_edges, annotations = [], []

    if is_volcano:
        comps = _volcano_components(isoclass, l)
        pos = _pack_grid([(c["positions"], c["radius"]) for c in comps])
        for c in comps:
            v_edges += c["v_edges"]
    else:
        comps = _cycle_components(isoclass, l)
        groups_map = {}
        for ec, cyc in comps:
            pos_c, r = _ngon_positions(cyc)
            groups_map.setdefault(ec, []).append((pos_c, r, cyc))
        groups = []
        for ec in sorted(groups_map):
            groups.append((f"End. conductor {ec}",
                           [(p, r) for (p, r, _c) in groups_map[ec]]))
        pos, annotations = _pack_grouped(groups)

    # Horizontal edges and self-loops come from the multiplicity data, so a
    # double/triple edge is drawn once with an "×N" label rather than as
    # overlapping parallel lines.
    pair_mult, self_mult = _horiz_edge_mult(isoclass, l)
    h_edges = [tuple(pair) for pair in pair_mult]

    # ── Nodes ────────────────────────────────────────────────────────────────
    qfs = [qf for qf in isoclass.qfs_ordered if qf in pos]
    n = len(qfs)
    size = max(12, min(26, int(420 / max(n, 16))))
    fs   = max(7, min(12, int(220 / max(n, 16))))

    xs, ys, node_colors, txt_colors, labels, hover, cdata = [], [], [], [], [], [], []
    sizes, line_colors, line_widths = [], [], []
    for qf in qfs:
        x, y = pos[qf]
        xs.append(x); ys.append(y)
        row = qf_to_row[qf]
        col = colors[row]
        node_colors.append(col)
        txt_colors.append(text_color_for(col))
        labels.append(str(row))
        hover.append(qf_to_hover.get(qf, ""))
        cdata.append(row)
        sel = (selected_row == row)
        sizes.append(size + 8 if sel else size)
        line_colors.append(_RED if sel else "black")
        line_widths.append(3.0 if sel else 0.8)

    node_trace = go.Scatter(
        x=xs, y=ys, mode="markers+text",
        marker=dict(color=node_colors, size=sizes,
                    line=dict(color=line_colors, width=line_widths)),
        text=labels, textposition="middle center",
        textfont=dict(size=fs, color=txt_colors),
        customdata=cdata,
        hovertext=hover, hoverinfo="text",
        selected=dict(marker=dict(opacity=1.0)),
        unselected=dict(marker=dict(opacity=1.0)),
    )

    # ── Edges, self-loops, multiplicity labels ────────────────────────────────
    shapes = (_edge_shapes(pos, h_edges, _H_EDGE)
              + _edge_shapes(pos, v_edges, _V_EDGE))
    mult_labels = []   # (x, y, "×N")

    # Direction "outward" from the layout centroid, for placing self-loops.
    cx0 = sum(p[0] for p in pos.values()) / max(len(pos), 1)
    cy0 = sum(p[1] for p in pos.values()) / max(len(pos), 1)
    loop_r = 0.34
    loop_extents = []
    for v, m in self_mult.items():
        if v not in pos:
            continue
        x, y = pos[v]
        dx, dy = x - cx0, y - cy0
        d = math.hypot(dx, dy)
        ux, uy = (dx / d, dy / d) if d > 1e-9 else (0.0, 1.0)
        lcx, lcy = x + ux * (loop_r + 0.16), y + uy * (loop_r + 0.16)
        shapes.append(dict(type="circle", xref="x", yref="y",
                           x0=lcx - loop_r, y0=lcy - loop_r,
                           x1=lcx + loop_r, y1=lcy + loop_r,
                           line=dict(color=_H_EDGE, width=1.5), layer="below"))
        loop_extents += [(lcx - loop_r, lcy - loop_r), (lcx + loop_r, lcy + loop_r)]
        if m > 1:
            mult_labels.append((x + ux * (2 * loop_r + 0.45),
                                y + uy * (2 * loop_r + 0.45), f"×{m}"))

    for pair, m in pair_mult.items():
        if m > 1:
            a, b = tuple(pair)
            if a in pos and b in pos:
                mx = 0.5 * (pos[a][0] + pos[b][0])
                my = 0.5 * (pos[a][1] + pos[b][1])
                mult_labels.append((mx, my, f"×{m}"))

    group_anns = [dict(x=ax, y=ay, text=lbl, showarrow=False,
                       xanchor="left", font=dict(size=11, color="gray"))
                  for (ax, ay, lbl) in annotations]
    mult_anns = [dict(x=mx, y=my, text=t, showarrow=False,
                      font=dict(size=12, color=_H_EDGE),
                      bgcolor="rgba(255,255,255,0.75)")
                 for (mx, my, t) in mult_labels]

    a_tr, p_ch = isoclass.ap
    title = f"{l}-isogeny graph for χ(x) = {_chi_str(a_tr, p_ch)}"

    fig = go.Figure(node_trace)
    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor="center",
                   font=dict(size=16, color="black")),
        shapes=shapes,
        annotations=group_anns + mult_anns,
        margin=dict(l=10, r=10, t=44, b=10),
        showlegend=False, plot_bgcolor="white", height=560,
    )
    # Equal aspect, no axes (include self-loop extents in the range).
    xs_all = [p[0] for p in pos.values()] + [e[0] for e in loop_extents]
    ys_all = [p[1] for p in pos.values()] + [e[1] for e in loop_extents]
    pad = 1.0
    fig.update_xaxes(visible=False,
                     range=[min(xs_all) - pad, max(xs_all) + pad])
    fig.update_yaxes(visible=False, scaleanchor="x", scaleratio=1,
                     range=[min(ys_all) - pad, max(ys_all) + pad])
    return fig
