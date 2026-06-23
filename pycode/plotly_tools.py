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

import numpy as np
import plotly.graph_objects as go

from ecqf_tools import abc_to_tau
from palette import row_colors

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
