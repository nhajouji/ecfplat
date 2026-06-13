"""
graphic_tools.py
================
All matplotlib plotting functions for ecfplat.

Convention
----------
- Every function takes plain data (numpy arrays, dicts, lists) and returns (fig, ax).
- All drawing is done through `ax.*` — never `plt.scatter()`, `plt.xlim()`, etc.
- Callers are responsible for calling `plt.close(fig)` after displaying.

Colour palette
--------------
BLUE   = "steelblue"   — curve, parallelogram, lattice points, horizontal isogeny edges
PURPLE = "purple"      — affine Fp points, classical picture
OLIVE  = "olive"       — distinguished points (origin, point at infinity)
ORANGE = "orange"      — vertical isogeny edges, chord-tangent sum
RED    = "red"         — highlights
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
from matplotlib.patches import Polygon, FancyArrow
from collections import defaultdict

from ecqf_tools import abc_to_tau, frob_to_mw_gens, mw_arr_from_gens

# ── Colour / style constants ──────────────────────────────────────────────────
_BLUE   = "steelblue"
_PURPLE = "purple"
_OLIVE  = "olive"
_ORANGE = "orange"
_RED    = "red"

_PARA_FACE = [0.85, 0.85, 0.95, 0.4]   # parallelogram fill
_PARA_EDGE = _BLUE
_PAD       = 0.1                         # axis padding around plots


# ─────────────────────────────────────────────────────────────────────────────
# Data helpers (no matplotlib)
# ─────────────────────────────────────────────────────────────────────────────

def fp_to_sym_fd(x: int, p: int) -> int:
    """Return x mod p in the symmetric range (-(p-1)//2, ..., p//2)."""
    xr = x % p
    return xr if 2 * xr < p else xr - p


def ecfp_aff_pts(fg: tuple, p: int) -> list:
    """Return all affine points on y^2 = x^3 + fx + g over F_p,
    with coordinates in the symmetric fundamental domain."""
    sqrt_dic = {y2: [] for y2 in range(p)}
    for y in range(p):
        sqrt_dic[pow(y, 2, p)].append(y)
    f, g = fg
    pts = []
    for x in range(p):
        xs = fp_to_sym_fd(x, p)
        y2 = (pow(xs, 3) + f * xs + g) % p
        pts += [(xs, fp_to_sym_fd(y, p)) for y in sqrt_dic[y2]]
    return pts


def qf_parallelogram_verts(qf: tuple) -> list:
    """Return the four vertices [0, 1, 1+τ, τ] of the fundamental parallelogram."""
    tau = abc_to_tau(qf)
    one = np.array([1.0, 0.0])
    return [np.zeros(2), one, one + tau, tau]


# ─────────────────────────────────────────────────────────────────────────────
# Lattice / parallelogram plots
# ─────────────────────────────────────────────────────────────────────────────

def qf_parallelogram_plot(qf: tuple):
    """Draw the fundamental parallelogram for the quadratic form qf.

    Returns (fig, ax).
    """
    verts = qf_parallelogram_verts(qf)
    xs = [v[0] for v in verts]
    ys = [v[1] for v in verts]
    xl, xr, ym = min(xs), max(xs), max(ys)

    fig, ax = plt.subplots(figsize=(5, 5 * ym))
    ax.add_patch(Polygon(verts, facecolor=_PARA_FACE, edgecolor=_PARA_EDGE, lw=1.5))
    ax.set_xlim(xl - _PAD, xr + _PAD)
    ax.set_ylim(-_PAD, ym + _PAD)
    ax.set_aspect("equal")
    ax.set_frame_on(False)
    return fig, ax


def ecqf_mw_lattice_plot(ecdata: dict, k: int):
    """Draw F_{p^k}-rational points inside the fundamental parallelogram.

    Points are computed from the MW generators derived from the Frobenius matrix.
    Returns (fig, ax, pts_array).
    """
    qf      = ecdata["qf"]
    frm     = ecdata["FrobMat"]
    mw_gens = frob_to_mw_gens(frm, k)
    pts     = mw_arr_from_gens(qf, mw_gens)
    pts_arr = np.array(pts) if len(pts) > 0 else np.empty((0, 2))

    fig, ax = qf_parallelogram_plot(qf)
    if len(pts_arr) > 0:
        ax.scatter(pts_arr[:, 0], pts_arr[:, 1], s=18, color=_BLUE, zorder=3)
    ax.scatter([0], [0], color=_OLIVE, s=30, zorder=4)
    tau = abc_to_tau(qf)
    ax.set_title(f"τ = {ecdata['tau_str']},  k = {k},  {len(pts_arr)} point(s)",
                 fontsize=10)
    return fig, ax, pts_arr


def ecqf_mwgen_arw_plot(ecdata: dict, k: int = 1):
    """Draw arrows representing the MW generators inside the fundamental parallelogram.

    Each generator is shown as an arrow from the origin pointing in the direction
    of the corresponding lattice vector. When there are two generators of distinct
    orders the arrow lengths are scaled to visually distinguish them (shorter arrow
    for the larger-order generator).

    Returns (fig, ax).
    """
    qf      = ecdata["qf"]
    frm     = ecdata["FrobMat"]
    mw_gens = frob_to_mw_gens(frm, k)

    verts = qf_parallelogram_verts(qf)
    one_arr, tau_arr = verts[1], verts[3]

    fig, ax = qf_parallelogram_plot(qf)

    # Compute unit-direction arrow for each generator
    gen_arrows = {}
    for vg in mw_gens:
        x, y     = vg
        vec      = x * one_arr + y * tau_arr
        norm     = np.linalg.norm(vec)
        gen_arrows[vg] = (vec / norm) * 0.5  # half-unit length

    # If two distinct orders, scale to visually differentiate
    orders = list(mw_gens.values())
    if len(set(orders)) == 2:
        sorted_gens = sorted(mw_gens, key=lambda v: mw_gens[v])
        gen_arrows[sorted_gens[0]] *= 1.2   # smaller order → longer arrow
        gen_arrows[sorted_gens[1]] *= 0.7   # larger order  → shorter arrow

    arrow_w = 0.04
    for vg, arr in gen_arrows.items():
        dx, dy = arr
        ax.add_patch(
            FancyArrow(
                0, 0, dx, dy,
                width=arrow_w,
                length_includes_head=True,
                head_width=arrow_w * 3,
                head_length=arrow_w * 3,
                color=_BLUE,
            )
        )

    ax.set_title(f"τ = {ecdata['tau_str']},  k = {k}  (generator directions)",
                 fontsize=10)
    return fig, ax


# ─────────────────────────────────────────────────────────────────────────────
# Classical F_p plot
# ─────────────────────────────────────────────────────────────────────────────

def ecfp_classic_plot(ecdata: dict):
    """Scatter plot of affine F_p points of an elliptic curve.

    Coordinates are placed in the symmetric fundamental domain [-(p-1)/2, p/2].
    The point at infinity is marked in olive in the top-right corner.

    Returns (fig, ax).
    """
    fg  = ecdata["coefs"]
    p   = ecdata["ap"][1]
    pts = ecfp_aff_pts(fg, p)
    xs  = [pt[0] for pt in pts]
    ys  = [pt[1] for pt in pts]

    amb = [(x, y)
           for x in range(-(p // 2), (p // 2) + 1)
           for y in range(-(p // 2), (p // 2) + 1)]
    ax_s, ay_s = zip(*amb)

    fig, ax = plt.subplots()
    ax.scatter(ax_s, ay_s, color="gray", alpha=0.3, s=8)
    ax.scatter(xs, ys, color=_PURPLE, s=18, zorder=3)
    ax.scatter([(p + 1) // 2], [(p + 1) // 2], color=_OLIVE, s=30, zorder=4)
    ax.set_xticks(range(-(p // 2), (p + 1) // 2))
    ax.set_yticks(range(-(p // 2), (p + 1) // 2))
    ax.set_xlim(-(p + 3) // 2, (p + 3) // 2)
    ax.set_ylim(-(p + 3) // 2, (p + 3) // 2)
    ax.set_aspect("equal")
    ax.set_frame_on(False)
    ax.set_title(ecdata["ec_eq"])
    return fig, ax


# ─────────────────────────────────────────────────────────────────────────────
# Upper half-plane plots
# ─────────────────────────────────────────────────────────────────────────────

def uhmodgam_fd_plot(height: float):
    """Draw the standard fundamental domain of SL(2,Z) up to the given height.

    Returns (fig, ax).
    """
    fig, ax = plt.subplots(figsize=(3, max(3, height)))
    ax.add_patch(mpl.patches.Rectangle(
        (-0.5, 0), width=1, height=height,
        facecolor=[0.7, 0.7, 0.7, 0.5], zorder=1
    ))
    ax.add_patch(mpl.patches.Circle(
        (0, 0), 1, facecolor="white", edgecolor="black", zorder=2
    ))
    ax.set_xlim(-1, 1)
    ax.set_ylim(0, height)
    ax.set_aspect("equal")
    return fig, ax


def qfs_uhfd_plot(qflist: list):
    """Plot a list of quadratic forms as points τ in the upper half-plane.

    Returns (fig, ax).
    """
    taus = np.array([abc_to_tau(qf) for qf in qflist])
    y_max = max(taus[:, 1])
    fig, ax = uhmodgam_fd_plot(y_max + _PAD)
    ax.scatter(taus[:, 0], taus[:, 1], s=30, color=_BLUE, zorder=3)
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(0, y_max + _PAD)
    return fig, ax


def qfs_uhfd_plot_w_highlights(qflist: list, qf0: tuple):
    """Like qfs_uhfd_plot but highlights qf0 in red.

    Returns (fig, ax).
    """
    tau0    = abc_to_tau(qf0)
    fig, ax = qfs_uhfd_plot(qflist)
    ax.scatter([tau0[0]], [tau0[1]], color=_RED, s=50, zorder=4)
    return fig, ax


# ─────────────────────────────────────────────────────────────────────────────
# Isogeny graph
# ─────────────────────────────────────────────────────────────────────────────

def _ring_layout(isoclass, l: int):
    """Compute (x, y) positions for every node in the degree-l isogeny graph.

    Layout: concentric rings by conductor (floor = innermost).
    Each ring is laid out independently over [0, 2π]; nodes with ascending
    parents are sorted by their parent's angle so children land near parents.

    Returns dict {qf: (x, y)} and dict {qf: angle}.
    """
    horz_data = isoclass.get_isog_neighbors_horz(l)
    asc_data  = isoclass.get_isog_neighbors_asc(l)
    qfs       = list(isoclass.qfs_ordered)
    qf_idx    = {qf: i for i, qf in enumerate(qfs)}
    cond_dict = isoclass.endo_cond_dict

    conds_sorted = sorted(set(cond_dict.values()))
    ring_radii   = {c: 1.8 + 2.2 * i for i, c in enumerate(conds_sorted)}

    node_angle = {}
    positions  = {}

    for cond in conds_sorted:
        nodes = [qf for qf in qfs if cond_dict[qf] == cond]
        m     = len(nodes)
        r     = ring_radii[cond]

        def sort_key(qf, _asc=asc_data, _ang=node_angle, _idx=qf_idx):
            if qf in _asc:
                return (0, _ang.get(_asc[qf], 0), _idx[qf])
            return (1, 0, _idx[qf])

        for i, qf in enumerate(sorted(nodes, key=sort_key)):
            angle          = 2 * np.pi * i / m
            node_angle[qf] = angle
            positions[qf]  = (r * np.cos(angle), r * np.sin(angle))

    return positions, node_angle


def isogeny_graph_figure(isoclass, l: int, qf_to_label: dict):
    """Draw the degree-l isogeny graph with a concentric-ring layout.

    - Lowest conductor (floor) at centre, higher conductors on outer rings.
    - Horizontal (cycle) edges in steelblue, vertical (tree) edges in orange.
    - Self-loops drawn as arcs near the node.
    - Node labels supplied by qf_to_label dict.

    Returns (fig, ax).
    """
    horz_data = isoclass.get_isog_neighbors_horz(l)
    asc_data  = isoclass.get_isog_neighbors_asc(l)
    qfs       = list(isoclass.qfs_ordered)
    n         = len(qfs)

    positions, _ = _ring_layout(isoclass, l)

    node_r = max(0.12, min(0.38, 4.5 / max(n, 12)))
    fs     = max(5, min(9, int(36 / max(n, 8))))

    fig, ax = plt.subplots(figsize=(8, 8))

    # Vertical edges
    for child, parent in asc_data.items():
        x1, y1 = positions[child]
        x2, y2 = positions[parent]
        ax.plot([x1, x2], [y1, y2], color=_ORANGE, lw=1.5, zorder=1)

    # Horizontal edges (deduplicate; self-loops as arcs)
    drawn = set()
    for qf, nbrs in horz_data.items():
        if qf in nbrs:
            x, y = positions[qf]
            offset = node_r * 2.2
            ax.add_patch(mpatches.Arc(
                (x + offset, y + offset),
                width=node_r * 2.8, height=node_r * 2.8,
                angle=45, theta1=0, theta2=330,
                color=_BLUE, lw=1.5, zorder=1,
            ))
        for nbr in nbrs:
            if nbr == qf:
                continue
            edge = frozenset([qf, nbr])
            if edge not in drawn:
                drawn.add(edge)
                x1, y1 = positions[qf]
                x2, y2 = positions[nbr]
                ax.plot([x1, x2], [y1, y2], color=_BLUE, lw=1.5, zorder=1)

    # Nodes and labels
    for qf in qfs:
        x, y = positions[qf]
        ax.add_patch(plt.Circle((x, y), node_r, color="white", ec="black",
                                lw=1.0, zorder=2))
        ax.text(x, y, qf_to_label[qf], ha="center", va="center",
                fontsize=fs, zorder=3)

    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"Degree-ℓ = {l} isogeny graph", fontsize=11)

    h_line = plt.Line2D([0], [0], color=_BLUE,   lw=2, label="Horizontal (cycle)")
    v_line = plt.Line2D([0], [0], color=_ORANGE, lw=2, label="Vertical (tree)")
    ax.legend(handles=[h_line, v_line], loc="upper right", fontsize=9)

    plt.tight_layout()
    return fig, ax
