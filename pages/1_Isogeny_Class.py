import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "pycode"))

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict

from ecqf_tools import ECQFIsogenyClass, ap_in_pc_data
from nt import primeQ

# ── Isogeny graph layout & plot ───────────────────────────────────────────────

def isogeny_graph_figure(isoclass, l, qf_to_label):
    """
    Draw the degree-l isogeny graph with a concentric-ring layout:
    - lowest conductor (floor) nearest the centre
    - higher conductors on outer rings
    - horizontal (cycle) edges in steelblue, vertical (tree) edges in orange

    Layout: each conductor ring is laid out independently over 0→2π.
    Nodes with ascending parents are sorted by their parent's angle so that
    children land angularly near their parent; roots keep qfs_ordered order.
    """
    horz_data = isoclass.get_isog_neighbors_horz(l)
    asc_data  = isoclass.get_isog_neighbors_asc(l)

    qfs = list(isoclass.qfs_ordered)
    n   = len(qfs)
    qf_idx = {qf: i for i, qf in enumerate(qfs)}
    cond_dict = isoclass.endo_cond_dict

    # ── Ring radii: one per distinct conductor value, floor at centre ─────────
    conds_sorted = sorted(set(cond_dict.values()))
    ring_radii   = {c: 1.8 + 2.2 * i for i, c in enumerate(conds_sorted)}

    # ── Assign positions ring by ring, inside out ─────────────────────────────
    node_angle = {}   # qf -> assigned angle (radians)
    positions  = {}

    for cond in conds_sorted:
        nodes = [qf for qf in qfs if cond_dict[qf] == cond]
        m     = len(nodes)
        r     = ring_radii[cond]

        # Sort: nodes with ascending parents first (by parent angle),
        # then roots in their original qfs_ordered position.
        def sort_key(qf):
            if qf in asc_data:
                return (0, node_angle.get(asc_data[qf], 0), qf_idx[qf])
            return (1, 0, qf_idx[qf])

        nodes_sorted = sorted(nodes, key=sort_key)

        for i, qf in enumerate(nodes_sorted):
            angle          = 2 * np.pi * i / m
            node_angle[qf] = angle
            positions[qf]  = (r * np.cos(angle), r * np.sin(angle))

    # ── Figure ────────────────────────────────────────────────────────────────
    node_r  = max(0.12, min(0.38, 4.5 / max(n, 12)))
    fs      = max(5, min(9, int(36 / max(n, 8))))
    fig, ax = plt.subplots(figsize=(8, 8))

    # Vertical edges
    for child, parent in asc_data.items():
        x1, y1 = positions[child]
        x2, y2 = positions[parent]
        ax.plot([x1, x2], [y1, y2], color="orange", lw=1.5, zorder=1)

    # Horizontal edges (deduplicate; draw self-loops as arcs)
    drawn = set()
    for qf, nbrs in horz_data.items():
        n_self = nbrs.count(qf)
        if n_self:
            x, y = positions[qf]
            offset = node_r * 2.2
            loop = mpatches.Arc((x + offset, y + offset),
                                width=node_r * 2.8, height=node_r * 2.8,
                                angle=45, theta1=0, theta2=330,
                                color="steelblue", lw=1.5, zorder=1)
            ax.add_patch(loop)
        for nbr in nbrs:
            if nbr == qf:
                continue
            edge = frozenset([qf, nbr])
            if edge not in drawn:
                drawn.add(edge)
                x1, y1 = positions[qf]
                x2, y2 = positions[nbr]
                ax.plot([x1, x2], [y1, y2], color="steelblue", lw=1.5, zorder=1)

    # Nodes and labels
    for qf in qfs:
        x, y = positions[qf]
        circle = plt.Circle((x, y), node_r, color="white", ec="black",
                             lw=1.0, zorder=2)
        ax.add_patch(circle)
        ax.text(x, y, qf_to_label[qf], ha="center", va="center",
                fontsize=fs, zorder=3)

    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"Degree-ℓ = {l} isogeny graph", fontsize=11)

    h_line = plt.Line2D([0], [0], color="steelblue", lw=2, label="Horizontal (cycle)")
    v_line = plt.Line2D([0], [0], color="orange",    lw=2, label="Vertical (tree)")
    ax.legend(handles=[h_line, v_line], loc="upper right", fontsize=9)

    plt.tight_layout()
    return fig


# ── Session state ─────────────────────────────────────────────────────────────
if "isoclass" not in st.session_state:
    st.session_state.isoclass = None
if "selected_row" not in st.session_state:
    st.session_state.selected_row = None

# ── Pre-fill from EC Search navigation ───────────────────────────────────────
prefill = st.session_state.pop("ic_prefill", None)

# ── Sidebar: (a, p) input ─────────────────────────────────────────────────────
with st.sidebar:
    st.title("Isogeny Class")
    st.markdown("Enter a pair *(a, p)* with *p* prime and *a² < 4p*.")
    a_input = st.number_input("a", value=prefill["a"] if prefill else -4, step=1)
    p_input = st.number_input("p", value=prefill["p"] if prefill else 5, step=1, min_value=2)
    load = st.button("Load isogeny class", use_container_width=True)

    if prefill:
        load = True

    if load:
        a, p = int(a_input), int(p_input)
        if not ap_in_pc_data((a, p)):
            st.error(
                f"(a, p) = ({a}, {p}) is not in the precomputed data.\n\n"
                "Check that p is prime, a² < 4p, and the pair is within "
                "the precomputed range (4 ≤ p ≤ 1024)."
            )
            st.session_state.isoclass = None
            st.session_state.selected_row = None
        else:
            with st.spinner("Loading…"):
                st.session_state.isoclass = ECQFIsogenyClass(a, p)
                st.session_state.selected_row = None
            st.success(f"Loaded (a, p) = ({a}, {p})")

# ── Main area ─────────────────────────────────────────────────────────────────
isoclass = st.session_state.isoclass

if isoclass is None:
    st.info("Enter a pair (a, p) in the sidebar and click **Load isogeny class**.")
    st.stop()

a, p = isoclass.ap
st.header(f"Isogeny class  (a, p) = ({a}, {p})")
st.caption(f"Discriminant: {isoclass.disc}   |   # lattice classes: {len(isoclass.qfs_all)}")

# ── Tabs ──────────────────────────────────────────────────────────────────────
tab_df, tab_isograph = st.tabs(["Isogeny class table", "Isogeny graph"])

# ── Tab 1: DataFrame ──────────────────────────────────────────────────────────
with tab_df:
    st.markdown("Click a row to select it and view the corresponding elliptic curve.")
    df = isoclass.ecqf_df()

    event = st.dataframe(
        df,
        use_container_width=True,
        hide_index=True,
        on_select="rerun",
        selection_mode="single-row",
    )

    selected_rows = event.selection.rows
    if selected_rows:
        st.session_state.selected_row = selected_rows[0]

    if st.session_state.selected_row is not None:
        row = df.iloc[st.session_state.selected_row]
        f_sel, g_sel = tuple(row["EC_coefs"])
        st.success(
            f"Selected row {st.session_state.selected_row}:  "
            f"ec_invs = {row['ec_invs']},  EC_coefs = {row['EC_coefs']},  qf_coefs = {row['qf_coefs']}"
        )
        if st.button("View this curve in EC Search →", use_container_width=True):
            st.session_state.ec_prefill = {"f": int(f_sel), "g": int(g_sel), "p": p}
            st.switch_page("pages/2_EC_Search.py")

# ── Tab 2: Isogeny graph ──────────────────────────────────────────────────────
with tab_isograph:
    st.markdown(
        "Enter a prime ℓ to compute the degree-ℓ isogeny graph. "
        "Rows and columns of the adjacency matrix are in the same order as the table above."
    )

    l_input = st.number_input("ℓ (prime)", min_value=2, value=2, step=1)
    compute = st.button("Compute", use_container_width=True)

    if compute:
        l = int(l_input)
        if not primeQ(l):
            st.error(f"{l} is not prime. Please enter a prime ℓ.")
        else:
            with st.spinner(f"Computing degree-{l} isogeny graph…"):
                mat    = isoclass.adjacency_matrix(l)
                labels = [str(qf) for qf in isoclass.qfs_ordered]
                adj_df = pd.DataFrame(mat, index=labels, columns=labels)
            st.session_state["iso_l"]      = l
            st.session_state["iso_adj_df"] = adj_df

    if "iso_adj_df" in st.session_state and st.session_state.get("iso_l") is not None:
        l      = st.session_state["iso_l"]
        adj_df = st.session_state["iso_adj_df"]

        st.subheader(f"Degree-{l} adjacency matrix")
        st.dataframe(adj_df, use_container_width=True)

        st.subheader("Graph picture")

        # ── Label toggle ──────────────────────────────────────────────────────
        label_choice = st.radio(
            "Node labels",
            ["Index (0, 1, …)", "Quadratic form (a, b, c)", "ec_invs"],
            horizontal=True,
            key="iso_label_mode",
        )

        qfs_ordered = list(isoclass.qfs_ordered)
        if label_choice == "Index (0, 1, …)":
            qf_to_label = {qf: str(i) for i, qf in enumerate(qfs_ordered)}
        elif label_choice == "Quadratic form (a, b, c)":
            qf_to_label = {qf: str(qf) for qf in qfs_ordered}
        else:  # ec_invs
            df_local = isoclass.ecqf_df()
            qf_to_label = {
                tuple(row["qf_coefs"]): str(row["ec_invs"])
                for _, row in df_local.iterrows()
            }
            # Fallback for any qf not in df (shouldn't happen, but just in case)
            for qf in qfs_ordered:
                if qf not in qf_to_label:
                    qf_to_label[qf] = str(qfs_ordered.index(qf))

        fig = isogeny_graph_figure(isoclass, l, qf_to_label)
        st.pyplot(fig)
        plt.close(fig)
