import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "pycode"))

import streamlit as st
import pandas as pd

from ecqf import ECQFIsogenyClass, ap_in_pc_data, ec_eq_str_base
from plotly_tools import fd_points_figure, isogeny_graph_figure
from palette import row_colors
from nt import primeQ


def _plotly_selected_index(event):
    """Extract the selected point's df-row index from a plotly_chart event.

    The selection payload shape varies slightly across versions, so read the
    point index defensively. Points carry customdata == their df-row index.
    """
    sel = getattr(event, "selection", None)
    if not sel:
        return None
    points = sel.get("points") if isinstance(sel, dict) else getattr(sel, "points", None)
    if not points:
        return None
    pt = points[0]
    for key in ("customdata", "point_index", "point_number", "pointIndex", "pointNumber"):
        val = pt.get(key) if isinstance(pt, dict) else getattr(pt, key, None)
        if isinstance(val, (list, tuple)):
            val = val[0] if val else None
        if val is not None:
            return int(val)
    return None



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
    a_input = st.number_input("a", value=prefill["a"] if prefill else 0, step=1)
    p_input = st.number_input("p", value=prefill["p"] if prefill else 307, step=1, min_value=2, max_value=1021)
    load = st.button("Load isogeny class", width="stretch")

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

# Selection events from each widget are gathered across both tabs, then resolved
# once at the end so a click on the table, the FD picture, or a graph node all
# drive the same st.session_state.selected_row.
graph_event = None

# ── Tab 1: DataFrame + fundamental-domain picture ─────────────────────────────
with tab_df:
    st.markdown(
        "Each row pairs an isomorphism class of **CM lattices** with an "
        "associated isomorphism class of **elliptic curves over $\\mathbb{F}_p$**. "
        "The two sides are matched by an equivalence of categories, which lets "
        "us use the lattices — concretely, reduced binary quadratic forms of "
        "discriminant *d = a² − 4p* — to carry out computations about the "
        "curves they correspond to. The picture plots each lattice class as its "
        "CM point *τ* in the fundamental domain of SL(2,ℤ)\\ℍ; the **Unit "
        "disc** view shows the image under *τ ↦ −1/τ*, where the trivial class "
        "sits near the origin. Each class has its own colour. **Click a point "
        "or a row** to select it — the table and picture stay in sync. The "
        "columns are explained at the bottom of the page."
    )
    df = isoclass.ecqf_df()
    sel = st.session_state.selected_row

    # Per-row colours, shared with the FD picture and the isogeny graph.
    colors = row_colors(len(df))

    # Highlight the selected row in the table, and show each row's colour swatch.
    def _highlight_row(s, _sel=sel):
        return [
            "background-color: #ffd6d6" if s.name == _sel else ""
            for _ in s
        ]

    styled = (
        df.style
        .apply(_highlight_row, axis=1)
        .apply(lambda s: [f"color: {colors[s.name]}; font-weight: 900" for _ in s],
               axis=1, subset=["ec_invs"])
        .hide(axis="index")
    )

    col_tbl, col_fig = st.columns([1, 1])

    with col_tbl:
        df_event = st.dataframe(
            styled,
            width="stretch",
            on_select="rerun",
            selection_mode="single-row",
            key="ic_table",
        )

    with col_fig:
        view_label = st.radio(
            "Domain view",
            ["Standard", "Unit disc (τ ↦ −1/τ)"],
            index=1,  # default to the unit disc — always bounded and tidy
            horizontal=True,
            key="ic_fd_view",
            help=(
                "The trivial (principal) class sits high near the cusp and can "
                "fall outside the standard view. The τ ↦ −1/τ view maps it to "
                "the origin and packs every class into the unit disc, so the "
                "picture is bounded and a fixed size."
            ),
        )
        view = "inverted" if view_label.startswith("Unit") else "standard"
        fig = fd_points_figure(df, selected_row=sel, view=view, colors=colors)
        plot_event = st.plotly_chart(
            fig,
            width="stretch",
            on_select="rerun",
            selection_mode="points",
            key=f"ic_fd_plot_{view}",
        )

    # ── Selected-curve summary + navigation ───────────────────────────────────
    if st.session_state.selected_row is not None:
        row = df.iloc[st.session_state.selected_row]
        f_sel, g_sel = tuple(row["EC_coefs"])
        st.success(
            f"Selected row {st.session_state.selected_row}:  "
            f"ec_invs = {row['ec_invs']},  EC_coefs = {row['EC_coefs']},  qf_coefs = {row['qf_coefs']}"
        )
        if st.button("View this curve in EC Search →", width="stretch"):
            st.session_state.ec_prefill = {"f": int(f_sel), "g": int(g_sel), "p": p}
            st.switch_page("pages/2_EC_Search.py")

    # ── Column legend ─────────────────────────────────────────────────────────
    st.divider()
    st.markdown(r"""
**The column labels mean the following:**

* **`ec_invs`** contains a *signature* associated to each row. For ordinary classes (i.e. if $a \neq 0$), `ec_invs` simply contains the $j$-invariant. For supersingular classes, this column contains pairs $(j, \pm 1)$, where $j$ is the $j$-invariant, and the sign $\pm 1$ is equal to $1$ if the constant term is a square (or the constant term is $0$ and the linear term is a square) and $-1$ otherwise. The signature (together with the trace of Frobenius, if the class is ordinary) completely determines the $\mathbb{F}_p$ isomorphism class of each row.
* **`j_inv`**: This column contains the $j$-invariants.
* **`EC_coefs`**: This column contains the coefficients $(c_4, c_6)$ of the short Weierstrass model $y^2 = x^3 + c_4 x + c_6$.
* **`qf_coefs`**: This column contains the coefficients $(a, b, c)$ of the minimal polynomial $a x^2 + b x y + c y^2$ that determines the lattice class.
* **`endo_disc`**: This column contains the discriminant of the endomorphism ring of the object in that row.
* **`endo_cond`**: This contains the conductor of the endomorphism ring.
* **`endo_cocond`**: The integer in this column is the index $[\mathrm{End}(E) : \mathbb{Z}[\phi_E]]$.
* **`frobmat`**: This column contains the matrices that represent the action of the lift of Frobenius on the lattice with respect to the ordered basis $1, \tau$.
* **`tau_s`**: This is a string that represents $\tau$.
* **`tau_xy`**: These are the $x, y$ coordinates of $\tau$.
""")

# ── Tab 2: Isogeny graph ──────────────────────────────────────────────────────
with tab_isograph:
    st.markdown(
        "Enter a prime ℓ to compute the degree-ℓ isogeny graph. "
        "Rows and columns of the adjacency matrix are in the same order as the table above."
    )

    l_input = st.number_input("ℓ (prime)", min_value=2, value=2, step=1)
    compute = st.button("Compute", width="stretch")

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

        # ── Graph picture (shown first) ───────────────────────────────────────
        # Per-class colours, df-row index, and hover text (τ + the curve's
        # equation) — same scheme as the table and the FD picture, so a node
        # matches its row.
        df_local = isoclass.ecqf_df()
        graph_colors = row_colors(len(df_local))
        p_char = isoclass.ap[1]

        def _centered(c, _p=p_char):
            c %= _p
            return c - _p if 2 * c > _p else c

        qf_to_row, qf_to_hover = {}, {}
        for i, (_, row) in enumerate(df_local.iterrows()):
            qf = tuple(row["qf_coefs"])
            qf_to_row[qf] = i
            fg = tuple(_centered(int(v)) for v in row["EC_coefs"])
            qf_to_hover[qf] = f"τ = {row['tau_s']}<br>{ec_eq_str_base(fg)}"

        if isoclass.cond % l == 0:
            st.caption(
                "ℓ ∣ conductor → **volcanoes**: each rim cycle (steelblue, "
                "horizontal isogenies) descends through trees (orange, vertical "
                "isogenies) to the floor."
            )
        else:
            st.caption(
                "ℓ ∤ conductor → disjoint **cycles** drawn as regular n-gons "
                "(steelblue, horizontal isogenies), grouped by endomorphism ring."
            )
        st.caption("Nodes are numbered by table row. Hover for details; click to select (syncs with the table and picture).")

        graph_fig = isogeny_graph_figure(
            isoclass, l, graph_colors, qf_to_row, qf_to_hover,
            selected_row=st.session_state.selected_row,
        )
        graph_event = st.plotly_chart(
            graph_fig, width="stretch",
            on_select="rerun", selection_mode="points",
            key=f"ic_graph_{l}",
        )

        # ── Adjacency matrix (below the picture) ──────────────────────────────
        st.subheader(f"Degree-{l} adjacency matrix")
        st.dataframe(adj_df, width="stretch")


# ── Resolve selection across the table, FD picture, and graph ─────────────────
# Each widget reports its current selection; the Plotly figures are rebuilt fresh
# every run (so they report a selection only on the run of a fresh click), while
# the dataframe persists its selection. To tell a *new* click from a persisted
# one, fire a widget only when its reported selection changed since last run.
def _fired(state_key, current):
    previous = st.session_state.get(state_key)
    st.session_state[state_key] = current
    return current if (current is not None and current != previous) else None

_df_rows = df_event.selection.rows if df_event and df_event.selection else []
_table_sel = _df_rows[0] if _df_rows else None
_fd_sel    = _plotly_selected_index(plot_event)
_graph_sel = _plotly_selected_index(graph_event)

_new = next(
    (s for s in (
        _fired("_sel_prev_graph", _graph_sel),
        _fired("_sel_prev_fd", _fd_sel),
        _fired("_sel_prev_table", _table_sel),
    ) if s is not None),
    None,
)
if _new is not None and _new != st.session_state.selected_row:
    st.session_state.selected_row = _new
    st.rerun()
