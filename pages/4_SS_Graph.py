import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "pycode"))

import streamlit as st
import pandas as pd

from nt import primeQ
from ss_graph import ss_isogeny_graph, available_ls, jstr
from plotly_tools import ss_graph_figure


def _plotly_selected_index(event):
    """Extract the selected point's vertex index from a plotly_chart event.

    The selection payload shape varies slightly across versions, so read the
    point index defensively. Points carry customdata == their vertex index.
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


@st.cache_data(show_spinner=False, max_entries=32)
def _compute_graph(p: int, l: int):
    return ss_isogeny_graph(p, l)


# ── Session state ─────────────────────────────────────────────────────────────
if "ss_graph" not in st.session_state:
    st.session_state.ss_graph = None
if "ss_selected" not in st.session_state:
    st.session_state.ss_selected = None

# ── Sidebar: (p, ℓ) input ─────────────────────────────────────────────────────
with st.sidebar:
    st.title("Supersingular Graph")
    st.markdown(
        "The **full** supersingular ℓ-isogeny graph in characteristic *p*: "
        "vertices are *all* supersingular *j*-invariants in $\\mathbb{F}_{p^2}$ "
        "(not just the $\\mathbb{F}_p$-rational ones), edges counted with "
        "multiplicity."
    )
    p_input = st.number_input("p (prime)", value=101, step=1,
                              min_value=5, max_value=509)
    l_options = [l for l in available_ls() if l != int(p_input)]
    l_input = st.selectbox("ℓ (isogeny degree)", l_options, index=0)
    load = st.button("Compute graph", width="stretch")

    if load:
        p, l = int(p_input), int(l_input)
        if not primeQ(p):
            st.error(f"{p} is not prime.")
        else:
            try:
                with st.spinner("Computing…"):
                    st.session_state.ss_graph = _compute_graph(p, l)
                    st.session_state.ss_selected = None
                st.success(f"Computed (p, ℓ) = ({p}, {l})")
            except ValueError as exc:
                st.error(str(exc))
                st.session_state.ss_graph = None
                st.session_state.ss_selected = None

# ── Main area ─────────────────────────────────────────────────────────────────
graph = st.session_state.ss_graph

if graph is None:
    st.info("Choose a prime p and an isogeny degree ℓ in the sidebar and click "
            "**Compute graph**.")
    st.stop()

p, l, s = graph["p"], graph["l"], graph["s"]
verts, rational, frob, adj = (graph["vertices"], graph["rational"],
                              graph["frob"], graph["adj"])
n = len(verts)
n_rat = sum(rational)
n_pairs = (n - n_rat) // 2
n_fixed = sum(1 for e in graph["edges"] if e["frob_class"] == "fixed")

st.header(f"Supersingular ℓ-isogeny graph  (p, ℓ) = ({p}, {l})")

c1, c2, c3, c4 = st.columns(4)
c1.metric("Vertices (ss j in 𝔽_p²)", n)
c2.metric("𝔽_p-rational vertices", n_rat)
c3.metric("Conjugate pairs", n_pairs, help="Equals the genus of X₀(p)⁺.")
c4.metric("Distinct edges", len(graph["edges"]),
          help=f"Counted with multiplicity, every vertex has exactly "
               f"ℓ + 1 = {l + 1} outgoing kernels.")

st.markdown(
    "Every supersingular curve is defined over $\\mathbb{F}_{p^2}$, and the "
    "graph is $(\\ell+1)$-regular once edges are counted with multiplicity: "
    "the multiplicity of the edge $j \\to j'$ is the number of order-$\\ell$ "
    "subgroups $C \\subset E_j[\\ell]$ with $j(E_j/C) = j'$, i.e. the "
    "multiplicity of $j'$ as a root of $\\Phi_\\ell(j, Y)$. The layout makes "
    "the $p$-power **Frobenius** literal: it is the reflection across the "
    "horizontal axis, so the gold $\\mathbb{F}_p$-rational vertices — the ones "
    "the volcano pages see — lie on the axis, and each conjugate pair "
    "$\\{j, j^p\\}$ is mirrored on the circle. The number of conjugate pairs "
    "equals the **genus of $X_0(p)^+$**. **Click a vertex** to inspect it."
)

col_a, col_b, col_c = st.columns(3)
with col_a:
    mark_rational = st.checkbox("Mark 𝔽_p-rational vertices", value=True,
                                key="ss_mark_rat")
with col_b:
    show_pairs = st.checkbox("Show Frobenius pairing", value=False,
                             key="ss_show_pairs",
                             help="Dotted connector between each conjugate "
                                  "pair {j, jᵖ}.")
with col_c:
    frob_edges = st.checkbox("Color edges by Frobenius action", value=False,
                             key="ss_frob_edges",
                             help="Steelblue: the edge is fixed by Frobenius. "
                                  "Orange: it is swapped with its mirror "
                                  "image.")

fig = ss_graph_figure(
    graph,
    selected=st.session_state.ss_selected,
    mark_rational=mark_rational,
    show_frob_pairs=show_pairs,
    frob_edge_colors=frob_edges,
)
plot_event = st.plotly_chart(
    fig,
    width="stretch",
    on_select="rerun",
    selection_mode="points",
    key=f"ss_plot_{p}_{l}",
)

if frob_edges:
    st.caption(
        f"Frobenius fixes {n_fixed} of the {len(graph['edges'])} edges "
        "(steelblue) and swaps the rest in mirror pairs (orange). An edge "
        "with both ends on the axis, or joining a conjugate pair, is "
        "automatically fixed."
    )

# ── Selected-vertex summary ──────────────────────────────────────────────────
sel = st.session_state.ss_selected
if sel is not None and 0 <= sel < n:
    nbrs = ", ".join(
        f"#{k} ({jstr(verts[k], s)})" + (f" ×{adj[sel][k]}" if adj[sel][k] > 1 else "")
        for k in range(n) if adj[sel][k]
    )
    tag = ("𝔽_p-rational" if rational[sel]
           else f"Frobenius-conjugate of #{frob[sel]} = {jstr(verts[frob[sel]], s)}")
    st.success(
        f"Selected vertex #{sel}:  j = {jstr(verts[sel], s)}  ({tag}).  "
        f"ℓ-isogenous to: {nbrs}"
    )

# ── Adjacency matrix ──────────────────────────────────────────────────────────
with st.expander("Adjacency matrix (rows sum to ℓ + 1)"):
    labels = [f"{i}: {jstr(j, s)}" for i, j in enumerate(verts)]
    st.dataframe(pd.DataFrame(adj, index=labels, columns=labels),
                 width="stretch")
    st.caption(
        "The matrix is symmetric except at j = 0 and j = 1728, where the "
        "extra automorphisms fold kernels together: Eichler's weighted "
        "symmetry e(j′)·m(j → j′) = e(j)·m(j′ → j) holds with "
        "e = |Aut(E)|/2 (3 at j = 0, 2 at j = 1728, else 1). On the picture "
        "an asymmetric edge shows one ×m per end — the number of kernels "
        "leaving that endpoint."
    )

# ── Resolve selection from the plot ───────────────────────────────────────────
# The figure is rebuilt fresh every run, so it reports a selection only on the
# run of a fresh click; fire only when the reported selection changed since
# the last run (same pattern as the Isogeny Class page).
_sel_now = _plotly_selected_index(plot_event)
_prev = st.session_state.get("_ss_sel_prev")
st.session_state["_ss_sel_prev"] = _sel_now
if _sel_now is not None and _sel_now != _prev and _sel_now != st.session_state.ss_selected:
    st.session_state.ss_selected = _sel_now
    st.rerun()
