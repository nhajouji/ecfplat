"""The full supersingular ℓ-isogeny graph over F_{p²} — canvas edition.

This used to be a standalone sidebar page; it now lives inside the blog post
that explains it (``render`` is called from ``pages/6_Blog.py``).  State still
lives in query params (?p=..&l=..) so every graph is a shareable URL, matching
the Explorer.  The picture is the canvas widget in ``explorer_viz``: Frobenius
acts literally as the reflection across the horizontal axis.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import streamlit as st
import streamlit.components.v1 as components
import pandas as pd

from nt import primeQ
from ss_graph import ss_isogeny_graph, available_ls, jstr
import explorer_viz

P_MIN, P_CAP = 5, 509      # the graph is recomputed live; keep p small


@st.cache_data(show_spinner=False, max_entries=32)
def _compute(p: int, l: int):
    graph = ss_isogeny_graph(p, l)
    return graph, explorer_viz.ss_graph_descriptor(graph)


def _qp_int(key, default):
    try:
        return int(st.query_params[key])
    except (KeyError, ValueError):
        return default


def render(show_header: bool = True):
    """Render the classical (F_{p²}, (ℓ+1)-regular) supersingular graph widget.

    Pass ``show_header=False`` when embedding under a post's own heading."""
    if show_header:
        st.header("The supersingular ℓ-isogeny graph")
        st.markdown(
            "Vertices are **all** supersingular $j$-invariants in $\\mathbb{F}_{p^2}$ "
            "— not just the $\\mathbb{F}_p$-rational ones the volcano pages see — "
            "with edges counted with multiplicity, so the graph is exactly "
            "$(\\ell+1)$-regular."
        )

    c_p, c_l, c_sp = st.columns([1, 1, 2])
    with c_p:
        p_raw = st.number_input("p (prime)", min_value=P_MIN, max_value=P_CAP,
                                value=_qp_int("p", 101), step=2)
    with c_l:
        l_options = [l for l in available_ls() if l != int(p_raw)]
        l_qp = _qp_int("l", 2)
        l = st.selectbox("ℓ (isogeny degree)", l_options,
                         index=l_options.index(l_qp) if l_qp in l_options else 0)

    p = int(p_raw)
    if not primeQ(p):
        st.error(f"{p} is not prime.")
        st.stop()

    # keep the URL shareable
    if st.query_params.get("p") != str(p) or st.query_params.get("l") != str(l):
        st.query_params["p"] = str(p)
        st.query_params["l"] = str(l)

    with st.spinner("Computing the graph…"):
        graph, desc = _compute(p, int(l))

    verts, rational, adj = graph["vertices"], graph["rational"], graph["adj"]
    n = len(verts)
    n_rat = sum(rational)
    n_pairs = (n - n_rat) // 2
    n_fixed = sum(1 for e in graph["edges"] if e["frob_class"] == "fixed")

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Vertices (ss j in 𝔽_p²)", n)
    m2.metric("𝔽_p-rational vertices", n_rat)
    m3.metric("Conjugate pairs", n_pairs, help="Equals the genus of X₀(p)⁺.")
    m4.metric("Distinct edges", len(graph["edges"]),
              help=f"With multiplicity every vertex has exactly ℓ + 1 = {int(l) + 1} "
                   "outgoing kernels.")

    components.html(explorer_viz.ss_graph_html(desc), height=790, scrolling=False)

    st.markdown(
        "The multiplicity of the edge $j \\to j'$ is the number of order-$\\ell$ "
        "subgroups $C \\subset E_j[\\ell]$ with $j(E_j/C) = j'$ — the multiplicity "
        "of $j'$ as a root of $\\Phi_\\ell(j, Y)$. The layout makes the $p$-power "
        "**Frobenius** literal: it is the reflection across the horizontal axis, "
        "so the gold $\\mathbb{F}_p$-rational vertices lie on the axis and each "
        "conjugate pair $\\{j, j^p\\}$ is mirrored on the circle. The number of "
        "conjugate pairs equals the **genus of $X_0(p)^+$**. Frobenius fixes "
        f"**{n_fixed}** of the **{len(graph['edges'])}** edges (toggle *edges by "
        "Frobenius*: steelblue = fixed, orange = swapped with its mirror image)."
    )

    st.caption(f"The 𝔽_p-rational vertices on the axis are the supersingular "
               f"classes of the [Explorer's (0, {p}) view](/Explorer?a=0&p={p}) — "
               f"note that view sees them up to 𝔽_p-isomorphism (twists), "
               f"this one up to 𝔽̄_p-isomorphism (j-invariants).")

    with st.expander("Adjacency matrix (rows sum to ℓ + 1)"):
        labels = [f"{i}: {jstr(j, graph['s'])}" for i, j in enumerate(verts)]
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
