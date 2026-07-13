import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "pycode"))

import streamlit as st
import streamlit.components.v1 as components

LOGO = Path(__file__).parent / "ec_logo.jpg"

title_col, logo_col = st.columns([3, 1])
with title_col:
    st.markdown("<style>h1 { text-align: center; }</style>", unsafe_allow_html=True)
    st.markdown("# Elliptic Curves over $\\mathbb{F}_p$")
    st.markdown("# via CM Lattices")
with logo_col:
    st.image(str(LOGO), width="stretch")
st.markdown(
    "Tools for visualizing categories of elliptic curves over $\\mathbb{F}_p$ "
    "through precomputed equivalences with categories of CM lattices. "
    "Touch the pictures — the math answers live."
)
st.caption("Created by Nadir Hajouji · Work in progress")

st.divider()

# ── Live showcase: a volcano that drills into the Explorer ────────────────────
st.markdown("### The picture this site is about")
st.markdown(
    "The 2-isogeny graph of the curves over $\\mathbb{F}_{101}$ with trace "
    "$a = 6$ — a **volcano**: the crater cycle is the class group acting on "
    "the maximal order, the trees below descend through the suborders. "
    "Hover a node, click it, and follow the link — you are in the Explorer."
)

from ecqf import ECQFIsogenyClass, class_graph_descriptor
import explorer_viz


@st.cache_resource(show_spinner=False)
def _showcase_html():
    cls = ECQFIsogenyClass(6, 101)
    descs = {l: class_graph_descriptor(cls, l) for l in (2, 3)}
    return explorer_viz.isogeny_graph_html(
        descs, link_base="/Explorer?a=6&p=101", height_px=470)


components.html(_showcase_html(), height=660, scrolling=False)

# ── Entry points ──────────────────────────────────────────────────────────────
st.divider()
c1, c2, c3 = st.columns(3)

with c1:
    st.subheader("Explorer")
    st.markdown(
        "The main journey: pick a **prime** or a **discriminant**, open an "
        "isogeny class, walk its ℓ-isogeny volcano, inspect each curve. "
        "Every view is a shareable URL."
    )
    st.page_link("pages/1_Explorer.py", label="Open the Explorer →")
    st.markdown(
        "Jump straight in: [(a, p) = (6, 101)](/Explorer?a=6&p=101) · "
        "[supersingular (0, 101)](/Explorer?a=0&p=101) · "
        "[discriminant −368](/Explorer?d=-368)"
    )

with c2:
    st.subheader("Background")
    st.markdown(
        "Three interactive chapters: elliptic curves from scratch, the "
        "analytic picture (tori, lattices, CM), and modular curves — "
        "57 canvas applets, no formulas required to start."
    )
    st.page_link("pages/3_Background.py", label="Read the Background →")

with c3:
    st.subheader("Supersingular Graph")
    st.markdown(
        "The full $(\\ell+1)$-regular supersingular graph over "
        "$\\mathbb{F}_{p^2}$, laid out so Frobenius is literally the "
        "reflection across the axis."
    )
    st.page_link("pages/4_SS_Graph.py", label="Open the SS Graph →")

# ── Gallery teaser ────────────────────────────────────────────────────────────
st.divider()
st.markdown("### Gallery — coming soon")
st.markdown(
    "Curated superlatives (first failures, record class numbers, the primes "
    "where the easy algorithms break) and hi-res prints from the physical "
    "gallery. Meanwhile, the shader-rendered artwork built from this site's "
    "lattice data — by Nadir Hajouji and Steve Trettel — lives at "
    "[elliptic-curves.art](https://elliptic-curves.art/)."
)

# ── Coverage ──────────────────────────────────────────────────────────────────
st.divider()
st.markdown("### Dataset coverage")
m1, m2, m3 = st.columns(3)
m1.metric("Ordinary classes computed", "117,155")
m2.metric("Coverage, $4 \\leq p \\leq 8192$", "99.86%")
m3.metric("Supersingular primes", "1,026 / 1,026")
st.caption(
    "The elliptic-curve ↔ CM-lattice bijection has been computed for **117,155** ordinary "
    "isogeny classes with $4 \\leq p \\leq 8192$ — **99.86%** of the range, up from **97.04%** "
    "using modular polynomials alone. The bootstrapping method fills in the harder, "
    "high-class-number discriminants, cutting the blocked classes roughly 20-fold. The "
    "supersingular correspondence is complete for all **1,026** primes in the same range. "
    "The interactive curve tables currently load $p \\leq 1024$; graph structure is exact "
    "for all $p < 8192$."
)

st.divider()

st.markdown("### Updates")
st.markdown(
    "**July 2026** — Full remodel: one drill-down **Explorer** (replacing the "
    "EC Search and Isogeny Class pages) with canvas isogeny volcanoes, "
    "clickable Hasse intervals and shareable URLs; the supersingular graph "
    "and the entire Background are now 100% canvas.\n\n"
    "**June 2026** — Public beta release. Isogeny class browser, EC search, "
    "isogeny graph visualisation, and Background lessons."
)

st.divider()

st.markdown("### Acknowledgements")
st.markdown(
    "This website was built with substantial assistance from "
    "[Claude](https://claude.ai) (Anthropic), which contributed to the codebase, "
    "the app architecture, and the background lessons throughout. "
    "Precomputed data was obtained and verified in part using "
    "[SageMath](https://www.sagemath.org/); the existence of Sage as an open-source "
    "mathematical computing library was part of the inspiration for making this project "
    "publicly accessible.\n\n"
    "The shader-rendered artwork is the work of "
    "[Steve Trettel](https://stevejtrettel.site/). "
    "This project grew out of conversations at the conference "
    "*[Integrating Research and Illustration in Number Theory](https://indico.math.cnrs.fr/event/16261/)*, "
    "and I am grateful to the illustrating mathematics community for the inspiration."
)
