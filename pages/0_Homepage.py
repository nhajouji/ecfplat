import streamlit as st
from pathlib import Path

LOGO = Path(__file__).parent / "ec_logo.jpg"

title_col, logo_col = st.columns([3, 1])
with title_col:
    st.markdown("<style>h1 { text-align: center; }</style>", unsafe_allow_html=True)
    st.markdown("# Elliptic Curves over $\\mathbb{F}_p$")
    st.markdown("# via CM Lattices")
with logo_col:
    st.image(str(LOGO), width="stretch")
st.markdown(
    "Tools for visualizing various aspects of categories of elliptic curves over $\\mathbb{F}_p$ "
    "using precomputed equivalences with a category of CM lattices."
)
st.caption("Created by Nadir Hajouji · Work in progress")

st.markdown(
    "The lattice point data produced here is used as input for shader-rendered artwork "
    "by Nadir Hajouji and Steve Trettel. You can view the results at "
    "[elliptic-curves.art](https://elliptic-curves.art/)."
)


st.divider()

col1, col2 = st.columns(2)

with col1:
    st.subheader("Isogeny Class")
    st.markdown(
        "Start from a pair *(a, p)* to load an entire isogeny class. "
        "Browse the bijection table and inspect lattice pictures for each curve."
    )
    st.page_link("pages/1_Isogeny_Class.py", label="Go to Isogeny Class →")

with col2:
    st.subheader("EC Search")
    st.markdown(
        "Start from a specific elliptic curve *y² = x³ + fx + g* over **F**_p. "
        "Look up its trace of Frobenius, associated lattice, and view classical and lattice pictures."
    )
    st.page_link("pages/2_EC_Search.py", label="Go to EC Search →")

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
    "The interactive tools above currently load the $p \\leq 1024$ tables."
)

st.divider()

st.markdown("### Updates")
st.markdown(
    "**June 2026** — Public beta release. "
    "Isogeny class browser, EC search, isogeny graph visualisation, "
    "and Background lessons (algebraic curves, real and $\\mathbb{F}_p$ elliptic curves, "
    "analytic viewpoint, endomorphisms and CM). "
    "Interactive features currently available for primes $p \\leq 1024$."
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
