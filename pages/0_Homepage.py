import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "pycode"))

import streamlit as st

LOGO = Path(__file__).parent / "ec_logo.jpg"
HERO = Path(__file__).parent / "gallery_img" / "full" / "FrontCover.jpg"

title_col, logo_col = st.columns([3, 1])
with title_col:
    st.markdown("<style>h1 { text-align: center; }</style>", unsafe_allow_html=True)
    st.markdown("# Elliptic Curves")
    st.markdown(
        '<p style="text-align:center; font-size:1.05rem; opacity:.85;">'
        "pictures, tools, and theory — from curves over finite fields to "
        "the lattices that draw them</p>",
        unsafe_allow_html=True,
    )
with logo_col:
    st.image(str(LOGO), width="stretch")

st.caption("Created by Nadir Hajouji and Steve Trettel")

st.divider()

# ── the picture this site is about ────────────────────────────────────────────
st.markdown("### The picture this site is about")
st.image(
    str(HERO),
    caption="Two elliptic curves over finite fields — y² = x³ + 3x (mod 5) "
            "with its 640 points over 𝔽₆₂₅, and y² = x³ + 3 (mod 7) with its "
            "2379 points over 𝔽₂₄₀₁ — lifted to characteristic zero and "
            "drawn on tori in ℝ³.",
    width="stretch",
)

st.markdown(
    "- **[Background](/Background)** — how a picture like this *is* the "
    "classical theory of elliptic curves over $\\mathbb{F}_p$, one "
    "interactive lesson at a time.\n"
    "- **[Explorer](/Explorer)** — the precomputed curve ↔ lattice "
    "equivalences behind it: open any isogeny class and see its curves, "
    "graphs, and lattices live.\n"
    "- **[Gallery](/Gallery)** — the pictures we've made, each with its "
    "story.\n"
    "- **[…and more](/SS_Graph)** — the supersingular isogeny graph, for a "
    "start."
)

# ── coverage ──────────────────────────────────────────────────────────────────
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
    "**July 2026** — Full remodel: one drill-down **Explorer** with canvas "
    "isogeny volcanoes and shareable URLs; the supersingular graph and the "
    "entire Background are 100% canvas; the **Gallery** opens with the ICM "
    "and Bridges renders; the site moves to **elliptic-curves.info**.\n\n"
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
    "[Steve Trettel](https://stevejtrettel.site/), and more of it lives at "
    "[elliptic-curves.art](https://elliptic-curves.art/). "
    "This project grew out of conversations at the conference "
    "*[Integrating Research and Illustration in Number Theory](https://indico.math.cnrs.fr/event/16261/)*, "
    "and I am grateful to the illustrating mathematics community for the inspiration."
)
