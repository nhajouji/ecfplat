"""The Gallery — path-traced renders of lifted elliptic curves.

The images are Nadir Hajouji & Steve Trettel's shader-rendered artwork
(see elliptic-curves.art and the Bridges 2025 paper, arXiv:2505.09627).
The prose here follows this site's narrative — lift Frobenius, read the
volcano — and every artwork links to the Explorer view where the same
curve is alive.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "pycode"))

import streamlit as st

IMG = Path(__file__).parent / "gallery_img"


def _pic(name: str, caption: str):
    st.image(str(IMG / name), caption=caption, width="stretch")


st.header("Gallery")
st.markdown(
    "Everything else on this site draws flat: fundamental domains, lattices, "
    "graphs. The renders below are what those pictures *want* to be — the "
    "torus ℂ/Λ embedded **conformally in ℝ³** (via the Hopf fibration), with "
    "the points of $E(\\mathbb{F}_{p^n})$ living on its surface as the fixed "
    "points of a lifted Frobenius. They are joint work with "
    "[Steve Trettel](https://stevejtrettel.site/): the mathematics of the "
    "lift is the same one the [Explorer](/Explorer) computes, the light is "
    "his. More at [elliptic-curves.art](https://elliptic-curves.art/); the "
    "story in print is [arXiv:2505.09627](https://arxiv.org/abs/2505.09627) "
    "(Bridges 2025)."
)

_pic("FrontCover.jpg",
     "The 640 points of y² = x³ + 3x (mod 5) over 𝔽₆₂₅ (orange), and the "
     "2379 points of y² = x³ + 3 (mod 7) over 𝔽₂₄₀₁ (blue). The marked "
     "point of each curve is highlighted.")

st.divider()

# ── 1. the warm-up ────────────────────────────────────────────────────────────
st.markdown("### 1 · The warm-up: a group you already know")
st.markdown(
    "Before lifting an elliptic curve, lift something simpler: the "
    "multiplicative group $\\mathbb{F}_{p^n}^{\\times}$. Plotting "
    "$\\mathbb{F}_9$ or $\\mathbb{F}_{27}$ as a vector space over "
    "$\\mathbb{F}_p$ (below left) scrambles everything we care about — the "
    "group law and the Galois action are invisible."
)
_pic("Mult_VS.jpg",
     "𝔽₉× and 𝔽₂₇× drawn as vector spaces over 𝔽₃. Blue: the Cayley graph "
     "of the group (ℤ/8 and ℤ/26); green: the Frobenius action. Neither "
     "fits the picture.")
st.markdown(
    "Lift to characteristic zero instead: $\\mathbb{F}_{p^n}^{\\times}$ "
    "becomes the $(p^n{-}1)$-th roots of unity in $\\mathbb{C}^{\\times}$, "
    "the group law becomes rotation, and Frobenius becomes $z \\mapsto z^p$ "
    "— the circle's $p$-fold covering map. One picture serves every $n$."
)
_pic("Mult_Lift.jpg",
     "The same two groups lifted: roots of unity on the circle. Blue: group "
     "structure; green: Frobenius. This is the picture the vector-space "
     "plot was hiding.")
st.caption("Interactive twin: the 𝔽_pⁿ-multiplicative-group applet in "
           "[Background §2](/Background), where you can click a point and "
           "watch its Frobenius orbit.")

st.divider()

# ── 2. lifting Frobenius ──────────────────────────────────────────────────────
st.markdown("### 2 · Lifting Frobenius: y² = x³ + 3x (mod 5)")
st.markdown(
    "The same game for an elliptic curve. Plotting solutions "
    "$(x, y) \\in \\mathbb{F}_{p^n} \\times \\mathbb{F}_{p^n}$ has the "
    "familiar problems — one point off at infinity, the dimension exploding "
    "with $n$:"
)
_pic("VS_Embedding.jpg",
     "y² = x³ + 3x (mod 5) over 𝔽₅ (left) and 𝔽₂₅ (right), plotted as "
     "solution pairs, with the Cayley graphs of ℤ/10 and ℤ/2 × ℤ/10 "
     "struggling to fit.")
st.markdown(
    "Instead, **lift**: the same equation over ℂ is the torus "
    "$\\mathbb{C}/\\Lambda$ for the square lattice $\\Lambda = \\mathbb{Z} "
    "\\oplus i\\mathbb{Z}$, and the Frobenius of the mod-5 curve lifts to an "
    "*endomorphism* of the torus. In Weierstrass coordinates the lift is a "
    "monstrous rational map of degree 5 — but uniformized on the lattice it "
    "collapses to"
)
st.latex(r"\tilde{\varphi}:\; z \;\longmapsto\; (-2+i)\,z \pmod{\Lambda},")
st.markdown(
    "and the points over $\\mathbb{F}_{5^n}$ are simply the fixed points of "
    "$\\tilde{\\varphi}^n$: a scaled-down copy of the lattice, "
    "$\\tfrac{1}{(-2+i)^n - 1}\\Lambda / \\Lambda$. This is exactly the "
    "picture the Explorer's curve view draws — flat there, embedded in ℝ³ "
    "here."
)
_pic("4_Lattice.jpg",
     "The 𝔽₅ points (left pair) and 𝔽₂₅ points (right pair) as scaled "
     "lattices in a fundamental domain of ℂ/(ℤ ⊕ iℤ), with Cayley graphs.")
_pic("4_Torus.jpg",
     "The same points after conformally embedding the torus in ℝ³; the "
     "cutaway shows the Cayley graph wrapping the surface.")
_pic("4_Gallery.jpg",
     "Going up the tower: the fixed points of φ̃³, φ̃⁴, φ̃⁵ are the curve's "
     "points over 𝔽₁₂₅, 𝔽₆₂₅ and 𝔽₃₁₂₅.")
st.markdown(
    "**This exact curve, alive:** "
    "[isogeny class (a, p) = (−4, 5)](/Explorer?a=-4&p=5) · "
    "[discriminant −4 view](/Explorer?d=-4) — the ×α panel there is "
    "multiplication by −2+i, drawn flat."
)

st.divider()

# ── 3. a gallery by discriminant ─────────────────────────────────────────────
st.markdown("### 3 · A gallery, organized the way this site thinks: by discriminant")
st.markdown(
    "Every curve below is a Deuring lift: CM theory hands us the lattice "
    "$\\Lambda$ and the α that Frobenius becomes, and the tower over "
    "$\\mathbb{F}_p$ is α's powers fixing finer and finer copies of "
    "$\\Lambda$. The renders are named by the **discriminant** of "
    "$\\mathbb{Z}[\\alpha]$ — the same label the Explorer files everything "
    "under."
)

st.markdown("#### disc −3 · y² = x³ + 3 (mod 7) — the hexagonal curve")
_pic("3_Intro.jpg",
     "Two conformal embeddings of the hexagonal torus ℂ/(ℤ ⊕ ωℤ), and the "
     "𝔽₄₉ points (≅ ℤ/39) on the fundamental domain and on the surface. "
     "Frobenius lifts to multiplication by −2+ω.")
_pic("3_Gallery.jpg",
     "y² = x³ + 3 (mod 7) over 𝔽₃₄₃, 𝔽₂₄₀₁ and 𝔽₁₆₈₀₇. The sixfold "
     "symmetry of ℤ[ω] — invisible in coordinates — is unmissable on the "
     "torus.")
st.markdown("Live: [(a, p) = (−5, 7)](/Explorer?a=-5&p=7) · "
            "[disc −3 view](/Explorer?d=-3)")

st.markdown("#### disc −7 (as −28 = −7·2²) · y² = x³ + 5x + 7 (mod 11)")
_pic("7_Gallery.jpg",
     "y² = x³ + 5x + 7 (mod 11) over 𝔽₁₂₁, 𝔽₁₃₃₁ and 𝔽₁₄₆₄₁. The Frobenius "
     "order ℤ[α] has discriminant −28; its endomorphism ring sits one level "
     "up the 2-volcano, at field discriminant −7.")
st.markdown("Live: [(a, p) = (−4, 11)](/Explorer?a=-4&p=11) · "
            "[disc −28 view](/Explorer?d=-28) — a two-level volcano")

st.markdown("#### disc −8 · y² = x³ + x + 3 (mod 11)")
_pic("8_Gallery.jpg",
     "y² = x³ + x + 3 (mod 11) over 𝔽₁₂₁, 𝔽₁₃₃₁ and 𝔽₁₄₆₄₁ — the ℤ[√−2] "
     "lattice.")
st.markdown("Live: [(a, p) = (−6, 11)](/Explorer?a=-6&p=11) · "
            "[disc −8 view](/Explorer?d=-8)")

st.markdown("#### disc −11 · y² = x³ + x + 1 (mod 5)")
_pic("11_Gallery.jpg",
     "y² = x³ + x + 1 (mod 5) over 𝔽₁₂₅, 𝔽₆₂₅ and 𝔽₃₁₂₅. Class number 1, "
     "so a single lattice class carries the whole story.")
st.markdown("Live: [(a, p) = (−3, 5)](/Explorer?a=-3&p=5) · "
            "[disc −11 view](/Explorer?d=-11)")

st.markdown("#### Bonus: the real points")
_pic("RealLifts.jpg",
     "The real curves y² = x³ + 3x (red), y² = x³ + 1 (green) and "
     "y² = x³ − x (orange) in the plane, and lifted to the uniformized "
     "complex curve — the real locus is a circle (or two) on the torus, and "
     "the point at infinity is just a point.")

st.divider()

# ── 4. the machine ────────────────────────────────────────────────────────────
st.markdown("### 4 · The machine: Hopf fibration + Pinkall's tori")
st.markdown(
    "Why do these embeddings exist at all? A conformal map "
    "$\\mathbb{C}/\\Lambda \\to \\mathbb{R}^3$ can't come from the usual "
    "torus of revolution (wrong conformal class, except for one lattice). "
    "The trick is the **Hopf fibration** $\\eta: S^3 \\to S^2$: the "
    "preimage of any simple closed curve on $S^2$ is a *flat torus* in "
    "$S^3$, and Pinkall showed every conformal class "
    "$\\mathbb{C}/\\Lambda$ arises this way — the curve's length and "
    "enclosed area read off the lattice. Stereographic projection then "
    "lands it in ℝ³."
)
_pic("hopf.jpeg",
     "The Hopf fibration under stereographic projection: point preimages "
     "are pairwise-linked circles (left); the preimage of a curve on S² is "
     "an annulus (right).")
_pic("HopfTori.jpg",
     "Curves on S² and the flat tori they bound in S³, colored to match "
     "the curves used above. Aesthetic freedom lives here: many curves "
     "realize the same lattice.")

st.divider()

# ── 5. coming next ────────────────────────────────────────────────────────────
st.markdown("### 5 · Coming next")
st.markdown(
    "- **Superlatives** — record-setting primes and first failures, e.g. "
    "**15073**, the smallest prime where every class-number-1 integer "
    "j-invariant is ordinary (the first prime an early algorithm of ours "
    "could not start from).\n"
    "- **New renders** — a wishlist is forming; the hardest part is the "
    "light, and that is Steve's art."
)
