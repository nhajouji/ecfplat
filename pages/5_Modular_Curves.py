import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "pycode"))

import streamlit as st
import streamlit.components.v1 as components

import modular_viz

st.header("Modular Curves")

st.markdown(
    "This page begins a thread on **modular curves** — the geometry behind how "
    "the project classifies and computes isogenies. We build it up one moduli "
    "space at a time, and always **analytically first**. The starting point, the "
    "one everything else rests on, is the space whose points are the possible "
    "*shapes* of a lattice."
)

st.divider()

st.subheader("The moduli space of lattices")

st.markdown(
    "A **lattice** in $\\mathbb{C}$ is the set of all integer combinations "
    "$\\Lambda = \\mathbb{Z}\\,\\omega_1 + \\mathbb{Z}\\,\\omega_2$ of two "
    "$\\mathbb{R}$-independent complex numbers — an infinite, evenly spaced grid. "
    "What we care about is a lattice's **shape**: not its size, not its "
    "orientation, and not where we happen to put the origin."
)

st.markdown(
    "Two operations leave the shape untouched. We may slide the grid so that one "
    "point sits at $0$, and we may **rotate and rescale** the whole plane — that "
    "is, multiply every point by a single nonzero complex number "
    "$\\zeta \\in \\mathbb{C}^{\\times}$. Using that freedom, divide the basis "
    "through by $\\omega_1$:"
)
st.latex(r"\Lambda \;\sim\; \tfrac{1}{\omega_1}\Lambda \;=\; \mathbb{Z} + \mathbb{Z}\,\tau,"
         r"\qquad \tau = \omega_2/\omega_1.")
st.markdown(
    "Swapping $\\omega_1$ and $\\omega_2$ if necessary so the basis is "
    "positively oriented, $\\tau$ lands in the **upper half-plane** "
    "$\\mathcal{H} = \\{\\tau : \\operatorname{Im}\\tau > 0\\}$. So every lattice, "
    "up to rotation and scaling, is $\\mathbb{Z} + \\mathbb{Z}\\tau$ for a single "
    "$\\tau \\in \\mathcal{H}$ — that one complex number records the shape."
)

st.markdown(
    "But $\\tau$ is not quite unique, because the basis wasn't. Any other "
    "positively-oriented basis of the *same* lattice comes from an integer change "
    "of basis of determinant $1$ — an element "
    "$\\gamma = \\left(\\begin{smallmatrix} a & b \\\\ c & d \\end{smallmatrix}\\right)$ "
    "of $\\mathrm{SL}(2,\\mathbb{Z})$ — and it moves $\\tau$ to"
)
st.latex(r"\gamma \cdot \tau \;=\; \frac{a\tau + b}{c\tau + d}.")
st.markdown(
    "Composing two changes of basis matches multiplying the matrices, so this is "
    "a genuine **left action** of $\\mathrm{SL}(2,\\mathbb{Z})$ on $\\mathcal{H}$; "
    "since $-I$ acts trivially, the group that really acts is "
    "$\\Gamma = \\mathrm{PSL}(2,\\mathbb{Z})$."
)

st.info(
    "**A word on sides.** The action on the *points* of $\\mathcal{H}$ is on the "
    "**left**, which is exactly why the quotient below is written "
    "$\\Gamma \\backslash \\mathcal{H}$ (group on the left). The **right** action "
    "— the one on spaces of modular forms, via the slash operator "
    "$f \\mapsto f|_k\\gamma$ — is a different action we will meet later. You will "
    "often see $\\mathcal{H}/\\Gamma$ written loosely for the same space; we keep "
    "the careful convention."
)

st.markdown(
    "Putting it together, the set of lattice shapes is the orbit space"
)
st.latex(r"\{\text{lattices}\}\,/\,\mathbb{C}^{\times} \;\cong\; \Gamma \backslash \mathcal{H}.")
st.markdown(
    "Every shape has exactly one representative in the standard **fundamental "
    "domain** $\\{\\tau \\in \\mathcal{H} : |\\operatorname{Re}\\tau| \\le "
    "\\tfrac12,\\ |\\tau| \\ge 1\\}$ (shaded in the applet below). This quotient "
    "is our first modular curve. Giving it an honest *coordinate* — the "
    "$j$-invariant, next — will turn it into a concrete object we can compute "
    "with."
)

st.markdown("#### Explore")
st.markdown(
    "Two views of one fact. In **drive the shape τ**, drag $\\tau$ around "
    "$\\Gamma\\backslash\\mathcal{H}$ and watch the lattice $\\mathbb{Z} + "
    "\\mathbb{Z}\\tau$ redraw; the gold ring shows the canonical representative in "
    "the fundamental domain. In **drag a basis**, move $\\omega_1, \\omega_2$ "
    "around $\\mathbb{C}$ and watch their ratio $\\tau = \\omega_2/\\omega_1$ "
    "appear — and reduce — in the domain. Different bases, and shapes outside the "
    "domain, all collapse to the same canonical point of $\\Gamma\\backslash"
    "\\mathcal{H}$."
)
components.html(modular_viz.moduli_applet_html(), height=540, scrolling=False)
