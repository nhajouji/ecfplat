import sys
import json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "pycode"))

import streamlit as st
import streamlit.components.v1 as components
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb

import modular_viz

# First few q-expansion coefficients of j = 1/q + 744 + c(1)q + ...  (enough for
# a domain-coloring picture; in the fundamental domain |q| is tiny so the series
# converges after a handful of terms).
_JC = [196884, 21493760, 864299970, 20245856256, 333202640600,
       4252023300096, 44656994071935, 401490886656000]


@st.cache_data(show_spinner=False)
def _j_domain_coloring(w: int = 440, h: int = 300):
    """RGB domain-coloring of j over a window of ℍ. Each τ is first reduced to
    the fundamental domain (j is Γ-invariant, so j(τ) = j(reduced τ) exactly),
    which both makes the q-series converge and reveals the Γ-tiling."""
    xs = np.linspace(-2.0, 2.0, w)
    ys = np.linspace(2.0, 0.05, h)
    Z = xs[None, :] + 1j * ys[:, None]
    z = Z.copy()
    for _ in range(60):                       # reduce to fundamental domain
        z -= np.round(z.real)
        m = np.abs(z) < 1 - 1e-9
        if not m.any():
            break
        z[m] = -1.0 / z[m]
    q = np.exp(2j * np.pi * z)                 # j via q-expansion
    j = 1.0 / q + 744.0
    qn = q.copy()
    for c in _JC:
        j += c * qn
        qn *= q
    mag = np.abs(j)
    hue = (np.angle(j) / (2 * np.pi)) % 1.0
    val = 0.34 + 0.46 * (0.5 + 0.5 * np.sin(np.log(mag + 1e-9) * 1.7))  # log-|j| bands
    hsv = np.stack([hue, np.full_like(mag, 0.58), np.clip(val, 0, 1)], axis=-1)
    return hsv_to_rgb(hsv)


def _j_coloring_figure():
    rgb = _j_domain_coloring()
    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    ax.imshow(rgb, extent=[-2, 2, 0.05, 2], origin="upper",
              aspect="auto", interpolation="bilinear")
    th = np.linspace(np.pi / 3, 2 * np.pi / 3, 120)   # fundamental-domain outline
    ax.plot(np.cos(th), np.sin(th), color="white", lw=1.1, alpha=0.6)
    ax.plot([-0.5, -0.5], [np.sqrt(3) / 2, 2], color="white", lw=1.1, alpha=0.6)
    ax.plot([0.5, 0.5], [np.sqrt(3) / 2, 2], color="white", lw=1.1, alpha=0.6)
    ax.scatter([0], [1], s=22, color="white", zorder=5)
    ax.annotate("i  (j = 1728)", (0, 1), color="white", fontsize=8,
                xytext=(6, -3), textcoords="offset points")
    ax.scatter([0.5], [np.sqrt(3) / 2], s=22, color="white", zorder=5)
    ax.annotate("ρ  (j = 0)", (0.5, np.sqrt(3) / 2), color="white", fontsize=8,
                xytext=(6, 3), textcoords="offset points")
    ax.set_xlim(-2, 2)
    ax.set_ylim(0.05, 2)
    ax.set_xlabel(r"$\mathrm{Re}\,\tau$")
    ax.set_ylabel(r"$\mathrm{Im}\,\tau$")
    ax.set_title(r"$j$ on $\mathcal{H}$ — hue = phase, bands = $\log|j|$; "
                 r"the pattern repeats over each $\Gamma$-translate", fontsize=9)
    fig.tight_layout()
    return fig

@st.cache_data(show_spinner=False)
def _hilbert_applet_data(dmax: int = 100, coeff_digits: int = 15):
    """CM points + Hilbert class polynomials for the singular-moduli applet.

    Self-contained (loads data/hilbpolys_crt.json, computes reduced forms here)
    so the page never pulls in the heavy bijection modules. Pruned to keep the
    fundamental domain readable: bound |D| (so Im τ stays in view) and drop
    polynomials with monstrous coefficients."""
    from math import gcd, sqrt

    lib_path = Path(__file__).parent.parent / "pycode" / "data" / "hilbpolys_crt.json"
    with open(lib_path) as f:
        lib = json.load(f)

    def reduced_forms(D):                 # primitive reduced (A,B,C), B^2-4AC = D < 0
        forms, B = [], D % 2
        while B * B <= -D / 3 + 1e-9:
            AC = (B * B - D) // 4
            A = max(B, 1)
            while A * A <= AC:
                if AC % A == 0:
                    Cc = AC // A
                    if gcd(gcd(A, B), Cc) == 1:
                        forms.append((A, B, Cc))
                        if B > 0 and B < A and A < Cc:
                            forms.append((A, -B, Cc))
                A += 1
            B += 2
        return forms

    def poly_html(cs):                    # H_D(x) = sum cs[k] x^k, monic, low-to-high
        h, parts = len(cs) - 1, []
        for k in range(h, -1, -1):
            c = cs[k]
            if c == 0:
                continue
            mag = abs(c)
            sign = ("−" if c < 0 else "") if not parts else (" − " if c < 0 else " + ")
            coef = "" if (k > 0 and mag == 1) else str(mag)
            xp = "" if k == 0 else ("x" if k == 1 else f"x<sup>{k}</sup>")
            parts.append(sign + coef + xp)
        return "".join(parts) or "0"

    entries = []
    for ds, cs in lib.items():
        d = int(ds)
        if d >= 0 or -d > dmax:
            continue
        if max(abs(c) for c in cs) >= 10 ** coeff_digits:
            continue
        sq = sqrt(-d)
        pts = [[-B / (2 * A), sq / (2 * A)] for (A, B, C) in reduced_forms(d)]
        entries.append({"d": d, "h": len(cs) - 1, "pts": pts,
                        "poly": poly_html(cs),
                        "j0": (-cs[0] if len(cs) - 1 == 1 else None)})
    entries.sort(key=lambda e: e["d"], reverse=True)
    return entries


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


st.divider()

st.subheader("The $j$-invariant")

st.markdown(
    "The space $\\Gamma\\backslash\\mathcal{H}$ is still an abstract quotient. To "
    "compute with it we need a **coordinate** — a function that hands each shape a "
    "number. That function is the **$j$-invariant**, one of the most storied "
    "objects in mathematics: it runs like a thread from Gauss, Hermite and Klein's "
    "study of elliptic and modular functions, through the theory of complex "
    "multiplication, all the way to Monstrous Moonshine."
)
st.markdown(
    "By construction $j$ is a function on $\\mathcal{H}$ that is **invariant under "
    "$\\Gamma$**: $j(\\gamma\\tau) = j(\\tau)$ for every $\\gamma \\in "
    "\\mathrm{PSL}(2,\\mathbb{Z})$ (weight $0$, so the slash action is just "
    "precomposition). It is the **Hauptmodul** — the *principal modulus* — for "
    "$\\mathrm{SL}(2,\\mathbb{Z})$: every $\\Gamma$-invariant function is a "
    "rational function of $j$."
)
st.caption(
    "A modular *function*, not a modular *form*. Forms like $E_4, E_6, \\Delta$ are "
    "holomorphic everywhere, including at the cusp $i\\infty$; $j$ has a **pole** "
    "there — that is the $1/q$ below. So it is fair to call $j$ the most famous "
    "modular *function*."
)
st.markdown(
    "Being $\\Gamma$-invariant, $j$ is in particular invariant under $\\tau "
    "\\mapsto \\tau + 1$, so it has a Fourier expansion in $q = e^{2\\pi i \\tau}$:"
)
st.latex(r"j(\tau) \;=\; \frac{1}{q} + 744 + 196884\,q + 21493760\,q^2 + 864299970\,q^3 + \cdots")
st.markdown(
    "Those coefficients are one of the most famous sequences in mathematics. The "
    "first, $196884$, is $1 + 196883$ — and $196883$ is the dimension of the "
    "smallest nontrivial irreducible representation of the **Monster**, the "
    "largest sporadic finite simple group. McKay's observation of that "
    "coincidence launched **Monstrous Moonshine**."
)

st.markdown("#### Coloring $\\mathcal{H}$ by $j$")
st.markdown(
    "The picture colors each $\\tau$ by the value $j(\\tau)$ — hue from its phase, "
    "brightness banded by $\\log|j|$. Because $j$ is $\\Gamma$-invariant, the "
    "pattern **repeats over every copy of the fundamental domain**: you are "
    "looking at $\\Gamma\\backslash\\mathcal{H}$ tiled across the plane. The two "
    "corners are the special shapes — $\\tau = i$ (the square lattice, $j = 1728$) "
    "and $\\tau = \\rho = e^{i\\pi/3}$ (the hexagonal lattice, $j = 0$)."
)
st.pyplot(_j_coloring_figure())

st.divider()
st.markdown("#### What $j$ does for us")
st.markdown("For this project, $j$ earns its keep in three ways.")
st.markdown(
    "1. **It separates shapes.** Two lattices are homothetic **iff** they have the "
    "same $j$ — equivalently, $j : \\mathcal{H} \\to \\mathbb{C}$ *is* the quotient "
    "map $\\mathcal{H} \\to \\Gamma\\backslash\\mathcal{H}$, now realized as an "
    "honest coordinate. Knowing $j$ is knowing the shape.\n"
    "2. **It pins down the algebraic curve.** Over any field, the $j$-invariant "
    "determines an elliptic curve up to isomorphism (over the algebraic closure), "
    "and conversely. So $j$ is the dictionary word translating between the "
    "analytic shape and the isomorphism class of the algebraic model.\n"
    "3. **CM shapes have algebraic $j$.** In general $j(\\tau)$ is "
    "**transcendental**, so we cannot write it down exactly. The exception is "
    "decisive: when $\\Lambda$ has **complex multiplication**, $j(\\tau)$ is an "
    "*algebraic integer*."
)
st.markdown(
    "Those exceptional values are the **singular moduli**, and the monic integer "
    "polynomial they satisfy is the **Hilbert class polynomial** $H_D(x)$: its "
    "roots are exactly the $j$-invariants of the CM lattices of discriminant $D$. "
    "This is the doorway from the analytic picture to *exact* computation — and "
    "it is where the whole project's CM machinery begins."
)
st.caption(
    "Historically these are called *singular moduli*; the terminology descends "
    "from the 19th-century theory of elliptic integrals, where the CM values are "
    "the “singular” ones carrying extra symmetry (Kronecker, Weber)."
)

st.markdown("#### Singular moduli and Hilbert class polynomials")
st.markdown(
    "The applet plots the **CM points** — the shapes $\\tau$ with complex "
    "multiplication — inside the fundamental domain, each colour a single "
    "discriminant $D$. The lattices of one discriminant share an endomorphism "
    "ring, hence a single **Hilbert class polynomial** $H_D$; the number of them "
    "is the class number $h(D)$. Click a point to light up its siblings and load "
    "$H_D(x)$. For $h(D) = 1$ there is a lone point and its $j$ is an honest "
    "integer — the singular modulus."
)
components.html(
    modular_viz.hilbert_applet_html(_hilbert_applet_data()),
    height=530, scrolling=False,
)
st.caption(
    "Shown: discriminants down to $|D| \\le 100$ whose class polynomials stay "
    "legible. The library in the repository goes much further — 772 polynomials, "
    "down to $D = -86227$ — but the coefficients grow enormous, so this view is "
    "deliberately pruned by discriminant and coefficient size."
)
