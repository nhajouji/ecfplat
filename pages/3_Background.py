import sys
import json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "pycode"))

import streamlit as st
import streamlit.components.v1 as components
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.patches import Arc as MplArc

import slide_viz
import modular_viz
import basics_viz


# ── Chapter 3 (Modular Curves) helpers ────────────────────────────────────────
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


def _click_index(event):
    """df-row / point index of a clicked Plotly point (customdata), or None."""
    sel = getattr(event, "selection", None)
    if not sel:
        return None
    points = sel.get("points") if isinstance(sel, dict) else getattr(sel, "points", None)
    if not points:
        return None
    cd = points[0].get("customdata") if isinstance(points[0], dict) else None
    if isinstance(cd, (list, tuple)):
        cd = cd[0] if cd else None
    return int(cd) if cd is not None else None


def _two_point_pick(prefix, clicked_pt):
    """Maintain a ≤2-point click selection in session_state.

    Clicks fill P then Q; clicking the same point twice gives P = Q (so the
    doubling / tangent case is reachable). Once two points are chosen, a third
    click clears the selection (a quick reset without the Clear button). Any
    change clears the display mode ('line' / 'sum').
    """
    sel = st.session_state.setdefault(prefix + "_sel", [])
    if clicked_pt is not None:
        if len(sel) >= 2:
            sel = []                       # third click clears the pair
        else:
            sel.append(clicked_pt)
        st.session_state[prefix + "_mode"] = None
        st.session_state[prefix + "_sel"] = sel
    return sel


_OINF = "O"   # sentinel for the point at infinity 𝒪 in the click selections



st.header("Background")
st.markdown(
    "A crash course on the mathematics underlying this project. "
    "Each section pairs a short conceptual overview with an interactive applet."
)

st.header("Visualizing Elliptic Curves")

with st.expander("§1 — Elliptic curves: the basics", expanded=False):
    (tab0, tab1, tab2, tab3) = st.tabs([
        "Algebraic curves over ℝ",
        "Real Elliptic Curves",
        "Elliptic curves over $\\mathbb{F}_p$",
        "Elliptic Curves over ℂ",
    ])

with st.expander("§2 — CM Lattices and Frobenius", expanded=False):
    (tab4, tab_frob, tab_lift) = st.tabs([
        "Endomorphisms and Complex Multiplication",
        "The Frobenius endomorphism over $\\mathbb{F}_p$",
        "Lifting Frobenius (Deuring)",
    ])

with st.expander("§3 — Pictures from lifts of Frobenius", expanded=False):
    (tab_mult, tab_ell) = st.tabs([
        "The recipe, via the multiplicative group",
        "Elliptic curves",
    ])

st.header("Equivalence of Categories")

with st.expander("§4 — Isogenies", expanded=False):
    (tab_isog, tab_fold, tab_velu) = st.tabs([
        "Isogenies: kernels and degree",
        "Analytic pictures: folding the torus",
        "Vélu's formulas over $\\mathbb{F}_p$",
    ])

with st.expander("§5 — Modular tools", expanded=False):
    (tab_jfun, tab_modpoly) = st.tabs([
        "The $j$-function and modular curves",
        "Modular polynomials $\\Phi_\\ell$",
    ])

with st.expander("§6 — The CM bijection", expanded=False):
    (tab_setup, tab_rigid, tab_algo) = st.tabs([
        "The gallery problem and the equivalence",
        "Isogeny graphs and rigidity",
        "The algorithm and the gallery",
    ])


st.header("Modular Curves")
st.markdown(
    "The geometry behind how the project classifies and computes isogenies, built "
    "one moduli space at a time and always **analytically first**. Where "
    "Chapters 1–2 draw the pictures and establish the equivalence of categories, "
    "this chapter is the machinery underneath: the modular curves whose points "
    "*are* the objects we have been labelling."
)

with st.expander("§7 — The moduli space of elliptic curves", expanded=False):
    tab_lat, tab_jinv = st.tabs(["The moduli space of lattices", "The $j$-invariant"])

    # ── §7.1 — the moduli space of lattices ──────────────────────────────────
    with tab_lat:
        st.markdown(
            "A **lattice** in $\\mathbb{C}$ is the set of all integer combinations "
            "$\\Lambda = \\mathbb{Z}\\,\\omega_1 + \\mathbb{Z}\\,\\omega_2$ of two "
            "$\\mathbb{R}$-independent complex numbers — an infinite, evenly spaced "
            "grid. What we care about is a lattice's **shape**: not its size, not its "
            "orientation, and not where we happen to put the origin."
        )
        st.markdown(
            "Two operations leave the shape untouched. We may slide the grid so that "
            "one point sits at $0$, and we may **rotate and rescale** the whole plane "
            "— that is, multiply every point by a single nonzero complex number "
            "$\\zeta \\in \\mathbb{C}^{\\times}$. Using that freedom, divide the basis "
            "through by $\\omega_1$:"
        )
        st.latex(r"\Lambda \;\sim\; \tfrac{1}{\omega_1}\Lambda \;=\; \mathbb{Z} + \mathbb{Z}\,\tau,"
                 r"\qquad \tau = \omega_2/\omega_1.")
        st.markdown(
            "Swapping $\\omega_1$ and $\\omega_2$ if necessary so the basis is "
            "positively oriented, $\\tau$ lands in the **upper half-plane** "
            "$\\mathcal{H} = \\{\\tau : \\operatorname{Im}\\tau > 0\\}$. So every "
            "lattice, up to rotation and scaling, is $\\mathbb{Z} + \\mathbb{Z}\\tau$ "
            "for a single $\\tau \\in \\mathcal{H}$ — that one complex number records "
            "the shape."
        )
        st.markdown(
            "But $\\tau$ is not quite unique, because the basis wasn't. Any other "
            "positively-oriented basis of the *same* lattice comes from an integer "
            "change of basis of determinant $1$ — an element "
            "$\\gamma = \\left(\\begin{smallmatrix} a & b \\\\ c & d \\end{smallmatrix}\\right)$ "
            "of $\\mathrm{SL}(2,\\mathbb{Z})$ — and it moves $\\tau$ to"
        )
        st.latex(r"\gamma \cdot \tau \;=\; \frac{a\tau + b}{c\tau + d}.")
        st.markdown(
            "Composing two changes of basis matches multiplying the matrices, so this "
            "is a genuine **left action** of $\\mathrm{SL}(2,\\mathbb{Z})$ on "
            "$\\mathcal{H}$; since $-I$ acts trivially, the group that really acts is "
            "$\\Gamma = \\mathrm{PSL}(2,\\mathbb{Z})$."
        )
        st.info(
            "**A word on sides.** The action on the *points* of $\\mathcal{H}$ is on "
            "the **left**, which is exactly why the quotient below is written "
            "$\\Gamma \\backslash \\mathcal{H}$ (group on the left). The **right** "
            "action — the one on spaces of modular forms, via the slash operator "
            "$f \\mapsto f|_k\\gamma$ — is a different action we will meet later. You "
            "will often see $\\mathcal{H}/\\Gamma$ written loosely for the same space; "
            "we keep the careful convention."
        )
        st.markdown("Putting it together, the set of lattice shapes is the orbit space")
        st.latex(r"\{\text{lattices}\}\,/\,\mathbb{C}^{\times} \;\cong\; \Gamma \backslash \mathcal{H}.")
        st.markdown(
            "Every shape has exactly one representative in the standard **fundamental "
            "domain** $\\{\\tau \\in \\mathcal{H} : |\\operatorname{Re}\\tau| \\le "
            "\\tfrac12,\\ |\\tau| \\ge 1\\}$ (shaded in the applet below). This "
            "quotient is our first modular curve. Giving it an honest *coordinate* — "
            "the $j$-invariant, next — will turn it into a concrete object we can "
            "compute with."
        )
        st.markdown("#### Explore")
        st.markdown(
            "Two views of one fact. In **drive the shape τ**, drag $\\tau$ around "
            "$\\Gamma\\backslash\\mathcal{H}$ and watch the lattice $\\mathbb{Z} + "
            "\\mathbb{Z}\\tau$ redraw; the gold ring shows the canonical "
            "representative in the fundamental domain. In **drag a basis**, move "
            "$\\omega_1, \\omega_2$ around $\\mathbb{C}$ and watch their ratio "
            "$\\tau = \\omega_2/\\omega_1$ appear — and reduce — in the domain. "
            "Different bases, and shapes outside the domain, all collapse to the same "
            "canonical point of $\\Gamma\\backslash\\mathcal{H}$."
        )
        components.html(modular_viz.moduli_applet_html(), height=540, scrolling=False)

    # ── §7.2 — the j-invariant ───────────────────────────────────────────────
    with tab_jinv:
        st.markdown(
            "The space $\\Gamma\\backslash\\mathcal{H}$ is still an abstract quotient. "
            "To compute with it we need a **coordinate** — a function that hands each "
            "shape a number. That function is the **$j$-invariant**, one of the most "
            "storied objects in mathematics: it runs like a thread from Gauss, Hermite "
            "and Klein's study of elliptic and modular functions, through the theory "
            "of complex multiplication, all the way to Monstrous Moonshine."
        )
        st.markdown(
            "By construction $j$ is a function on $\\mathcal{H}$ that is **invariant "
            "under $\\Gamma$**: $j(\\gamma\\tau) = j(\\tau)$ for every $\\gamma \\in "
            "\\mathrm{PSL}(2,\\mathbb{Z})$ (weight $0$, so the slash action is just "
            "precomposition). It is the **Hauptmodul** — the *principal modulus* — for "
            "$\\mathrm{SL}(2,\\mathbb{Z})$: every $\\Gamma$-invariant function is a "
            "rational function of $j$."
        )
        st.caption(
            "A modular *function*, not a modular *form*. Forms like $E_4, E_6, "
            "\\Delta$ are holomorphic everywhere, including at the cusp $i\\infty$; "
            "$j$ has a **pole** there — that is the $1/q$ below. So it is fair to call "
            "$j$ the most famous modular *function*."
        )
        st.markdown(
            "Being $\\Gamma$-invariant, $j$ is in particular invariant under $\\tau "
            "\\mapsto \\tau + 1$, so it has a Fourier expansion in "
            "$q = e^{2\\pi i \\tau}$:"
        )
        st.latex(r"j(\tau) \;=\; \frac{1}{q} + 744 + 196884\,q + 21493760\,q^2 + 864299970\,q^3 + \cdots")
        st.markdown(
            "Those coefficients are one of the most famous sequences in mathematics. "
            "The first, $196884$, is $1 + 196883$ — and $196883$ is the dimension of "
            "the smallest nontrivial irreducible representation of the **Monster**, "
            "the largest sporadic finite simple group. McKay's observation of that "
            "coincidence launched **Monstrous Moonshine**."
        )
        st.markdown("#### Coloring $\\mathcal{H}$ by $j$")
        st.markdown(
            "The picture colors each $\\tau$ by the value $j(\\tau)$ — hue from its "
            "phase, brightness banded by $\\log|j|$. Because $j$ is $\\Gamma$-"
            "invariant, the pattern **repeats over every copy of the fundamental "
            "domain**: you are looking at $\\Gamma\\backslash\\mathcal{H}$ tiled "
            "across the plane. The two corners are the special shapes — $\\tau = i$ "
            "(the square lattice, $j = 1728$) and $\\tau = \\rho = e^{i\\pi/3}$ (the "
            "hexagonal lattice, $j = 0$). **Drag $\\tau$** to read off $j(\\tau)$; "
            "the gold ring is the canonical representative of $\\tau$ in the "
            "fundamental domain — same colour, same $j$."
        )
        components.html(modular_viz.j_coloring_html(), height=475, scrolling=False)
        st.divider()
        st.markdown("#### What $j$ does for us")
        st.markdown("For this project, $j$ earns its keep in three ways.")
        st.markdown(
            "1. **It separates shapes.** Two lattices are homothetic **iff** they "
            "have the same $j$ — equivalently, $j : \\mathcal{H} \\to \\mathbb{C}$ "
            "*is* the quotient map $\\mathcal{H} \\to \\Gamma\\backslash\\mathcal{H}$, "
            "now realized as an honest coordinate. Knowing $j$ is knowing the shape.\n"
            "2. **It pins down the algebraic curve.** Over any field, the "
            "$j$-invariant determines an elliptic curve up to isomorphism (over the "
            "algebraic closure), and conversely. So $j$ is the dictionary word "
            "translating between the analytic shape and the isomorphism class of the "
            "algebraic model.\n"
            "3. **CM shapes have algebraic $j$.** In general $j(\\tau)$ is "
            "**transcendental**, so we cannot write it down exactly. The exception is "
            "decisive: when $\\Lambda$ has **complex multiplication**, $j(\\tau)$ is "
            "an *algebraic integer*."
        )
        st.markdown(
            "Those exceptional values are the **singular moduli**, and the monic "
            "integer polynomial they satisfy is the **Hilbert class polynomial** "
            "$H_D(x)$: its roots are exactly the $j$-invariants of the CM lattices of "
            "discriminant $D$. This is the doorway from the analytic picture to "
            "*exact* computation — and it is where the whole project's CM machinery "
            "begins."
        )
        st.caption(
            "Historically these are called *singular moduli*; the terminology "
            "descends from the 19th-century theory of elliptic integrals, where the "
            "CM values are the “singular” ones carrying extra symmetry (Kronecker, "
            "Weber)."
        )
        st.markdown("#### Singular moduli and Hilbert class polynomials")
        st.markdown(
            "The applet plots the **CM points** — the shapes $\\tau$ with complex "
            "multiplication — inside the fundamental domain, each colour a single "
            "discriminant $D$. The lattices of one discriminant share an endomorphism "
            "ring, hence a single **Hilbert class polynomial** $H_D$; the number of "
            "them is the class number $h(D)$. Click a point to light up its siblings "
            "and load $H_D(x)$. For $h(D) = 1$ there is a lone point and its $j$ is an "
            "honest integer — the singular modulus."
        )
        components.html(
            modular_viz.hilbert_applet_html(_hilbert_applet_data()),
            height=530, scrolling=False,
        )
        st.caption(
            "Shown: discriminants down to $|D| \\le 100$ whose class polynomials stay "
            "legible. The library in the repository goes much further — 772 "
            "polynomials, down to $D = -86227$ — but the coefficients grow enormous, "
            "so this view is deliberately pruned by discriminant and coefficient size."
        )

with st.expander("§8 — The modular curve $X_0(\\ell)$", expanded=False):
    st.markdown(
        "In §4 an isogeny was determined, up to isomorphism, by its **kernel**. "
        "Here we turn that fact into geometry: a modular curve whose points *are* "
        "the cyclic $\\ell$-isogenies. We keep writing $X(1) = \\Gamma\\backslash"
        "\\mathcal{H}$ for the modular curve of §7 (level $1$; its coordinate "
        "$j$ makes it the $j$-line). Throughout, $\\ell$ is prime and we work "
        "analytically first."
    )
    (tab_x0a, tab_x0j, tab_x0alg) = st.tabs([
        "The analytic model",
        "The $j$-map and Fricke involution",
        "Algebraic models",
    ])

    # ── §8.1 — the analytic model ─────────────────────────────────────────────
    with tab_x0a:
        st.subheader("The Analytic Model of $X_0(\\ell)$")
        st.markdown(
            "Fix $E = \\mathbb{C}/\\Lambda$. A cyclic subgroup $C \\subset E$ of "
            "order $\\ell$ pulls back to a **superlattice** $\\Lambda' \\supset "
            "\\Lambda$ with $\\Lambda'/\\Lambda = C$ cyclic of order $\\ell$; the "
            "isogeny is the quotient map $\\mathbb{C}/\\Lambda \\to \\mathbb{C}/"
            "\\Lambda'$. So a cyclic $\\ell$-isogeny out of $E$ is the same data "
            "as an **index-$\\ell$ superlattice** of $\\Lambda$."
        )
        st.markdown(
            "For prime $\\ell$ we can choose bases of both lattices at once: there "
            "is always a basis $\\omega_1, \\omega_2$ of $\\Lambda$ for which "
            "$\\tfrac{1}{\\ell}\\omega_1, \\omega_2$ is a basis of $\\Lambda'$. "
            "Dividing through by $\\omega_1$ as in §7,"
        )
        st.latex(r"\Lambda \sim \langle 1, \tau\rangle \;\subset\; "
                 r"\langle \tfrac{1}{\ell},\, \tau\rangle \sim \Lambda',"
                 r"\qquad \tau = \omega_2/\omega_1 \in \mathcal{H}.")
        st.markdown(
            "Now $\\tau$ records the **whole inclusion** — the domain *and* its "
            "chosen order-$\\ell$ subgroup $C = \\langle\\tfrac{1}{\\ell}\\rangle$ "
            "— not just the shape of $E$."
        )
        st.markdown(
            "How much freedom is left? Exactly the changes of basis of $\\Lambda$ "
            "that preserve the flag, i.e. that fix the subgroup $C$. A matrix "
            "$\\left(\\begin{smallmatrix} a & b \\\\ c & d\\end{smallmatrix}"
            "\\right) \\in \\mathrm{SL}(2,\\mathbb{Z})$ does so precisely when "
            "$c \\equiv 0 \\pmod{\\ell}$ — which *defines*"
        )
        st.latex(r"\Gamma_0(\ell) = \left\{ \left(\begin{smallmatrix} a & b \\ "
                 r"c & d\end{smallmatrix}\right) \in \mathrm{SL}(2,\mathbb{Z}) : "
                 r"c \equiv 0 \pmod{\ell} \right\}.")
        st.markdown(
            "The definition falls out of the geometry rather than being imposed. "
            "The moduli space of cyclic $\\ell$-isogenies is therefore"
        )
        st.latex(r"X_0(\ell) = \Gamma_0(\ell)\backslash\mathcal{H}.")
        st.markdown(
            "How big is it over $X(1)$? The $\\ell$-torsion $E[\\ell] \\cong "
            "(\\mathbb{Z}/\\ell)^2$ has exactly $\\ell + 1$ cyclic subgroups of "
            "order $\\ell$ — the lines in $\\mathbb{F}_\\ell^{\\,2}$ — matching the "
            "index $[\\Gamma(1) : \\Gamma_0(\\ell)] = \\ell + 1$. Choosing "
            "$\\ell+1$ coset representatives $\\gamma_0, \\dots, \\gamma_\\ell$, a "
            "fundamental domain for $X_0(\\ell)$ is the union of the $\\ell+1$ "
            "translates $\\gamma_i \\cdot \\mathcal{F}$ of the level-$1$ domain "
            "$\\mathcal{F}$ — a connected region, symmetric about the imaginary "
            "axis when $\\ell$ is odd."
        )
        st.markdown(
            "The applet makes the dictionary concrete. A point of $X_0(\\ell)$ is "
            "**a curve $E$ plus one of its $\\ell+1$ order-$\\ell$ subgroups**. Drag "
            "the marker across the $\\ell+1$ tiles of the $\\Gamma_0(\\ell)$ "
            "fundamental domain on the left: each tile is one subgroup, so a single "
            "point fixes the curve *and* the subgroup at once. On the right, "
            "$E = \\mathbb{C}/\\langle 1, \\tau\\rangle$ shows its $\\ell^2$ points "
            "of order dividing $\\ell$, coloured by which cyclic subgroup they "
            "generate — you can also click a point there to pick its subgroup "
            "directly, and the whole subgroup lights up."
        )
        components.html(modular_viz.x0_subgroup_html(), height=470, scrolling=False)

    # ── §8.2 — the j-map and Fricke involution ────────────────────────────────
    with tab_x0j:
        st.subheader("The $j$-map $j_\\ell$ and the Fricke Involution $\\mathfrak{F}_\\ell$")
        st.markdown(
            "A point of $X_0(\\ell)$ is a whole isogeny $[E \\to E']$. Two "
            "functions pull the arithmetic back out of it."
        )
        st.markdown(
            "**The $j$-map $j_\\ell$.** Forgetting the subgroup is the map"
        )
        st.latex(r"j_\ell : X_0(\ell) \to X(1), \qquad j_\ell\big([E \to E']\big) = [E].")
        st.markdown(
            "Since $X(1)$ is the $j$-line (§7), $j_\\ell$ records the **domain**'s "
            "$j$-invariant; in the $\\tau$-model it is simply $\\tau \\mapsto "
            "j(\\tau)$. (We keep the subscript because on the algebraic models "
            "there are several maps down to $X(1)$ and we must say which one.)"
        )
        st.markdown(
            "**The Fricke involution $\\mathfrak{F}_\\ell$.** Sending an isogeny to "
            "its dual gives an involution $\\mathfrak{F}_\\ell : X_0(\\ell) \\to "
            "X_0(\\ell)$, $[E \\to E'] \\mapsto [E' \\to E]$. On $\\mathcal{H}$ it "
            "is the Atkin–Lehner matrix $\\left(\\begin{smallmatrix} 0 & -1 \\\\ "
            "\\ell & 0\\end{smallmatrix}\\right)$, i.e."
        )
        st.latex(r"\mathfrak{F}_\ell : \tau \mapsto -\frac{1}{\ell\,\tau}.")
        st.markdown(
            "It normalizes $\\Gamma_0(\\ell)$ (so it descends to $X_0(\\ell)$), "
            "squares to the identity ($\\mathfrak{F}_\\ell^2\\tau = \\tau$), and "
            "carries the flag $\\langle 1,\\tau\\rangle \\subset "
            "\\langle\\tfrac1\\ell,\\tau\\rangle$ to the dual flag."
        )
        st.markdown(
            "Now the payoff. The codomain is $E' = \\mathbb{C}/\\Lambda'$ with "
            "$\\Lambda' = \\langle\\tfrac1\\ell,\\tau\\rangle \\sim \\langle 1, "
            "\\ell\\tau\\rangle$, so its $j$-invariant is $j(\\ell\\tau)$. And "
            "since $j(-1/z) = j(z)$,"
        )
        st.latex(r"j_\ell\big(\mathfrak{F}_\ell\, \tau\big) = j\!\left(-\tfrac{1}{\ell\tau}\right) "
                 r"= j(\ell\tau) = j(E').")
        st.markdown(
            "So $j_\\ell$ and $j_\\ell \\circ \\mathfrak{F}_\\ell$ hand you the "
            "**two endpoints** of the isogeny — the domain and codomain "
            "$j$-invariants. This is exactly the pair $(j, j')$ that the modular "
            "polynomial $\\Phi_\\ell$ relates (§5). With *any* model of "
            "$X_0(\\ell)$ plus these two functions, a point tells you both."
        )
        st.markdown("#### The two $j$-maps for a $2$-isogeny")
        st.markdown(
            "The same domain-colouring as in §7, now for both endpoints. On the "
            "left, $\\tau \\mapsto j_2(\\tau) = j(\\tau)$ colours each point of "
            "$X_0(2)$ by its **domain**; on the right, $\\tau \\mapsto "
            "j_2(\\mathfrak{F}_2\\tau) = j(2\\tau)$ colours it by its **codomain**. "
            "The involution $\\mathfrak{F}_2$ interchanges the two pictures: drag "
            "$\\tau$ and compare it with its gold partner $\\mathfrak{F}_2\\tau$ — "
            "the left colour at one is the right colour at the other."
        )
        components.html(modular_viz.j_fricke_pair_html(), height=445, scrolling=False)

        st.markdown(
            "**Where the two agree.** The isogeny is an endomorphism ($E \\cong "
            "E'$) exactly where $j_\\ell(\\tau) = j_\\ell(\\mathfrak{F}_\\ell"
            "\\tau)$ — the $\\mathfrak{F}_\\ell$-fixed points, and the swapped "
            "pairs collapsing to a single $j$-value. These are the **CM points**: "
            "the singular moduli of §7, precisely the diagonal $\\Phi_\\ell(X, X) "
            "= 0$, and they lie on the circle $|\\tau| = 1/\\sqrt{\\ell}$ (the "
            "fixed locus of $\\mathfrak{F}_\\ell$). So $X_0(\\ell)$ sees the CM "
            "world again — the thread we pull to build the bijection of §6."
        )
        st.markdown(
            "Drag the point $\\tau$ around $X_0(\\ell)$ on the left. On the right, "
            "the isogeny it names is drawn as the two lattices $\\Lambda = \\langle "
            "1, \\tau\\rangle \\subset \\Lambda'$ — the domain in blue, the "
            "index-$\\ell$ superlattice (the codomain $E' = \\mathbb{C}/\\Lambda'$) "
            "in gold. The gold point $\\mathfrak{F}_\\ell(\\tau) = -1/(\\ell\\tau)$ "
            "is the **dual** isogeny $E' \\to E$. The **dashed curves** are the "
            "interior part of the real locus of $j$ — the images of "
            "$\\operatorname{Re}\\tau = 0$ under $\\Gamma_0(\\ell)$ — while the "
            "tile boundaries carry the rest of it (the arcs where $j \\in [0, "
            "1728]$ and the walls where $j \\le 0$)."
        )
        components.html(modular_viz.x0_fricke_html(), height=470, scrolling=False)

    # ── §8.3 — algebraic models ───────────────────────────────────────────────
    with tab_x0alg:
        st.subheader("Algebraic Models")
        st.markdown(
            "The analytic model says *what* $X_0(\\ell)$ is; to compute over "
            "$\\mathbb{F}_p$ we need it as an **algebraic curve** with explicit "
            "equations. The clean route runs through universal families."
        )
        st.markdown(
            "**From the universal curve to $X_1(\\ell)$.** Take the universal "
            "elliptic curve carrying a marked point $P$ — a Weierstrass family "
            "with a few parameters. Imposing “$P$ has exact order $\\ell$” is a "
            "polynomial condition on those parameters, and the relation it cuts "
            "out is a model of $X_1(\\ell)$, the curve parametrizing pairs "
            "$(E, P)$ with $P$ of order $\\ell$. (Smooth for small $\\ell$; we "
            "won't dwell on the resolution needed for larger $\\ell$.)"
        )
        st.markdown(
            "**Down to $X_0(\\ell)$.** A point of $X_0(\\ell)$ remembers only the "
            "*subgroup* $C = \\langle P\\rangle$, not the chosen generator, so"
        )
        st.latex(r"X_0(\ell) = X_1(\ell)\big/\langle m\rangle, \qquad "
                 r"(E, P) \mapsto (E, mP),")
        st.markdown(
            "where $m$ generates $(\\mathbb{Z}/\\ell)^\\times/\\{\\pm 1\\}$ (the "
            "diamond operator)."
        )
        st.markdown(
            "**The maps, explicitly.** Because we hold the universal curve in "
            "hand, $j$ is an explicit rational function of the parameters; it "
            "descends to the $j$-map $j_\\ell : X_0(\\ell) \\to X(1)$. For the "
            "*other* endpoint, apply **Vélu's formulas** (§4) to the universal "
            "curve to get the universal degree-$\\ell$ isogeny $E \\to E/C$; the "
            "$j$-invariant of its codomain is $j_\\ell \\circ \\mathfrak{F}_\\ell$ "
            "as a rational function. (Computing this universal isogeny can even "
            "shortcut the quotient above.)"
        )
        st.markdown(
            "**Worked example: $\\ell = 5$.** Here $X_0(5)$ has genus $0$, so a "
            "single **Hauptmodul** $t$ coordinatizes it and both maps are rational "
            "in $t$:"
        )
        st.latex(r"j_5 = \frac{(t^2 + 10\,t + 5)^3}{t}, \qquad "
                 r"j_5 \circ \mathfrak{F}_5 = \frac{(t^2 + 250\,t + 3125)^3}{t^5},")
        st.markdown(
            "with $\\mathfrak{F}_5$ acting as $t \\mapsto 125/t$ — indeed "
            "$j_5(125/t)$ reproduces the second formula. Dialling $t$ walks along "
            "$X_0(5)$, and the two rational functions read off the domain and "
            "codomain $j$-invariants of the corresponding $5$-isogeny. (Which "
            "formula is *domain* is pinned down by the universal curve below: at "
            "$t = -11$ the first formula returns $-4096/11$, the $j$-invariant of "
            "the curve `11a3` $= X_1(11)$, which is the one *carrying* the "
            "$5$-torsion point.)"
        )
        st.markdown(
            "**Why these levels.** The repository stores such $j$-maps for the "
            "prime levels where $X_0(\\ell)$ is genus $0$ — $\\ell \\in \\{2, 3, "
            "5, 7, 13\\}$ — with $\\mathfrak{F}_\\ell$ always of the tidy form "
            "$t \\mapsto m/t$, where"
        )
        st.latex(r"m = \ell^{\,12/(\ell - 1)} \in \{4096,\ 729,\ 125,\ 49,\ 13\}.")
        st.markdown(
            "The exponent $12/(\\ell-1)$ is an integer exactly when $\\ell - 1 "
            "\\mid 12$ — that is, for $\\ell \\in \\{2,3,5,7,13\\}$, which is "
            "*precisely* the list of genus-$0$ prime levels. The same arithmetic "
            "that keeps $t \\mapsto m/t$ clean is what makes the curve a "
            "$\\mathbb{P}^1$."
        )
        st.markdown(
            "**Endomorphism points.** The diagonal condition $j_\\ell(x) = "
            "j_\\ell(m/x)$ — a curve $\\ell$-isogenous to itself, i.e. an "
            "$\\ell$-endomorphism, i.e. CM — has $2\\ell$ solutions, and every "
            "one of them lies on the circle $|x| = \\sqrt{m}$, the fixed locus "
            "of $\\mathfrak{F}_\\ell$ (where $m/x = \\bar{x}$). This is the "
            "algebraic mirror of the analytic fact from §8.2 that the self-dual "
            "locus is $|\\tau| = 1/\\sqrt{\\ell}$."
        )
        st.markdown(
            "Here is the algebraic model as a picture. Where the analytic $X_0("
            "\\ell)$ was a fundamental domain with edges to glue, the genus-$0$ "
            "model is just the **real $x$-line** $\\cong \\mathbb{P}^1("
            "\\mathbb{R})$, sitting in its ambient plane — no gluing, cusps at "
            "$0$ and $\\infty$. The plane (dimly coloured by $j_\\ell$, the §7 "
            "scheme) is only the backdrop. Overlaid is the **real locus of "
            "$j_\\ell$**, weighted Belyi-style: the preimage of the segment "
            "$[0, 1728]$ between the two finite critical values — the *dessin* "
            "— is bright, while the preimages of $j \\le 0$ and $j \\ge 1728$ "
            "stay faint. Red dots are the $2\\ell$ endomorphism points on the "
            "dashed self-dual circle."
        )
        components.html(modular_viz.genus0_plane_html(), height=625, scrolling=False)
        st.caption(
            "Levels $\\{3, 5, 7, 13\\}$ — the odd genus-$0$ primes tabulated in "
            "`jcoefs.json`, with $j_\\ell(x) = a_1(x)\\,a_3(x)^3\\,a_{-1}(x)^{-1}$ "
            "and $m = \\ell^{12/(\\ell-1)}$."
        )
        st.markdown(
            "**From the point back to the curves.** A point $x$ of $X_0(5)$ is "
            "supposed to *be* a pair $(E, C)$ — so we should be able to hand back "
            "actual equations. We can. The universal curve over $X_1(5)$ is the "
            "Tate normal form"
        )
        st.latex(r"E_t:\; y^2 + (1-t)\,xy - t\,y = x^3 - t\,x^2,")
        st.markdown(
            "whose point $(0,0)$ has exact order $5$; the diamond quotient down "
            "to $X_0(5)$ is $x = t - 11 - 1/t$, and solving back gives the "
            "section $t(x) = \\tfrac{1}{2}\\big(11 + x - \\sqrt{x^2 + 22x + "
            "125}\\big)$. Plugging $t(x)$ into $E_t$ recovers the **domain** "
            "curve with its kernel $C = \\langle(0,0)\\rangle$; the same recipe "
            "at the Fricke partner $m/x$ recovers the **codomain** $E/C$. The "
            "radicand $x^2 + 22x + 125$ has negative discriminant, so it is "
            "positive on *all* of $\\mathbb{R}$: every real $x$ yields two "
            "honest real elliptic curves and a real $5$-isogeny between them. "
            "Walk the real points of $X_0(5)$ and watch the pair move:"
        )
        components.html(modular_viz.x05_isogeny_html(), height=490, scrolling=False)
        st.caption(
            "Equation recovery (the $t(x)$ section into the Tate form) is the "
            "$\\ell = 5$ worked example, from the `x05pic.nb` computation. At "
            "$x = \\pm\\sqrt{125}$ the marker snaps to the self-dual point, "
            "where domain and codomain coincide."
        )
        st.markdown(
            "The full computations — the universal isogeny, larger $\\ell$, and "
            "the singularity resolution glossed over here — are carried out in "
            "the appendices of the author's papers on Mordell–Weil torsion and "
            "F-theory ([arXiv:1910.04095](https://arxiv.org/abs/1910.04095), "
            "2019) and on computing supersingular isogeny graphs via modular "
            "curves ([arXiv:2303.09096](https://arxiv.org/abs/2303.09096), 2023)."
        )

with st.expander("§9 — Modular polynomials: the Atkin and classical models", expanded=False):
    st.markdown(
        "§8 closed with the prettiest model of $X_0(\\ell)$: a Hauptmodul $t$ "
        "with both $j$-maps as explicit rational functions. But that model only "
        "exists at the five genus-$0$ levels. This section describes the two "
        "algebraic models that carry the actual computations everywhere else — "
        "the **Atkin polynomial**, built from the quotient of $X_0(\\ell)$ by "
        "the Fricke involution, and the **classical modular polynomial**, which "
        "presents $X_0(\\ell)$ through its two maps to the $j$-line. All three "
        "models answer the same question — *which $j$-invariants are "
        "$\\ell$-isogenous?* — from different coordinates."
    )
    tab_atk, tab_phi, tab_cmp = st.tabs([
        "The Atkin polynomial",
        "Classical modular polynomials",
        "Comparing the three models",
    ])

    # ── §9.1 — the Atkin polynomial ───────────────────────────────────────────
    with tab_atk:
        st.subheader("The Atkin Polynomial")
        st.markdown(
            "The Fricke involution $\\mathfrak{F}_\\ell$ (§8.2) swaps an isogeny "
            "with its dual. Taking the quotient by it — the **Atkin–Lehner "
            "quotient** — gives a new modular curve,"
        )
        st.latex(r"X_0(\ell)^{+} \;=\; X_0(\ell)\,/\,\langle \mathfrak{F}_\ell\rangle,")
        st.markdown(
            "whose points are cyclic $\\ell$-isogenies *with the direction "
            "forgotten*: the unordered pair $\\{E \\to E',\\ E' \\to E\\}$."
        )
        st.markdown(
            "**Functions on the quotient.** A function on $X_0(\\ell)$ descends "
            "to $X_0(\\ell)^{+}$ exactly when it is $\\mathfrak{F}_\\ell$-"
            "invariant, and the two $j$-maps symmetrize into a natural pair:"
        )
        st.latex(r"a_\ell = j_\ell + j_\ell\!\circ\!\mathfrak{F}_\ell, \qquad "
                 r"b_\ell = j_\ell \cdot \big(j_\ell\!\circ\!\mathfrak{F}_\ell\big).")
        st.markdown(
            "Nothing is lost: the two endpoints of the isogeny are recovered "
            "from $(a_\\ell, b_\\ell)$ as the roots of the **Atkin quadratic**"
        )
        st.latex(r"X^2 - a_\ell\,X + b_\ell \;=\; (X - j_\ell)\,"
                 r"\big(X - j_\ell\!\circ\!\mathfrak{F}_\ell\big).")
        st.markdown(
            "**Genus 0 upstairs vs. downstairs.** $X_0(\\ell)$ itself has genus "
            "$0$ for only five primes, but the *quotient* $X_0(\\ell)^{+}$ has "
            "genus $0$ for **fifteen**: $\\ell \\in \\{2, 3, 5, 7, 11, 13, 17, "
            "19, 23, 29, 31, 41, 47, 59, 71\\}$ — by Ogg's celebrated "
            "observation, precisely the primes dividing the order of the "
            "Monster. For these, pick a Hauptmodul $u$ on $X_0(\\ell)^{+}$; "
            "since $a_\\ell$ and $b_\\ell$ have poles only at the single cusp, "
            "they become **polynomials** in $u$, of degrees"
        )
        st.latex(r"\deg a_\ell = \ell, \qquad \deg b_\ell = \ell + 1, \qquad "
                 r"b_\ell = b_1 \cdot b_3^{\,3} \ \text{(a perfect-cube factor)}.")
        st.markdown(
            "The resulting **Atkin model** of $X_0(\\ell)$ is the curve"
        )
        st.latex(r"A_\ell(u, X) \;:\; X^2 - a_\ell(u)\,X + b_\ell(u) \;=\; 0")
        st.markdown(
            "in the $(u, X)$-plane: quadratic over the $u$-line (the double "
            "cover $X_0(\\ell) \\to X_0(\\ell)^{+}$, the sheet choosing which "
            "endpoint is the domain) and of degree $\\ell + 1$ over the "
            "$X = j$ line (the cover $X_0(\\ell) \\to X(1)$)."
        )
        st.markdown(
            "**Worked example: $\\ell = 5$.** The Fricke involution acts on the "
            "§8.3 Hauptmodul by $t \\mapsto 125/t$, so the invariant coordinate "
            "(in the normalisation of the repository's tables) is"
        )
        st.latex(r"u \;=\; t + \frac{125}{t} + 6,")
        st.markdown("and pushing the two $j$-maps through this substitution gives")
        st.latex(r"a_5(u) = u^5 - 670\,u^3 - 3800\,u^2 + 73056\,u + 449408,")
        st.latex(r"b_5(u) = \big(u^2 + 248\,u + 3856\big)^{3}.")
        st.markdown(
            "Six-digit integers — compare the $48$-digit coefficients of the "
            "classical $\\Phi_5$ in the next tab. This economy is the whole "
            "point of the Atkin model."
        )
        st.markdown(
            "**How it is used.** Over $\\mathbb{F}_p$, to find the $\\ell$-"
            "isogeny neighbours of a given $j$: the expression $j^2 - "
            "a_\\ell(u)\\,j + b_\\ell(u)$ is a degree-$(\\ell+1)$ polynomial in "
            "$u$, and its roots in $\\mathbb{F}_p$ are the points of "
            "$X_0(\\ell)^{+}$ whose pair contains $j$. For each root $u_0$, the "
            "partner is read off the quadratic: $j' = a_\\ell(u_0) - j$. The "
            "repository tabulates $(a_\\ell, b_1, b_3)$ for all fifteen Atkin "
            "primes in `atkinpolys.json`."
        )

    # ── §9.2 — classical modular polynomials ──────────────────────────────────
    with tab_phi:
        st.subheader("Classical Modular Polynomials")
        st.markdown(
            "The Atkin model borrowed a coordinate from the quotient curve. The "
            "classical model borrows *both* coordinates from $X(1)$: it "
            "describes the degree-$(\\ell+1)$ cover $X_0(\\ell) \\to X(1)$ "
            "purely in terms of $j$. On function fields, that cover is a field "
            "extension of $\\mathbb{C}(j)$ of degree $\\ell + 1$, generated by "
            "the second $j$-map, and the **classical modular polynomial** "
            "$\\Phi_\\ell$ is its minimal polynomial:"
        )
        st.latex(r"\mathbb{C}\big(X_0(\ell)\big) \;\cong\; "
                 r"\mathbb{C}(j)[X]\,/\,\big(\Phi_\ell(X, j)\big),")
        st.markdown(
            "monic of degree $\\ell + 1$ in $X$, with **integer** coefficients. "
            "Equivalently and more geometrically: the pair of $j$-maps embeds "
            "$X_0(\\ell)$ into $X(1) \\times X(1)$, and $\\Phi_\\ell(X, Y) = 0$ "
            "is the (highly singular) plane model cut out by that image. It "
            "plays the role of a *division polynomial for the $j$-line*: the "
            "roots of $\\Phi_\\ell(X, j_0)$ are precisely the $\\ell + 1$ "
            "$j$-invariants $\\ell$-isogenous to $j_0$."
        )
        st.markdown("Three classical properties, each one line:")
        st.markdown(
            "1. **Symmetry**: $\\Phi_\\ell(X, Y) = \\Phi_\\ell(Y, X)$ — the "
            "dual isogeny, again.\n"
            "2. **Analytic incarnation**: $\\Phi_\\ell\\big(j(\\ell\\tau), "
            "j(\\tau)\\big) = 0$ identically on $\\mathcal{H}$.\n"
            "3. **Kronecker congruence**: $\\Phi_\\ell(X, Y) \\equiv "
            "(X^{\\ell} - Y)(X - Y^{\\ell}) \\bmod \\ell$ — the shadow of "
            "Frobenius, and the reason these polynomials know about "
            "characteristic $p$.\n"
        )
        st.markdown(
            "The **diagonal** $\\Phi_\\ell(X, X) = 0$ picks out the $j$-"
            "invariants $\\ell$-isogenous to themselves — the endomorphism "
            "points drawn as red dots in §8.3, now as roots of an integer "
            "polynomial (it factors into Hilbert class polynomials, §7)."
        )
        st.markdown("**The example everyone shows** (because it is the only small one):")
        st.latex(r"""\Phi_2(X,Y) = X^3 + Y^3 - X^2Y^2 + 1488\,(X^2Y + XY^2)
                 - 162000\,(X^2 + Y^2)""")
        st.latex(r"""\qquad\qquad +\; 40773375\,XY + 8748000000\,(X + Y)
                 - 157464000000000.""")
        st.markdown(
            "Coefficient growth is ferocious — from the repository's own tables "
            "(`classical_modpolys.json`, $\\ell \\le 29$), the largest "
            "coefficient has:"
        )
        st.markdown(
            "| $\\ell$ | 2 | 3 | 5 | 11 |\n"
            "|---|---|---|---|---|\n"
            "| digits | 15 | 22 | 48 | 127 |\n"
        )
        st.markdown(
            "This is why nobody evaluates $\\Phi_\\ell$ over $\\mathbb{Q}$ if "
            "they can help it: the polynomials are reduced mod $p$ once and "
            "used over $\\mathbb{F}_p$, where the size disappears. Computing "
            "the tables in the first place is its own art (and an active side "
            "project of this repository)."
        )

    # ── §9.3 — comparing the three models ─────────────────────────────────────
    with tab_cmp:
        st.subheader("Three Models, One Curve")
        st.markdown(
            "All three constructions describe the *same* curve $X_0(\\ell)$; "
            "they differ in **which curve lends its coordinates**: $X_0(\\ell)$ "
            "itself (Hauptmodul), the Atkin–Lehner quotient plus one $j$-line "
            "(Atkin), or two $j$-lines (classical)."
        )
        st.markdown(
            "| | Hauptmodul model (§8.3) | Atkin model (§9.1) | classical $\\Phi_\\ell$ (§9.2) |\n"
            "|---|---|---|---|\n"
            "| lives in | $X_0(\\ell) \\cong \\mathbb{P}^1$ | $X_0(\\ell)^{+} \\times X(1)$ | $X(1) \\times X(1)$ |\n"
            "| available for | $\\ell \\in \\{2,3,5,7,13\\}$ | the 15 Monster primes | **every** $\\ell$ |\n"
            "| coefficient size | single digits | small (6 digits at $\\ell{=}5$) | enormous (127 digits at $\\ell{=}11$) |\n"
            "| a point gives | the **ordered** pair $(j, j')$ — and via $X_1(\\ell)$, the curves and kernel themselves | the **unordered** pair $\\{j, j'\\}$ | the neighbours of $j$, as roots |\n"
        )
        st.markdown(
            "Reading the columns left to right: each step **widens the range of "
            "$\\ell$ and pays for it in explicitness**. The Hauptmodul model "
            "hands you the isogeny itself (down to equations, §8.3); the Atkin "
            "model hands you the pair of endpoints but forgets the direction; "
            "the classical polynomial only answers *\"who are the "
            "neighbours?\"* — but answers it for every prime."
        )
        st.markdown(
            "The three lists are nested — $\\{2,3,5,7,13\\} \\subset "
            "\\{\\text{15 Atkin primes}\\} \\subset \\{\\text{all primes}\\}$ — "
            "and so is their use here. For isogeny *graphs* the forgotten "
            "direction costs nothing (the dual isogeny makes every edge "
            "two-way), so the Atkin model, with its tiny coefficients, is the "
            "workhorse at the fifteen primes where it exists; the classical "
            "$\\Phi_\\ell$ is the general-purpose fallback beyond them; and the "
            "genus-$0$ Hauptmodul models are reserved for when we need the "
            "isogeny itself, not just its endpoints."
        )


# ── Tab 0: Algebraic curves over ℝ ───────────────────────────────────────────
with tab0:
    st.subheader("Algebraic Curves over ℝ")
    st.markdown(
        "An **algebraic curve** over ℝ is the set of points $(x, y) \\in \\mathbb{R}^2$ "
        "satisfying an equation"
    )
    st.latex(r"f(x, y) = 0")
    st.markdown(
        "for some polynomial $f$ in two variables. "
        "This is a very broad family — most of the curves you encountered in school fit this definition."
    )

    st.markdown("#### Familiar examples")
    st.markdown(
        "Graphs of functions $y = g(x)$ are algebraic curves: just rewrite as $y - g(x) = 0$. "
        "But algebraic curves don't have to be function graphs — circles and ellipses aren't, "
        "yet they're still zero sets of polynomials."
    )

    # ── Small example plots ───────────────────────────────────────────────────
    N = 300
    u = np.linspace(-2.2, 2.2, N)
    X0, Y0 = np.meshgrid(u, u)

    examples = [
        ("Line",      r"$y = x$",          Y0 - X0),
        ("Parabola",  r"$y = x^2$",        Y0 - X0**2),
        ("Hyperbola", r"$xy = 1$",         X0 * Y0 - 1),
        ("Circle",    r"$x^2+y^2=1$",      X0**2 + Y0**2 - 1),
        ("Ellipse",   r"$x^2/4+y^2=1$",   X0**2 / 4 + Y0**2 - 1),
    ]

    fig_ex, axes = plt.subplots(1, 5, figsize=(11, 2.2))
    for ax_ex, (name, eq, Z_ex) in zip(axes, examples):
        ax_ex.contour(X0, Y0, Z_ex, levels=[0], colors=["steelblue"], linewidths=2)
        ax_ex.axhline(0, color="k", lw=0.4)
        ax_ex.axvline(0, color="k", lw=0.4)
        ax_ex.set_xlim(-2.2, 2.2)
        ax_ex.set_ylim(-2.2, 2.2)
        ax_ex.set_aspect("equal")
        ax_ex.set_frame_on(False)
        ax_ex.set_xticks([])
        ax_ex.set_yticks([])
        ax_ex.set_title(f"{name}\n{eq}", fontsize=9, linespacing=1.6)

    fig_ex.tight_layout(pad=0.5)
    st.pyplot(fig_ex)
    plt.close(fig_ex)

    st.markdown(
        "| Curve | Familiar form | As $f(x,y) = 0$ |\n"
        "|---|---|---|\n"
        "| Line | $y = x$ | $y - x = 0$ |\n"
        "| Parabola | $y = x^2$ | $y - x^2 = 0$ |\n"
        "| Hyperbola | $y = 1/x$ | $xy - 1 = 0$ |\n"
        "| Circle | $x^2 + y^2 = 1$ | $x^2 + y^2 - 1 = 0$ |\n"
        "| Ellipse | $\\frac{x^2}{4} + y^2 = 1$ |"
        " $\\frac{x^2}{4} + y^2 - 1 = 0$ |"
    )

    st.markdown("#### Quadratic forms and ellipses")
    st.markdown(
        "A **quadratic form** in two variables is a degree-2 polynomial with no linear or constant terms:"
    )
    st.latex(r"q(x, y) = ax^2 + bxy + cy^2.")
    st.markdown(
        "When $q$ is **positive definite** — meaning $q(x,y) > 0$ for all $(x,y) \\neq (0,0)$ — "
        "the level set $q(x,y) = 1$ is an ellipse. "
        "Positive definiteness is equivalent to $a > 0$ and $4ac - b^2 > 0$."
        "\n\n"
        "Use the controls below to explore how the shape of the ellipse depends on $(a, b, c)$. "
        "This family of curves is central to the rest of the project."
    )

    st.divider()

    ctrl0, plot0 = st.columns([1, 2])
    with ctrl0:
        st.markdown("**Quadratic form coefficients**")
        a_qf = st.number_input("a", value=2.0, step=0.1, key="bg0_a")
        b_qf = st.number_input("b", value=1.0, step=0.1, key="bg0_b")
        c_qf = st.number_input("c", value=2.0, step=0.1, key="bg0_c")

        disc = 4 * a_qf * c_qf - b_qf ** 2
        if a_qf <= 0 or disc <= 0:
            st.error(
                "The form is not positive definite. "
                "You need $a > 0$ and $4ac - b^2 > 0$."
            )
            qf_ok = False
        else:
            st.success(f"Positive definite  ✓  ($4ac - b^2 = {disc:.2f}$)")
            qf_ok = True

    with plot0:
        R = 2.5
        grid_n = 400
        xs = np.linspace(-R, R, grid_n)
        ys = np.linspace(-R, R, grid_n)
        X, Y = np.meshgrid(xs, ys)
        Z = a_qf * X**2 + b_qf * X * Y + c_qf * Y**2

        fig0, ax0 = plt.subplots(figsize=(5, 5))

        if qf_ok:
            # Background family of level sets for context
            for level, alpha in [(0.25, 0.15), (0.5, 0.2), (2.0, 0.2), (4.0, 0.15)]:
                ax0.contour(X, Y, Z, levels=[level],
                            colors=["steelblue"], alpha=alpha, linewidths=1)
            # Main level set q = 1
            ax0.contour(X, Y, Z, levels=[1.0],
                        colors=["steelblue"], linewidths=2.5)
            ax0.set_title(
                f"$q(x,y) = {a_qf}x^2 + {b_qf}xy + {c_qf}y^2 = 1$",
                fontsize=10,
            )
        else:
            ax0.text(0, 0, "Not positive definite",
                     ha="center", va="center", fontsize=12, color="gray")
            ax0.set_title("q(x, y) = 1")

        ax0.axhline(0, color="k", lw=0.5)
        ax0.axvline(0, color="k", lw=0.5)
        ax0.set_xlim(-R, R)
        ax0.set_ylim(-R, R)
        ax0.set_aspect("equal")
        ax0.set_frame_on(False)
        st.pyplot(fig0)
        plt.close(fig0)


# ── Tab 1: Real Elliptic Curves ───────────────────────────────────────────────
with tab1:
    st.subheader("Real Elliptic Curves")

    # ── Not ellipses ──────────────────────────────────────────────────────────
    st.markdown(
        "Despite the name, **elliptic curves are not ellipses**. "
        "The name is historical: elliptic curves arise when computing the arc length "
        "of an ellipse, but the curves themselves are a different object entirely. "
        "What they share with ellipses is that they are a *family* of curves, each "
        "defined by an equation of a specific form. For elliptic curves, that form is"
    )
    st.latex(r"y^2 = x^3 + fx + g, \qquad f, g \in \mathbb{R}.")
    st.markdown(
        "Different choices of $f$ and $g$ give different curves — "
        "but they all belong to the same family."
    )

    # ── Singular curves ───────────────────────────────────────────────────────
    st.markdown("#### Singular curves")
    st.markdown(
        "Not every choice of $f$ and $g$ gives a well-behaved curve. "
        "When the cubic $x^3 + fx + g$ has a repeated root, the curve develops "
        "a **singularity** — a point where it crosses or pinches itself. "
        "Singular curves are excluded from the definition of an elliptic curve. "
        "There are two types:\n\n"
        "- **Node** ($f = -3,\\, g = 2$): the cubic has a double root at $x = 1$, "
        "and the curve self-intersects there.\n"
        "- **Cusp** ($f = 0,\\, g = 0$): the cubic has a triple root at $x = 0$, "
        "and the curve has a sharp pinch point at the origin."
    )

    def _plot_ec(ax, f, g, xlim, ylim, title, singular_pt=None):
        """Plot y^2 = x^3 + fx + g on ax."""
        xs = np.linspace(xlim[0], xlim[1], 2000)
        y2 = xs**3 + f*xs + g
        yp = np.where(y2 >= 0, np.sqrt(np.clip(y2, 0, None)), np.nan)
        ax.plot(xs,  yp, color="steelblue", lw=2)
        ax.plot(xs, -yp, color="steelblue", lw=2)
        ax.axhline(0, color="k", lw=0.4)
        ax.axvline(0, color="k", lw=0.4)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal")
        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title, fontsize=10)
        if singular_pt is not None:
            ax.scatter(*singular_pt, color="red", s=60, zorder=5)

    fig_sing, (ax_nd, ax_cu) = plt.subplots(1, 2, figsize=(8, 3.5))
    _plot_ec(ax_nd, -3, 2, (-2.3, 2.8), (-2.5, 2.5),
             "Node:  $y^2 = x^3 - 3x + 2$", singular_pt=(1, 0))
    _plot_ec(ax_cu,  0, 0, (-0.5, 2.5), (-2.5, 2.5),
             "Cusp:  $y^2 = x^3$",          singular_pt=(0, 0))
    fig_sing.tight_layout()
    st.pyplot(fig_sing)
    plt.close(fig_sing)

    # ── Two smooth flavors ────────────────────────────────────────────────────
    st.markdown("#### Two visual flavors")
    st.markdown(
        "When $f$ and $g$ are chosen so that the curve is smooth, it always comes "
        "in one of two shapes, determined by whether the cubic $x^3 + fx + g$ has "
        "one or three real roots:\n\n"
        "- **One component** (one real root, $\\Delta < 0$): a single unbounded branch.\n"
        "- **Two components** (three real roots, $\\Delta > 0$): an unbounded branch "
        "plus a compact oval.\n\n"
        "Here $\\Delta = -16(4f^3 + 27g^2)$ is the discriminant."
    )

    fig_smth, (ax_one, ax_two) = plt.subplots(1, 2, figsize=(8, 3.5))
    _plot_ec(ax_one, 0,  1, (-1.3, 2.5), (-2.5, 2.5),
             "One component:  $y^2 = x^3 + 1$")
    _plot_ec(ax_two, -1, 0, (-1.3, 2.5), (-2.5, 2.5),
             "Two components:  $y^2 = x^3 - x$")
    fig_smth.tight_layout()
    st.pyplot(fig_smth)
    plt.close(fig_smth)

    # ── Curve explorer ────────────────────────────────────────────────────────
    st.markdown("#### Explore")
    st.markdown(
        "The two shapes live in one family. **Drag the point in the "
        "$(f, g)$-plane**: the gold curve is the discriminant locus "
        "$4f^3 + 27g^2 = 0$, and crossing it is exactly where the oval pinches "
        "off. The marker snaps onto the locus nearby, so you can see the "
        "singular curves themselves — a node anywhere on it, and the cusp at "
        "$f = g = 0$."
    )
    components.html(basics_viz.cubic_family_html(), height=470, scrolling=False)

    st.divider()

    # ── Group law ─────────────────────────────────────────────────────────────
    st.markdown("#### The group law")
    st.markdown(
        "What makes elliptic curves remarkable is that the points on a smooth curve "
        "— together with one extra 'hidden' point called the **point at infinity**, "
        "written $\\mathcal{O}$ — form an **abelian group**.\n\n"
        "$\\mathcal{O}$ is the identity element: $P + \\mathcal{O} = P$ for every point $P$. "
        "Geometrically it lives 'at the end of every vertical line', infinitely far in "
        "the $y$-direction.\n\n"
        "The group law has a clean geometric description: "
        "**three points sum to $\\mathcal{O}$ if and only if they are collinear.** "
        "To compute $P + Q$:"
    )
    st.markdown(
        "1. Draw the line through $P$ and $Q$ "
        "(or the tangent to the curve at $P$ if $P = Q$).\n"
        "2. This line meets the curve in a third point $R$.\n"
        "3. Reflect $R$ over the $x$-axis to obtain $P + Q$.\n\n"
        "*(Reflection gives the inverse: since $P + Q + R = \\mathcal{O}$, "
        "we get $P + Q = -R$, and the inverse of $(x, y)$ is $(x, -y)$.)*"
    )
    st.markdown("Use the applet below to see this in action.")

    st.divider()

    # ── Chord-tangent applet: click two points, compute their sum ─────────────
    _XR = 3.3   # visible x half-range

    _YR = 4.5   # visible y half-range
    _O_XY1 = (0.0, 4.2)   # 𝒪 marker, top centre

    def _real_line(P, Q, f):
        """Two endpoints of the chord/tangent through P, Q (vertical if 𝒪/inverse)."""
        if P == _OINF and Q == _OINF:
            return None
        if P == _OINF or Q == _OINF:
            fin = Q if P == _OINF else P
            return [(fin[0], -_YR), (fin[0], _YR)]
        x1, y1 = P
        x2, y2 = Q
        double = abs(x1 - x2) < 1e-9 and abs(y1 - y2) < 1e-9
        if (double and abs(y1) < 1e-9) or (not double and abs(x1 - x2) < 1e-9):
            return [(x1, -_YR), (x1, _YR)]
        m = (3 * x1 * x1 + f) / (2 * y1) if double else (y2 - y1) / (x2 - x1)
        return [(-_XR, m * (-_XR - x1) + y1), (_XR, m * (_XR - x1) + y1)]

    def _real_sum(P, Q, f):
        """(R, sum, message). sum is _OINF (𝒪), a point, or None; R None at 𝒪."""
        if P == _OINF and Q == _OINF:
            return None, _OINF, "𝒪 + 𝒪 = 𝒪."
        if P == _OINF:
            return None, Q, "𝒪 + Q = Q (𝒪 is the identity)."
        if Q == _OINF:
            return None, P, "P + 𝒪 = P (𝒪 is the identity)."
        x1, y1 = P
        x2, y2 = Q
        double = abs(x1 - x2) < 1e-9 and abs(y1 - y2) < 1e-9
        if double and abs(y1) < 1e-9:
            return None, _OINF, "P is a 2-torsion point, so 2P = 𝒪."
        if (not double) and abs(x1 - x2) < 1e-9:
            return None, _OINF, "P and Q are inverses, so P + Q = 𝒪."
        m  = (3 * x1 * x1 + f) / (2 * y1) if double else (y2 - y1) / (x2 - x1)
        x3 = m * m - x1 - x2
        yR = m * (x3 - x1) + y1
        return (x3, yR), (x3, -yR), None

    def _real_fmt(pt):
        return "𝒪" if pt == _OINF else f"({pt[0]:.2f}, {pt[1]:.2f})"

    ctrl1, plot1 = st.columns([1, 2])

    with ctrl1:
        st.markdown("**Curve**  $y^2 = x^3 + fx + g$")
        f_val = st.slider("f", -5.0, 5.0, -1.0, 0.1, key="bg1_f")
        g_val = st.slider("g", -5.0, 5.0,  1.0, 0.1, key="bg1_g")
        disc  = -16 * (4 * f_val**3 + 27 * g_val**2)
        _ok1  = abs(disc) >= 1e-6
        if not _ok1:
            st.warning("Singular curve (Δ ≈ 0) — adjust f or g.")

    # Clickable sample points on the curve (both branches).
    _samples = []
    if _ok1:
        for x in np.round(np.linspace(-_XR, _XR, 140), 4):
            y2 = x**3 + f_val * x + g_val
            if y2 < -1e-9:
                continue
            y = float(np.sqrt(max(y2, 0.0)))
            if y < 1e-6:
                _samples.append((float(x), 0.0))
            else:
                _samples.append((float(x), round(y, 4)))
                _samples.append((float(x), round(-y, 4)))

    _samp_set = set(_samples)
    _o_idx1   = len(_samples)              # customdata index of 𝒪
    sel1 = [pt for pt in st.session_state.get("bg1_sel", [])
            if pt == _OINF or pt in _samp_set]
    st.session_state["bg1_sel"] = sel1
    mode1 = st.session_state.get("bg1_mode") if len(sel1) == 2 else None

    R1 = S1 = msg1 = None
    line1 = None
    if mode1 in ("line", "sum"):
        line1 = _real_line(sel1[0], sel1[1], f_val)
    if mode1 == "sum":
        R1, S1, msg1 = _real_sum(sel1[0], sel1[1], f_val)
    is_double1 = (len(sel1) == 2 and _OINF not in sel1
                  and abs(sel1[0][0]-sel1[1][0]) < 1e-9 and abs(sel1[0][1]-sel1[1][1]) < 1e-9)
    sum_lbl = "2P" if is_double1 else "P + Q"

    with ctrl1:
        st.markdown("**Pick two points**")
        st.caption("Click the curve — or **𝒪** — to choose **P** then **Q** "
                   "(same spot twice = **2P**; a third click clears).")
        for nm, pt in zip(("P", "Q"), sel1):
            st.markdown(f"- **{nm}** = {_real_fmt(pt)}")
        c1a, c1b, c1c = st.columns(3)
        if c1a.button("Clear", key="bg1_clear", width="stretch"):
            st.session_state["bg1_sel"] = []
            st.session_state["bg1_mode"] = None
            st.rerun()
        if c1b.button("Show line", key="bg1_line", width="stretch",
                      disabled=(len(sel1) != 2)):
            st.session_state["bg1_mode"] = "line"
            st.rerun()
        if c1c.button("Compute sum", key="bg1_compute",
                      width="stretch", disabled=(len(sel1) != 2)):
            st.session_state["bg1_mode"] = "sum"
            st.rerun()
        if mode1 == "sum":
            if msg1:
                st.info(msg1)
            else:
                st.success(f"{sum_lbl} = {_real_fmt(S1)}")

    with plot1:
        if not _ok1:
            st.info("Adjust f and g to get a smooth curve.")
        else:
            xs  = np.linspace(-_XR, _XR, 1200)
            ys2 = xs**3 + f_val * xs + g_val
            yp  = np.where(ys2 >= 0, np.sqrt(np.clip(ys2, 0, None)), np.nan)

            _role1 = {pt: nm for nm, pt in zip(("P", "Q"), sel1)}
            sum_is_inf1 = (mode1 == "sum" and S1 == _OINF)
            scols, ssizes = [], []
            for pt in _samples:
                if pt in _role1:
                    scols.append("red" if _role1[pt] == "P" else "green")
                    ssizes.append(15)
                else:
                    scols.append("steelblue"); ssizes.append(7)

            fig = go.Figure()
            # curve (visual, both branches)
            fig.add_trace(go.Scatter(x=xs, y=yp, mode="lines",
                          line=dict(color="steelblue", width=2),
                          hoverinfo="skip", showlegend=False))
            fig.add_trace(go.Scatter(x=xs, y=-yp, mode="lines",
                          line=dict(color="steelblue", width=2),
                          hoverinfo="skip", showlegend=False))
            # line through P, Q
            if line1 is not None:
                fig.add_trace(go.Scatter(
                    x=[line1[0][0], line1[1][0]], y=[line1[0][1], line1[1][1]],
                    mode="lines", line=dict(color="orange", width=1.5, dash="dash"),
                    hoverinfo="skip", showlegend=False))
            # full sum construction (R + reflected sum)
            if mode1 == "sum" and R1 is not None and S1 not in (None, _OINF):
                fig.add_trace(go.Scatter(x=[R1[0], S1[0]], y=[R1[1], S1[1]],
                              mode="lines", line=dict(color="gray", width=1, dash="dot"),
                              hoverinfo="skip", showlegend=False))
                fig.add_trace(go.Scatter(
                    x=[R1[0], S1[0]], y=[R1[1], S1[1]], mode="markers+text",
                    marker=dict(color=["mediumpurple", "orange"], size=[13, 16]),
                    text=["R", sum_lbl], textposition="top center",
                    textfont=dict(size=12), hoverinfo="skip", showlegend=False))
            # clickable curve samples
            fig.add_trace(go.Scatter(
                x=[s[0] for s in _samples], y=[s[1] for s in _samples],
                mode="markers", marker=dict(color=scols, size=ssizes),
                customdata=list(range(len(_samples))),
                hovertext=[f"({s[0]:.2f}, {s[1]:.2f})" for s in _samples],
                hoverinfo="text", showlegend=False,
                selected=dict(marker=dict(opacity=1.0)),
                unselected=dict(marker=dict(opacity=1.0))))
            # the point at infinity 𝒪 (clickable)
            if _OINF in _role1:
                o_color1, o_size1 = ("red" if _role1[_OINF] == "P" else "green"), 15
            elif sum_is_inf1:
                o_color1, o_size1 = "orange", 17
            else:
                o_color1, o_size1 = "olive", 12
            fig.add_trace(go.Scatter(
                x=[_O_XY1[0]], y=[_O_XY1[1]], mode="markers+text",
                marker=dict(color=o_color1, size=o_size1,
                            line=dict(color="white", width=0.5)),
                text=["𝒪"], textposition="top center",
                textfont=dict(size=13, color="olive"),
                customdata=[_o_idx1], hovertext=["point at infinity 𝒪"],
                hoverinfo="text", showlegend=False,
                selected=dict(marker=dict(opacity=1.0)),
                unselected=dict(marker=dict(opacity=1.0))))
            fig.update_layout(
                margin=dict(l=10, r=10, t=38, b=10), height=470,
                plot_bgcolor="white", showlegend=False,
                title=dict(text=f"y² = x³ + {f_val:.1f}x + {g_val:.1f}",
                           x=0.5, xanchor="center", font=dict(size=14)))
            fig.update_xaxes(range=[-_XR, _XR], zeroline=True, zerolinecolor="lightgray",
                             showgrid=False)
            fig.update_yaxes(range=[-_YR, _YR], zeroline=True, zerolinecolor="lightgray",
                             showgrid=False)
            ev = st.plotly_chart(fig, width="stretch", on_select="rerun",
                                 selection_mode="points", key="bg1_chart")

            idx = _click_index(ev)
            if idx is not None:
                clicked = _OINF if idx == _o_idx1 else (_samples[idx] if 0 <= idx < len(_samples) else None)
                if clicked is not None:
                    _two_point_pick("bg1", clicked)
                    st.rerun()


# ── Tab 2: Elliptic curves over 𝔽ₚ ──────────────────────────────────────────
with tab2:
    st.subheader("Elliptic Curves over $\\mathbb{F}_p$")

    # ── Explanation ───────────────────────────────────────────────────────────
    st.markdown(
        "Everything from the previous two lessons carries over to the finite field "
        "$\\mathbb{F}_p = \\{0, 1, \\ldots, p-1\\}$, with one change: "
        "instead of the real numbers, all coordinates and all arithmetic are taken **mod $p$**."
    )
    st.markdown(
        "The **ambient space** is no longer the continuous real plane — it becomes a "
        "discrete $p \\times p$ grid of points $\\mathbb{F}_p^2$. "
        "An algebraic curve over $\\mathbb{F}_p$ is the zero set"
    )
    st.latex(r"\{(x, y) \in \mathbb{F}_p^2 : f(x, y) \equiv 0 \pmod{p}\}.")
    st.markdown(
        "Every curve we studied over $\\mathbb{R}$ has an exact analogue here. "
        "We represent elements of $\\mathbb{F}_p$ using integers in "
        "$\\{-\\lfloor p/2 \\rfloor, \\ldots, \\lfloor p/2 \\rfloor\\}$ "
        "so the grid is centred at the origin."
    )

    # ── Familiar curves over F17 ──────────────────────────────────────────────
    st.markdown("#### Familiar curves over $\\mathbb{F}_{17}$")
    st.markdown(
        "Here are the same five curves from the first lesson, now over $\\mathbb{F}_{17}$. "
        "The full $17 \\times 17$ ambient grid is shown in gray."
    )

    def _fp_sym(x, p):
        xr = x % p
        return xr if 2 * xr < p else xr - p

    def _fp_curve_pts(poly_fn, p):
        h = p // 2
        return [(x, y)
                for x in range(-h, h + 1)
                for y in range(-h, h + 1)
                if poly_fn(x, y) % p == 0]

    _P17   = 17
    _AMB17 = [(_fp_sym(x, _P17), _fp_sym(y, _P17))
              for x in range(_P17) for y in range(_P17)]

    _EX17 = [
        ("Line",      r"$y = x$",          _fp_curve_pts(lambda x, y: y - x,           _P17)),
        ("Parabola",  r"$y = x^2$",        _fp_curve_pts(lambda x, y: y - x*x,         _P17)),
        ("Hyperbola", r"$xy = 1$",         _fp_curve_pts(lambda x, y: x*y - 1,         _P17)),
        ("Circle",    r"$x^2+y^2=1$",      _fp_curve_pts(lambda x, y: x*x + y*y - 1,  _P17)),
        ("Ellipse",   r"$x^2+4y^2=4$",     _fp_curve_pts(lambda x, y: x*x+4*y*y-4,    _P17)),
    ]

    fig17, axes17 = plt.subplots(1, 5, figsize=(11, 2.5))
    _ax17 = [a for a in axes17]
    _ax17 = list(axes17)
    for ax_i, (name, eq, pts) in zip(_ax17, _EX17):
        ax_i.scatter([p[0] for p in _AMB17], [p[1] for p in _AMB17],
                     color="gray", alpha=0.25, s=5, zorder=1)
        ax_i.scatter([p[0] for p in pts], [p[1] for p in pts],
                     color="steelblue", s=12, zorder=3)
        ax_i.set_xlim(-9, 9); ax_i.set_ylim(-9, 9)
        ax_i.set_aspect("equal"); ax_i.set_frame_on(False)
        ax_i.set_xticks([]); ax_i.set_yticks([])
        ax_i.set_title(f"{name}\n{eq}", fontsize=9, linespacing=1.6)
    fig17.tight_layout(pad=0.5)
    st.pyplot(fig17)
    plt.close(fig17)

    # ── Group law ─────────────────────────────────────────────────────────────
    st.markdown("#### Elliptic curves and the group law")
    st.markdown(
        "An elliptic curve over $\\mathbb{F}_p$ is the set of solutions to"
    )
    st.latex(r"y^2 \equiv x^3 + fx + g \pmod{p}")
    st.markdown(
        "together with a point at infinity $\\mathcal{O}$, provided the curve is smooth "
        "($\\Delta = -16(4f^3 + 27g^2) \\not\\equiv 0 \\pmod{p}$). "
        "The **group law is identical** to the real case: three points sum to "
        "$\\mathcal{O}$ if and only if they are collinear, where lines are now "
        "$\\mathbb{F}_p$-lines (solutions to $ax + by \\equiv c \\pmod{p}$). "
        "Use the applet below to explore."
    )

    st.divider()

    # ── Fp helpers ────────────────────────────────────────────────────────────
    def _ec_pts(f, g, p):
        """All affine points on y^2 = x^3+fx+g over Fp (symmetric coords)."""
        h   = p // 2
        pts = []
        for x in range(-h, h + 1):
            y2 = (pow(x, 3, p) + f * x + g) % p
            for y in range(p):
                if pow(y, 2, p) == y2:
                    pts.append((x, _fp_sym(y, p)))
        return pts

    def _chord_pts(x1, y1, x2, y2, p):
        """All Fp points on the line through (x1,y1) and (x2,y2)."""
        if (x2 - x1) % p == 0:                        # vertical line
            return [(_fp_sym(x1, p), _fp_sym(y, p)) for y in range(p)]
        m = (y2 - y1) * pow((x2 - x1) % p, -1, p) % p
        b = (y1 - m * x1) % p
        return [(_fp_sym(x, p), _fp_sym((m * x + b) % p, p)) for x in range(p)]

    def _tangent_pts(x0, y0, f, p):
        """All Fp points on the tangent to y^2=x^3+fx+g at (x0,y0), y0≠0."""
        m = (3 * pow(x0, 2, p) + f) * pow((2 * y0) % p, -1, p) % p
        b = (y0 - m * x0) % p
        return [(_fp_sym(x, p), _fp_sym((m * x + b) % p, p)) for x in range(p)]

    # ── Applet: click two points, compute their sum ───────────────────────────
    _PRIMES = [n for n in range(5, 72)
               if n > 1 and all(n % d != 0 for d in range(2, n))]

    def _fp_line(P, Q, f, p):
        """Fp points on the line through P, Q (tangent if P==Q; vertical if 𝒪)."""
        if P == _OINF and Q == _OINF:
            return []
        if P == _OINF or Q == _OINF:
            fin = Q if P == _OINF else P
            return [(_fp_sym(fin[0], p), _fp_sym(y, p)) for y in range(p)]
        x1, y1 = P
        x2, y2 = Q
        same = (x1 - x2) % p == 0 and (y1 - y2) % p == 0
        if same:
            if y1 % p == 0:                       # 2-torsion: tangent is vertical
                return [(_fp_sym(x1, p), _fp_sym(y, p)) for y in range(p)]
            return _tangent_pts(x1, y1, f, p)
        if (x1 - x2) % p == 0:                    # vertical chord (inverses)
            return [(_fp_sym(x1, p), _fp_sym(y, p)) for y in range(p)]
        return _chord_pts(x1, y1, x2, y2, p)

    def _fp_group_sum(P, Q, f, p):
        """(sum, R, message). sum is _OINF (𝒪), a point, or None; R None at 𝒪."""
        if P == _OINF and Q == _OINF:
            return _OINF, None, "𝒪 + 𝒪 = 𝒪."
        if P == _OINF:
            return Q, None, "𝒪 + Q = Q (𝒪 is the identity)."
        if Q == _OINF:
            return P, None, "P + 𝒪 = P (𝒪 is the identity)."
        x1, y1 = P
        x2, y2 = Q
        same = (x1 - x2) % p == 0 and (y1 - y2) % p == 0
        if same and y1 % p == 0:
            return _OINF, None, "P = Q is a 2-torsion point, so 2P = 𝒪."
        if (x1 - x2) % p == 0 and (y1 + y2) % p == 0 and not same:
            return _OINF, None, "P and Q are inverses, so P + Q = 𝒪."
        if same:
            m = (3 * pow(x1, 2, p) + f) * pow((2 * y1) % p, -1, p) % p
        else:
            m = (y2 - y1) * pow((x2 - x1) % p, -1, p) % p
        x3 = (m * m - x1 - x2) % p
        yR = (m * (x3 - x1) + y1) % p
        R  = (_fp_sym(x3, p), _fp_sym(yR, p))
        S  = (_fp_sym(x3, p), _fp_sym((-yR) % p, p))
        return S, R, None

    def _fp_fmt(pt):
        return "𝒪" if pt == _OINF else f"({pt[0]}, {pt[1]})"

    app_left, app_right = st.columns([1, 2])

    with app_left:
        st.markdown("**Field**")
        _p = st.selectbox("p", _PRIMES, index=_PRIMES.index(17), key="bg2_p")
        st.markdown("**Curve**  $y^2 = x^3 + fx + g$")
        _f = int(st.number_input("f", value=0, step=1, key="bg2_f"))
        _g = int(st.number_input("g", value=1, step=1, key="bg2_g"))
        _disc     = (-16 * (4 * pow(_f, 3) + 27 * pow(_g, 2))) % _p
        _curve_ok = (_disc != 0)
        if not _curve_ok:
            st.warning("Singular curve mod p — adjust f or g.")
        else:
            st.success("Smooth curve ✓")

    if not _curve_ok:
        with app_right:
            st.info("Adjust f and g to get a smooth curve.")
    else:
        _pts    = _ec_pts(_f, _g, _p)
        _pt_set = set(_pts)
        _h      = _p // 2
        _o_xy   = (_h + 1.4, _h + 1.4)            # 𝒪 marker, top-right corner
        _o_idx  = len(_pts)                       # customdata index of 𝒪
        # selection persists; finite points must stay on the curve, 𝒪 is always valid
        sel = [pt for pt in st.session_state.get("bg2_sel", [])
               if pt == _OINF or pt in _pt_set]
        st.session_state["bg2_sel"] = sel
        mode = st.session_state.get("bg2_mode") if len(sel) == 2 else None

        S = R = msg = None
        line = []
        if mode in ("line", "sum"):
            line = _fp_line(sel[0], sel[1], _f, _p)
        if mode == "sum":
            S, R, msg = _fp_group_sum(sel[0], sel[1], _f, _p)

        with app_left:
            st.markdown("**Pick two points**")
            st.caption("Click curve points — or **𝒪** — to choose **P** then **Q** "
                       "(same point twice = **2P**; a third click clears).")
            for nm, pt in zip(("P", "Q"), sel):
                st.markdown(f"- **{nm}** = {_fp_fmt(pt)}")
            c1, c2, c3 = st.columns(3)
            if c1.button("Clear", key="bg2_clear", width="stretch"):
                st.session_state["bg2_sel"] = []
                st.session_state["bg2_mode"] = None
                st.rerun()
            if c2.button("Show line", key="bg2_line", width="stretch",
                         disabled=(len(sel) != 2)):
                st.session_state["bg2_mode"] = "line"
                st.rerun()
            if c3.button("Compute sum", key="bg2_compute", width="stretch",
                         disabled=(len(sel) != 2)):
                st.session_state["bg2_mode"] = "sum"
                st.rerun()
            if mode == "sum":
                if msg:
                    st.info(msg)
                else:
                    st.success(f"P + Q = {_fp_fmt(S)}")

        with app_right:
            _amb = [(_fp_sym(x, _p), _fp_sym(y, _p))
                    for x in range(_p) for y in range(_p)]
            _role = {pt: nm for nm, pt in zip(("P", "Q"), sel)}
            sum_is_inf = (mode == "sum" and S == _OINF)

            colors, sizes = [], []
            for pt in _pts:
                if mode == "sum" and S not in (None, _OINF) and pt == S:
                    colors.append("orange"); sizes.append(17)
                elif mode == "sum" and R is not None and pt == R:
                    colors.append("mediumpurple"); sizes.append(14)
                elif pt in _role:
                    colors.append("red" if _role[pt] == "P" else "green")
                    sizes.append(16)
                else:
                    colors.append("steelblue"); sizes.append(9)

            if _OINF in _role:
                o_color, o_size = ("red" if _role[_OINF] == "P" else "green"), 16
            elif sum_is_inf:
                o_color, o_size = "orange", 18
            else:
                o_color, o_size = "olive", 13

            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=[q[0] for q in _amb], y=[q[1] for q in _amb],
                mode="markers", marker=dict(color="lightgray", size=5),
                hoverinfo="skip", showlegend=False))
            if line:
                fig.add_trace(go.Scatter(
                    x=[q[0] for q in line], y=[q[1] for q in line],
                    mode="markers", marker=dict(color="darkgray", size=7),
                    hoverinfo="skip", showlegend=False))
            fig.add_trace(go.Scatter(
                x=[q[0] for q in _pts], y=[q[1] for q in _pts],
                mode="markers",
                marker=dict(color=colors, size=sizes,
                            line=dict(color="white", width=0.5)),
                customdata=list(range(len(_pts))),
                hovertext=[f"({q[0]}, {q[1]})" for q in _pts], hoverinfo="text",
                showlegend=False,
                selected=dict(marker=dict(opacity=1.0)),
                unselected=dict(marker=dict(opacity=1.0))))
            fig.add_trace(go.Scatter(     # the point at infinity 𝒪 (clickable)
                x=[_o_xy[0]], y=[_o_xy[1]], mode="markers+text",
                marker=dict(color=o_color, size=o_size,
                            line=dict(color="white", width=0.5)),
                text=["𝒪"], textposition="top center",
                textfont=dict(size=13, color="olive"),
                customdata=[_o_idx], hovertext=["point at infinity 𝒪"],
                hoverinfo="text", showlegend=False,
                selected=dict(marker=dict(opacity=1.0)),
                unselected=dict(marker=dict(opacity=1.0))))
            fig.update_layout(
                margin=dict(l=10, r=10, t=38, b=10), height=470,
                plot_bgcolor="white", showlegend=False,
                title=dict(text=f"y² = x³ + {_f % _p}x + {_g % _p}  (mod {_p})",
                           x=0.5, xanchor="center", font=dict(size=14)))
            fig.update_xaxes(visible=False, range=[-_h - 0.6, _h + 2.0])
            fig.update_yaxes(visible=False, scaleanchor="x", scaleratio=1,
                             range=[-_h - 0.6, _h + 2.0])
            ev = st.plotly_chart(fig, width="stretch", on_select="rerun",
                                 selection_mode="points", key="bg2_chart")

        idx = _click_index(ev)
        if idx is not None:
            clicked = _OINF if idx == _o_idx else (_pts[idx] if 0 <= idx < len(_pts) else None)
            if clicked is not None:
                _two_point_pick("bg2", clicked)
                st.rerun()


# ── Tab 3: Elliptic curves over ℂ ────────────────────────────────────────────
with tab3:
    st.subheader("Elliptic Curves over ℂ")

    # ── Same definition, new problem ──────────────────────────────────────────
    st.markdown(
        "An elliptic curve over $\\mathbb{C}$ is defined by exactly the same equation"
    )
    st.latex(r"y^2 = x^3 + fx + g, \qquad f, g \in \mathbb{C},")
    st.markdown(
        "with the same smoothness condition $\\Delta \\neq 0$. "
        "But there is an immediate obstacle to making pictures: the solution set lives in "
        "$\\mathbb{C}^2 \\cong \\mathbb{R}^4$, which requires four real dimensions to draw. "
        "We need a different approach."
    )

    # ── The bijection ─────────────────────────────────────────────────────────
    st.markdown("#### A different perspective: complex tori")
    st.markdown(
        "Fortunately, there is a completely different way to think about elliptic curves "
        "over $\\mathbb{C}$. Every elliptic curve over $\\mathbb{C}$ is isomorphic, "
        "as a complex manifold, to a **complex torus**"
    )
    st.latex(r"\mathbb{C}/\Lambda, \qquad \Lambda = \mathbb{Z} + \tau\,\mathbb{Z},")
    st.markdown(
        "for some $\\tau$ in the upper half-plane (Im $\\tau > 0$). "
        "The torus $\\mathbb{C}/\\Lambda$ is the complex plane with points identified "
        "whenever they differ by an element of $\\Lambda$: every point has an equivalent "
        "representative in the fundamental parallelogram with vertices $0, 1, 1+\\tau, \\tau$.\n\n"
        "The bijection goes in both directions via the **Weierstrass $\\wp$-function**: "
        "given a lattice $\\Lambda$, the map"
    )
    st.latex(r"z \;\longmapsto\; \bigl(\wp(z),\; \wp'(z)\bigr)")
    st.markdown(
        "sends $\\mathbb{C}/\\Lambda$ isomorphically onto an algebraic elliptic curve, "
        "with specific coefficients $f$ and $g$ determined by the lattice. "
        "Conversely, every algebraic elliptic curve arises this way. "
        "The isomorphism class of the torus is determined by $\\tau$ up to the action of "
        "$\\mathrm{SL}(2,\\mathbb{Z})$ by Möbius transformations, so the space of "
        "isomorphism classes of elliptic curves over $\\mathbb{C}$ is parametrised by"
    )
    st.latex(r"\mathfrak{h} / \mathrm{SL}(2,\mathbb{Z}),")
    st.markdown("where $\\mathfrak{h}$ is the upper half-plane.")

    # ── Why we can't compute exactly ─────────────────────────────────────────
    st.markdown("#### Why we can't compute the bijection exactly")
    st.markdown(
        "The functions involved — the Weierstrass $\\wp$-function, the Eisenstein series "
        "$g_2$ and $g_3$ — are defined by infinite series and are **transcendental**: "
        "they cannot be expressed in terms of elementary functions. "
        "Going from the algebraic to the analytic picture requires computing integrals of the form"
    )
    st.latex(r"\int \frac{dx}{\sqrt{x^3 + fx + g}},")
    st.markdown(
        "known as **elliptic integrals**. These are famously impossible to evaluate "
        "in closed form using elementary functions — and this is precisely where elliptic "
        "curves get their name. The arc length of an ellipse leads to exactly such an "
        "integral, and the impossibility of evaluating it 'nicely' was one of the driving "
        "problems of 19th-century mathematics."
    )

    # ── Why analytic models are better ───────────────────────────────────────
    st.markdown("#### Why the analytic picture is so useful")
    st.markdown(
        "Even though we cannot go back and forth between the two pictures by hand, "
        "the complex torus $\\mathbb{C}/\\Lambda$ offers two major advantages over "
        "the algebraic model:\n\n"
        "**1. We can make pictures.** The torus is a genuine 2-dimensional surface, "
        "and its fundamental domain is just a parallelogram — something we can draw.\n\n"
        "**2. The group law becomes trivial.** On $\\mathbb{C}/\\Lambda$, the group "
        "law is simply complex addition mod $\\Lambda$:"
    )
    st.latex(
        r"z_1 + z_2 \;=\; z_1 + z_2 \pmod{\Lambda}."
    )
    st.markdown(
        "The identity is $0$, the inverse of $z$ is $-z$, and there is no need for "
        "chord-tangent constructions. The complicated geometry of the algebraic group law "
        "is an artefact of the algebraic presentation — it becomes invisible in the "
        "analytic picture."
    )
    st.markdown(
        "Use the applet below to explore the group law on a torus. A point of "
        "the parallelogram has fractional coordinates $(s, t)$: $z = s \\cdot 1 "
        "+ t \\cdot \\tau$ with $s, t \\in [0, 1)$. **Drag the gold $\\tau$ "
        "corner** to reshape the lattice, then **click** to place $z_1$ and "
        "$z_2$ — their sum is drawn live as honest vector addition, and when it "
        "falls outside the fundamental domain it is translated back in: the "
        "'wrapping around' of the torus."
    )
    components.html(basics_viz.torus_group_html(), height=560, scrolling=False)


# ── Tab 4: Endomorphisms and Complex Multiplication ───────────────────────────
with tab4:
    st.subheader("Endomorphisms and Complex Multiplication")

    # ── Definition ────────────────────────────────────────────────────────────
    st.markdown(
        "An **endomorphism** of an elliptic curve $E$ is a morphism $\\psi: E \\to E$ "
        "that sends the identity point $\\mathcal{O}$ to itself. "
        "Every elliptic curve has an endomorphism for each integer $n$: "
        "the **multiplication-by-$n$ map** $[n]: P \\mapsto P + P + \\cdots + P$ ($n$ times). "
        "The question is whether there are any others."
    )

    # ── Algebraic description ─────────────────────────────────────────────────
    st.markdown("#### On algebraic models")
    st.markdown(
        "On the algebraic model $y^2 = x^3 + fx + g$, an endomorphism takes the form"
    )
    st.latex(r"\psi(x, y) = \bigl(\psi_1(x,y),\; \psi_2(x,y)\bigr),")
    st.markdown(
        "where $\\psi_1$ and $\\psi_2$ are rational functions. "
        "The condition that $\\psi$ maps the curve to itself is that"
    )
    st.latex(
        r"\psi_2(x_0, y_0)^2 = \psi_1(x_0, y_0)^3 + f\,\psi_1(x_0, y_0) + g"
    )
    st.markdown(
        "holds whenever $y_0^2 = x_0^3 + fx_0 + g$, "
        "i.e. whenever $(x_0, y_0)$ lies on $E$."
    )
    st.markdown(
        "**Example — the duplication map $[2]$.** "
        "The chord-tangent construction gives an explicit formula for $[2]$. "
        "At a point $(x, y)$ on $E$, the tangent line has slope "
        "$m = (3x^2 + f)/(2y)$, and"
    )
    st.latex(
        r"[2](x, y) \;=\; \left(\frac{x^4 - 2fx^2 - 8gx + f^2}{4(x^3+fx+g)},\quad "
        r"m\!\left(x - \psi_1(x,y)\right) - y\right)."
    )
    st.markdown(
        "This is a rational function in $x$ and $y$, defined wherever $y \\neq 0$ "
        "(i.e. away from the 2-torsion points). It has **degree 2**: "
        "each point of $E$ has exactly 2 preimages under $[2]$, "
        "counted with multiplicity."
    )

    # ── Analytic description ──────────────────────────────────────────────────
    st.markdown("#### On analytic models")
    st.markdown(
        "The analytic picture makes the structure of endomorphisms completely transparent. "
        "A holomorphic group homomorphism $\\psi: \\mathbb{C}/\\Lambda \\to \\mathbb{C}/\\Lambda$ "
        "must lift to a holomorphic map on $\\mathbb{C}$, and one can show that any "
        "such map is multiplication by a fixed complex number. That is, every endomorphism "
        "has the form"
    )
    st.latex(
        r"[\alpha]_\Lambda : z \pmod{\Lambda} \;\longmapsto\; \alpha z \pmod{\Lambda}"
    )
    st.markdown(
        "for some $\\alpha \\in \\mathbb{C}$ satisfying $\\alpha\\Lambda \\subseteq \\Lambda$. "
        "The condition $\\alpha\\Lambda \\subseteq \\Lambda$ means that $\\alpha$ sends "
        "every lattice point to another lattice point — it must map the generators "
        "$1$ and $\\tau$ of $\\Lambda$ back into $\\Lambda$."
    )
    st.markdown(
        "The integer endomorphisms $[n]$ correspond to $\\alpha = n \\in \\mathbb{Z}$, "
        "which trivially satisfies $n\\Lambda \\subseteq \\Lambda$ for any lattice."
    )

    # ── Endomorphism ring and CM ──────────────────────────────────────────────
    st.markdown("#### The endomorphism ring and complex multiplication")
    st.markdown(
        "The set of all $\\alpha \\in \\mathbb{C}$ with $\\alpha\\Lambda \\subseteq \\Lambda$ "
        "forms a ring under addition and multiplication — the **endomorphism ring** "
        "$\\mathrm{End}(E)$. It always contains $\\mathbb{Z}$."
    )
    st.markdown(
        "For a 'generic' lattice $\\Lambda$, the only $\\alpha$ satisfying "
        "$\\alpha\\Lambda \\subseteq \\Lambda$ are the integers: $\\mathrm{End}(E) = \\mathbb{Z}$. "
        "When $\\mathrm{End}(E)$ is strictly larger than $\\mathbb{Z}$, the curve $E$ "
        "— and the lattice $\\Lambda$ — is said to have **complex multiplication** (CM). "
        "In this case one can show that any non-integer $\\alpha \\in \\mathrm{End}(E)$ "
        "must lie in an imaginary quadratic field $K = \\mathbb{Q}(\\sqrt{-d})$, and "
        "$\\mathrm{End}(E)$ is an order in $K$."
    )
    st.markdown(
        "The two conditions are equivalent: $\\Lambda$ has CM if and only if "
        "$\\mathbb{C}/\\Lambda$ has CM. The CM lattices are precisely those of the form "
        "$\\Lambda = \\mathbb{Z} + \\tau\\mathbb{Z}$ where $\\tau$ is a quadratic "
        "irrational — a root of a quadratic equation with integer coefficients."
    )
    st.divider()
    st.markdown(
        "The applet below illustrates the difference between an endomorphism and a "
        "non-endomorphism. Two lattices are shown: $\\mathbb{Z}[i]$ (square grid) and "
        "$\\mathbb{Z}[2i]$ (rectangular grid with even imaginary parts). "
        "Enter $\\alpha = a + bi$ to see where $\\alpha$ sends each lattice point. "
        "Because $\\mathbb{Z}[i]$ has CM by $\\mathbb{Z}[i]$ itself, every "
        "$\\alpha \\in \\mathbb{Z}[i]$ is an endomorphism — all images land in the lattice "
        "(blue). For $\\mathbb{Z}[2i]$, this only works when $b$ is even; "
        "when $b$ is odd some images fall outside the lattice (red)."
    )

    endo_col, _ = st.columns([1, 3])
    with endo_col:
        a_end = int(st.number_input("a", value=1, step=1, key="endo_a"))
        b_end = int(st.number_input("b", value=1, step=1, key="endo_b"))

    # ── Compute images ────────────────────────────────────────────────────────
    # Source: a patch of each lattice; images may land anywhere
    SRC = 3   # source range: |m|,|n| ≤ SRC

    # Z[i]: z = m + ni;  α·z = (am-bn) + (an+bm)i
    zi_src = [(m, n) for m in range(-SRC, SRC+1) for n in range(-SRC, SRC+1)]
    zi_img = list({(a_end*m - b_end*n, a_end*n + b_end*m) for m, n in zi_src})

    # Z[2i]: z = m + 2ni;  α·z = (am-2bn) + (bm+2an)i
    z2i_src = [(m, n) for m in range(-SRC, SRC+1) for n in range(-SRC//2, SRC//2+1)]
    z2i_good, z2i_bad = [], []
    for m, n in z2i_src:
        rx = a_end*m - b_end*(2*n)
        ry = b_end*m + a_end*(2*n)
        (z2i_good if ry % 2 == 0 else z2i_bad).append((rx, ry))
    z2i_good = list(set(z2i_good))
    z2i_bad  = list(set(z2i_bad))

    # Dynamic display range
    all_coords = [c for pt in zi_img + z2i_good + z2i_bad for c in pt]
    disp_R = max(5, min(14, max((abs(c) for c in all_coords), default=5) + 1))

    # Background lattice points within display window
    BG = int(disp_R) + 1
    zi_bg  = [(m, n)   for m in range(-BG, BG+1) for n in range(-BG, BG+1)
              if abs(m) <= disp_R and abs(n) <= disp_R]
    z2i_bg = [(m, 2*n) for m in range(-BG, BG+1) for n in range(-BG, BG+1)
              if abs(m) <= disp_R and abs(2*n) <= disp_R]

    # ── Plot ──────────────────────────────────────────────────────────────────
    sign_str = f"{b_end:+d}".replace("+", "+").replace("-", "-")
    alpha_str = f"{a_end}{sign_str}i"

    fig_endo, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(10, 5))

    for ax, bg, imgs_blue, imgs_red, title in [
        (ax_l, zi_bg,   zi_img,    [],        f"$\\mathbb{{Z}}[i]$"),
        (ax_r, z2i_bg,  z2i_good,  z2i_bad,   f"$\\mathbb{{Z}}[2i]$"),
    ]:
        ax.scatter([p[0] for p in bg], [p[1] for p in bg],
                   color="gray", alpha=0.25, s=8, zorder=1)
        if imgs_blue:
            ax.scatter([p[0] for p in imgs_blue], [p[1] for p in imgs_blue],
                       color="steelblue", s=30, zorder=3)
        if imgs_red:
            ax.scatter([p[0] for p in imgs_red], [p[1] for p in imgs_red],
                       color="red", s=30, zorder=3)
        ax.axhline(0, color="k", lw=0.4)
        ax.axvline(0, color="k", lw=0.4)
        ax.set_xlim(-disp_R, disp_R)
        ax.set_ylim(-disp_R, disp_R)
        ax.set_aspect("equal")
        ax.set_frame_on(False)
        ax.set_title(f"{title},  $\\alpha = {alpha_str}$", fontsize=11)

    fig_endo.tight_layout()
    st.pyplot(fig_endo)
    plt.close(fig_endo)

    if z2i_bad:
        st.caption(
            f"$\\alpha = {alpha_str}$ is **not** an endomorphism of $\\mathbb{{Z}}[2i]$: "
            f"{len(set(z2i_bad))} image(s) fall outside the lattice (shown in red)."
        )
    else:
        st.caption(
            f"$\\alpha = {alpha_str}$ **is** an endomorphism of both lattices: "
            "all images land back in the lattice."
        )


# ── Tab: The Frobenius endomorphism over 𝔽ₚ ───────────────────────────────────
with tab_frob:
    st.subheader("The Frobenius Endomorphism over $\\mathbb{F}_p$")

    # ── Definition ────────────────────────────────────────────────────────────
    st.markdown(
        "Over $\\mathbb{F}_p$ there is one endomorphism that no other field offers — "
        "and it turns out to control almost everything about the curve. "
        "The **Frobenius endomorphism** is the $p$-th power map on coordinates:"
    )
    st.latex(r"\phi : (x, y) \longmapsto (x^p,\; y^p), \qquad \phi(\mathcal{O}) = \mathcal{O}.")
    st.markdown(
        "Because raising to the $p$-th power respects addition and multiplication "
        "mod $p$ (the 'freshman's dream' $(u+v)^p \\equiv u^p + v^p$), $\\phi$ sends "
        "points of $E$ to points of $E$ and is a genuine endomorphism. "
        "Its **fixed points are exactly the $\\mathbb{F}_p$-rational points**: "
        "$x^p = x$ in $\\mathbb{F}_p$ if and only if $x \\in \\mathbb{F}_p$. So"
    )
    st.latex(r"E(\mathbb{F}_p) = \ker(\phi - 1) = \{P \in E : \phi(P) = P\}.")

    # ── Characteristic equation ───────────────────────────────────────────────
    st.markdown("#### The characteristic equation")
    st.markdown(
        "Frobenius satisfies a single quadratic relation in the endomorphism ring,"
    )
    st.latex(r"\phi^2 - a\,\phi + p = 0,")
    st.markdown(
        "where the integer $a$ is the **trace of Frobenius**. "
        "Counting fixed points of $\\phi - 1$ turns this into a point count: "
        "$\\#E(\\mathbb{F}_p) = \\deg(1 - \\phi) = (1-\\phi)(1-\\bar\\phi) = 1 - a + p$, so"
    )
    st.latex(r"\#E(\mathbb{F}_p) = p + 1 - a.")
    st.markdown(
        "**Hasse's theorem** bounds how far the count can stray from $p+1$:"
    )
    st.latex(r"|a| \le 2\sqrt{p}, \qquad\text{equivalently}\qquad |\,\#E(\mathbb{F}_p) - (p+1)\,| \le 2\sqrt{p}.")
    st.markdown(
        "The trace $a$ is the single most important invariant of the curve: by "
        "**Tate's theorem**, two curves over $\\mathbb{F}_p$ are isogenous if and "
        "only if they have the *same* $a$. So $a$ labels an entire **isogeny class** "
        "at once — exactly the classes this project studies."
    )

    # ── Ordinary vs supersingular ─────────────────────────────────────────────
    st.markdown("#### Frobenius as a quadratic integer")
    st.markdown(
        "Solving $\\phi^2 - a\\phi + p = 0$ gives"
    )
    st.latex(r"\phi = \frac{a + \sqrt{a^2 - 4p}}{2}, \qquad d := a^2 - 4p \le 0.")
    st.markdown(
        "Since $|a| \\le 2\\sqrt p$ the discriminant $d = a^2 - 4p$ is negative "
        "(or zero), so $\\phi$ is an **imaginary quadratic integer**: it generates an "
        "order $\\mathbb{Z}[\\phi]$ inside the imaginary quadratic field "
        "$K = \\mathbb{Q}(\\sqrt{d})$. This is the bridge to the lattice picture — "
        "the same CM orders from the previous tab now appear *forced on us* by "
        "Frobenius. Two cases:\n\n"
        "- **Ordinary** ($a \\not\\equiv 0 \\bmod p$, generically $a \\neq 0$): "
        "$\\mathbb{Z}[\\phi]$ is an order in an imaginary quadratic field and "
        "$\\mathrm{End}(E)$ is commutative.\n"
        "- **Supersingular** ($a = 0$ over $\\mathbb{F}_p$, so $\\phi = \\sqrt{-p}$): "
        "the full endomorphism ring is non-commutative, but restricting to "
        "endomorphisms defined over $\\mathbb{F}_p$ recovers an order in "
        "$\\mathbb{Q}(\\sqrt{-p})$ — so the same machinery still applies."
    )

    st.divider()

    # ── Point-counting applet ─────────────────────────────────────────────────
    st.markdown("#### Count points, read off the trace")
    st.markdown(
        "Pick a prime $p$ and a curve $y^2 = x^3 + fx + g$ (step $f$ and $g$ "
        "with the buttons). The applet counts $\\#E(\\mathbb{F}_p)$ by a "
        "character sum, reads off $a = p + 1 - \\#E$, and places it inside the "
        "Hasse interval. It also runs **backwards**: click any gray dot and a "
        "curve with that trace is hunted down for you — every admissible $a$ is "
        "realized by some curve (Deuring)."
    )
    components.html(basics_viz.hasse_count_html(), height=360, scrolling=False)
    st.caption(
        "Every dot is an integer $a$ with $|a| \\le 2\\sqrt p$; each "
        "corresponds to an isogeny class over $\\mathbb{F}_p$. The red dot is "
        "the class of the current curve."
    )


# ── Tab: Lifting Frobenius (Deuring) ──────────────────────────────────────────
with tab_lift:
    st.subheader("Lifting Frobenius (Deuring)")

    # ── The idea ──────────────────────────────────────────────────────────────
    st.markdown(
        "As we saw on the previous page, many of the most important arithmetic "
        "properties of an elliptic curve $E/\\mathbb{F}_p$ are determined by its "
        "**Frobenius endomorphism** $\\phi$."
    )
    st.markdown(
        "To obtain new pictures of $E$, we will **lift Frobenius** to an analytic "
        "model. Concretely, we want a pair $(\\Lambda, \\alpha)$, where $\\Lambda$ "
        "is a lattice and $\\alpha$ is an endomorphism of $\\Lambda$, such that "
        "$\\mathbb{C}/\\Lambda$ can be described (as an algebraic elliptic curve) by "
        "a model with the following properties:"
    )
    st.markdown(
        "1. the model **descends to $E/\\mathbb{F}_p$**; and\n"
        "2. the endomorphism $\\alpha$ **descends to Frobenius** $\\phi$ on the "
        "curve downstairs.\n\n"
        "We call such a pair $(\\Lambda, \\alpha)$ a **lift of Frobenius**."
    )
    st.markdown(
        "The results of **Deuring** show that a lift of Frobenius can always be "
        "found."
    )

    # ── An explicit example (y^2 = x^3 + 3x over F_5) ─────────────────────────
    st.markdown("#### A worked example: $y^2 = x^3 + 3x$ over $\\mathbb{F}_5$")
    st.markdown(
        "Let us carry out the whole construction once, by hand. (We will never "
        "have to do this again — but it is reassuring to see a lift of Frobenius "
        "explicitly.)"
    )

    st.markdown("**Reading off the pair $(\\Lambda, \\alpha)$.**")
    st.markdown(
        "$E : y^2 = x^3 + 3x$ over $\\mathbb{F}_5$ has **10 points**, so the trace "
        "of Frobenius is $a = p + 1 - \\#E = 6 - 10 = -4$, and Frobenius acts as "
        "multiplication by"
    )
    st.latex(r"\alpha = \frac{a + \sqrt{a^2 - 4p}}{2} = -2 + i, "
             r"\qquad \alpha\bar\alpha = (-2+i)(-2-i) = 5 = p.")
    st.markdown(
        "The lattice is forced to be $\\Lambda = \\mathbb{Z}[i]$: if $-2+i$ is an "
        "endomorphism of $\\Lambda$, then so is $(-2+i) + 2 = i$ — and $i$ is an "
        "**automorphism of order 4** ($i^4 = 1$). But $\\mathbb{Z}[i]$ is the only "
        "lattice (up to scaling) admitting an order-4 automorphism. So our "
        "candidate lift is $(\\Lambda, \\alpha) = (\\mathbb{Z}[i],\\, -2+i)$."
    )

    st.markdown("**The model.**")
    st.markdown(
        "The curve $y^2 = x^3 + 3x$ does the job. Over $\\mathbb{C}$ it is the "
        "analytic curve $\\mathbb{C}/\\mathbb{Z}[i]$ (it has $j = 1728$, with CM by "
        "$\\mathbb{Z}[i]$), and its equation **descends to $E/\\mathbb{F}_5$** by "
        "reducing the coefficients (already integers) mod $5$. It remains to check "
        "that $\\alpha = -2+i$ **descends to Frobenius**."
    )

    st.markdown("**Checking that $\\alpha$ descends to Frobenius.**")
    st.markdown(
        "1. Multiplication by $i$ on the model is $(x, y) \\mapsto (-x,\\, iy)$.\n"
        "2. Using the group law, $[-2+i] = [-2] + [i]$ gives an explicit (very "
        "large) rational map $F(x, y)$ — shown below.\n"
        "3. Reduce modulo the prime ideal $\\mathfrak{p} = (-2+i)$, for which "
        "$\\mathbb{Z}[i]/\\mathfrak{p} \\cong \\mathbb{F}_5$ and $i \\equiv 2$: "
        "replace every $i$ by $2$ and reduce coefficients mod $5$.\n"
        "4. After simplifying (using $y^2 = x^3 + 3x$), the map collapses to"
    )
    st.latex(r"F(x, y) \equiv (x^5,\; y^5) \pmod{\mathfrak{p}},")
    st.markdown(
        "which is exactly the Frobenius endomorphism. So $(\\mathbb{Z}[i],\\, -2+i)$ "
        "really is a lift of Frobenius for $E/\\mathbb{F}_5$."
    )

    with st.expander("The explicit lift map $F(x,y)$  (never actually needed)"):
        st.markdown(
            "It is striking how large the lift is. Even just the $x$-coordinate of "
            "$F(x, y)$ is"
        )
        st.latex(
            r"\scriptsize "
            r"-\frac{\left(x^2-3\right)^2 \left(x^8+(36-96 i) x^6+(342-576 i) x^4"
            r"+(324-864 i) x^2+81\right)^2}"
            r"{4 x \left(x^2+3\right) \left(x^2+(3-6 i)\right)^2 "
            r"\left((2+i) x^2+3 i\right)^2 \left(x^4-(6-24 i) x^2+9\right)^2}"
        )
        st.caption(
            "The $y$-coordinate is larger still (it carries a degree-16 factor), so "
            "we show only $x$. In practice we never compute this map at all — it "
            "just sits in the background, guaranteeing the lift exists."
        )

    st.markdown(
        "Going forward we will **not** carry out this computation. Finding the "
        "lattice–endomorphism pair $(\\Lambda, \\alpha)$ for each $E/\\mathbb{F}_p$ "
        "is done with more powerful tools, explained in the coming lessons."
    )

    st.divider()

    # ── Applet: Frobenius acting on the CM lattice ────────────────────────────
    st.markdown("#### Frobenius as multiplication by $\\alpha$")
    st.markdown(
        "Choose $(a, p)$ with $a^2 < 4p$ by **clicking a dot in the picker** — "
        "one row per small prime, the admissible traces fanned out under the "
        "Hasse parabola. We take the lattice $\\Lambda = \\mathbb{Z} + "
        "\\alpha\\,\\mathbb{Z}$ (the order $\\mathcal{O} = "
        "\\mathbb{Z}[\\alpha]$ itself) and watch multiplication by $\\alpha$ "
        "send its generators to other lattice points — the analytic shadow of "
        "Frobenius. On the basis $(1, \\alpha)$, multiplication by $\\alpha$ "
        "has the integer matrix (from $\\alpha\\cdot 1 = \\alpha$ and "
        "$\\alpha^2 = a\\alpha - p$):"
    )
    st.latex(r"[\alpha] = \begin{pmatrix} 0 & -p \\ 1 & a \end{pmatrix},")
    st.markdown(
        "and integer entries mean exactly that $\\alpha\\Lambda \\subseteq "
        "\\Lambda$: Frobenius really is an endomorphism of the lattice. (Small "
        "primes keep the picture compact — multiplication by $\\alpha$ scales "
        "by $\\sqrt p$, so the image cell flies off-screen for large $p$.)"
    )
    components.html(basics_viz.frobenius_lift_html(), height=530, scrolling=False)


# ── Tab: Pictures from lifts of Frobenius — the multiplicative group ───────────
with tab_mult:
    st.subheader("Pictures from Lifts of Frobenius")

    # ── The general recipe ────────────────────────────────────────────────────
    st.markdown(
        "A lift of Frobenius is exactly the data we need to **draw** a group over "
        "$\\mathbb{F}_p$. Suppose we have:"
    )
    st.markdown(
        "- an algebraic group $E/\\mathbb{F}_p$ we want to visualize;\n"
        "- an algebraic group $X/\\mathbb{C}$ we *can* visualize;\n"
        "- a model of $X$ over the algebraic integers that **descends to "
        "$E/\\mathbb{F}_p$** modulo a prime ideal $P$; and\n"
        "- an endomorphism $F : X \\to X$ that **descends to Frobenius** on "
        "$E/\\mathbb{F}_p$."
    )
    st.markdown(
        "Then we draw $E$ inside the ambient space $X(\\mathbb{C})$ by plotting the "
        "**fixed points of $F$**,"
    )
    st.latex(r"\mathrm{Fix}(F) = \{\,x \in X(\mathbb{C}) : F(x) = x\,\},")
    st.markdown(
        "which are precisely the $\\mathbb{F}_p$-rational points of $E$ (Frobenius "
        "fixes exactly what is defined over $\\mathbb{F}_p$). Better still, the "
        "fixed points of the **powers** $F^n$ are the points defined over the "
        "extension $\\mathbb{F}_{p^n}$:"
    )
    st.latex(r"\mathrm{Fix}(F^n) = E(\mathbb{F}_{p^n}).")
    st.markdown(
        "So one picture, drawn in $X(\\mathbb{C})$, can hold the points of $E$ over "
        "$\\mathbb{F}_p$ **and** over many extensions at once. And because $F$ "
        "*is* the Frobenius, we can **watch the Galois action**: applying $F$ "
        "permutes the points, and the length of a point's $F$-orbit is the degree "
        "of the smallest field over which it is defined."
    )

    # ── The multiplicative group instance ─────────────────────────────────────
    st.markdown("#### Warm-up: the multiplicative group")
    st.markdown(
        "Before elliptic curves, the same recipe draws the multiplicative group "
        "$\\mathbb{G}_m$ — that is, the groups $\\mathbb{F}_q^*$. Take"
    )
    st.latex(r"X = \mathbb{C}^*, \qquad F : x \mapsto x^p,")
    st.markdown(
        "which descends to the Frobenius $x \\mapsto x^p$ on $\\mathbb{F}_p$. Its "
        "fixed points solve $x^p = x$ with $x \\neq 0$ — the $(p-1)$-th roots of "
        "unity, i.e. $\\mathbb{F}_p^*$. More generally"
    )
    st.latex(r"\mathrm{Fix}(F^n) = \{x \neq 0 : x^{p^n} = x\} "
             r"= \{(p^n-1)\text{-th roots of unity}\} = \mathbb{F}_{p^n}^*,")
    st.markdown(
        "all sitting on the **unit circle**. Frobenius acts by $x \\mapsto x^p$ — "
        "it *multiplies the angle by $p$*. The $p-1$ points of $\\mathbb{F}_p^*$ "
        "are fixed; the rest fall into orbits whose length is the degree of the "
        "extension they generate."
    )

    st.divider()

    # ── Applet: roots of unity, coloured by degree, with Frobenius orbits ──────
    st.markdown("#### The picture")
    st.markdown(
        "Choose a prime $p$ and a top degree $n$. We plot $\\mathbb{F}_{p^n}^*$ as "
        "the $(p^n-1)$-th roots of unity, coloured by the field each point is "
        "defined over, and draw the full Frobenius action: an **arrow "
        "$x \\mapsto x^p$ from every point**, with a **loop on each fixed point** "
        "(the elements of $\\mathbb{F}_p^*$). Pick a point to emphasise its orbit — "
        "the number of hops back to the start is the degree of its field."
    )

    mult_ctrl, mult_plot = st.columns([1, 2])

    with mult_ctrl:
        mp = st.selectbox("p", [2, 3, 5, 7], index=1, key="mult_p")
        # Cap the number of points so the circle stays readable; n's range
        # depends on p.  (p-keyed slider so switching p never leaves a stale n.)
        n_max = 1
        while mp ** (n_max + 1) - 1 <= 160:
            n_max += 1
        mn = st.slider("top degree n  (show $\\mathbb{F}_{p^n}$)",
                       1, n_max, min(2, n_max), key=f"mult_n_{mp}")
        M = mp ** mn - 1
        st.caption(
            f"$\\#\\mathbb{{F}}_{{{mp}^{{{mn}}}}}^* = {mp}^{{{mn}}} - 1 = {M}$ "
            f"points; $\\#\\mathbb{{F}}_{{{mp}}}^* = {mp - 1}$ are fixed by $F$."
        )
        if M > 1:
            k_hi = st.slider(
                "highlight orbit of point $k$ (exponent of a generator)",
                0, M - 1, min(1, M - 1), key=f"mult_k_{mp}_{mn}")
        else:
            k_hi = 0
            st.caption("Only one point: $\\mathbb{F}_2^* = \\{1\\}$.")

    with mult_plot:
        # Each F_{p^n}^* element is g^k for a fixed generator g; we draw it at the
        # (p^n-1)-th root of unity ζ^k.  Frobenius x↦x^p acts as k ↦ p·k mod M.
        def _orbit(k0):
            seen = [k0]
            k = (mp * k0) % M
            while k != k0:
                seen.append(k)
                k = (mp * k) % M
            return seen

        deg = {}                       # k -> field-of-definition degree
        for k in range(M):
            deg[k] = len(_orbit(k))

        divisors = [d for d in range(1, mn + 1) if mn % d == 0]
        pal = {1: "crimson"}
        others = [d for d in divisors if d != 1]
        ocols = plt.cm.viridis(np.linspace(0.25, 0.85, max(1, len(others))))
        for i, d in enumerate(others):
            pal[d] = ocols[i]

        def _xy(k):
            t = 2 * np.pi * k / M
            return np.cos(t), np.sin(t)

        fig_m, ax_m = plt.subplots(figsize=(5.5, 5.5))
        circ = np.linspace(0, 2 * np.pi, 400)
        ax_m.plot(np.cos(circ), np.sin(circ), color="lightgray", lw=1, zorder=1)

        orb = _orbit(k_hi)
        orb_set = set(orb)

        # The full Frobenius action x ↦ x^p: an arrow from every point to its
        # image; fixed points (x^p = x) carry a loop.  The selected orbit is
        # emphasised in black, the rest drawn faint.
        for k in range(M):
            kk = (mp * k) % M
            hot = k in orb_set
            acol = "black" if hot else "0.6"
            alw = 1.5 if hot else 0.8
            az = 4 if hot else 2
            if kk == k:                                   # fixed point → loop
                t = 2 * np.pi * k / M
                rad = 0.09
                cx, cy = (1 + rad) * np.cos(t), (1 + rad) * np.sin(t)
                ax_m.add_patch(MplArc((cx, cy), 2 * rad, 2 * rad, angle=0,
                                      theta1=0, theta2=360,
                                      color=acol, lw=alw, zorder=az))
            else:
                xa, ya = _xy(k)
                xb, yb = _xy(kk)
                ax_m.annotate("", xy=(xb, yb), xytext=(xa, ya),
                              arrowprops=dict(arrowstyle="->", color=acol,
                                              lw=alw, alpha=0.85,
                                              connectionstyle="arc3,rad=0.18"),
                              zorder=az)

        # points coloured by degree (selected orbit ringed in black)
        for k in range(M):
            x, y = _xy(k)
            d = deg[k]
            ax_m.scatter(x, y, color=pal[d], s=(95 if d == 1 else 55),
                         edgecolor=("black" if k in orb_set else "white"),
                         linewidth=(1.5 if k in orb_set else 0.6), zorder=5)

        ax_m.set_xlim(-1.25, 1.25)
        ax_m.set_ylim(-1.25, 1.25)
        ax_m.set_aspect("equal")
        ax_m.set_frame_on(False)
        ax_m.set_xticks([]); ax_m.set_yticks([])
        # degree legend
        handles = [plt.Line2D([0], [0], marker="o", linestyle="",
                              markerfacecolor=pal[d], markeredgecolor="white",
                              markersize=9,
                              label=(f"$\\mathbb{{F}}_{{{mp}}}^*$ (deg 1, fixed)"
                                     if d == 1 else
                                     f"deg {d}  ($\\mathbb{{F}}_{{{mp}^{{{d}}}}}$)"))
                   for d in divisors]
        ax_m.legend(handles=handles, loc="upper left", fontsize=8,
                    frameon=False, bbox_to_anchor=(-0.02, 1.02))
        ax_m.set_title("$\\mathbb{F}_{%d^%d}^*$ on the unit circle" % (mp, mn),
                       fontsize=11)
        st.pyplot(fig_m)
        plt.close(fig_m)

        d_sel = deg[k_hi]
        st.caption(
            f"Highlighted orbit has length {d_sel}, so $\\zeta^{{{k_hi}}}$ is "
            f"defined over $\\mathbb{{F}}_{{{mp}^{{{d_sel}}}}}$ "
            + ("(a fixed point — it lies in $\\mathbb{F}_p^*$)." if d_sel == 1
               else f"and nowhere smaller. Applying $F$ {d_sel} times returns it "
                    "to the start.")
        )


# ── Tab: Pictures from lifts of Frobenius — elliptic curves ───────────────────
with tab_ell:
    st.subheader("Pictures from Lifts of Frobenius: Elliptic Curves")

    st.markdown(
        "We now have both halves of the story. In **§3.1** we built the recipe — "
        "plot the fixed points of $F$ (and its powers) inside an ambient "
        "$X(\\mathbb{C})$ — and in **§2** we obtained the ingredients for an "
        "elliptic curve: a CM lattice $\\Lambda$ together with the complex number "
        "$\\alpha$ that Frobenius lifts to. Putting them together gives a new kind "
        "of picture, drawn inside a **fundamental domain of the torus "
        "$\\mathbb{C}/\\Lambda$**."
    )
    st.markdown(
        "Here $X = \\mathbb{C}/\\Lambda$ and $F$ is multiplication by $\\alpha$, so "
        "the points of $E(\\mathbb{F}_{p^n})$ are the solutions of "
        "$\\alpha^n z \\equiv z$, i.e. $(\\alpha^n - 1)\\,z \\in \\Lambda$ — a "
        "finite set of $|\\alpha^n - 1|^2 = \\#E(\\mathbb{F}_{p^n})$ points in the "
        "parallelogram. These are exactly the **lattice pictures** on the "
        "*Elliptic Curve Search* page."
    )

    st.divider()

    st.markdown("#### Classical picture vs. lattice picture")
    st.markdown(
        "Below, the same curve is drawn both ways. We use the four curves with "
        "$j = 1728$ over $\\mathbb{F}_5$ — the twists $y^2 = x^3 + a x$ — which all "
        "have $\\Lambda = \\mathbb{Z}[i]$ and Frobenius $\\alpha$ an associate of "
        "$-2+i$."
    )

    # ── Gaussian-integer helpers ──────────────────────────────────────────────
    def _gmul(z, w):
        return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])

    def _gpow(z, n):
        r = (1, 0)
        for _ in range(n):
            r = _gmul(r, z)
        return r

    def _fixed_pts(beta):
        """Points z in [0,1)^2 of C/Z[i] with beta*z in Z[i]; returns N of them."""
        br, bi = beta
        N = br * br + bi * bi
        seen = set()
        for a in range(N):
            for b in range(N):
                seen.add(((a * br + b * bi) % N, (b * br - a * bi) % N))
        return [(rx / N, ry / N) for rx, ry in seen], N

    def _fmt_gauss(z):
        r, i = z
        if i == 0:
            return str(r)
        sign = "+" if i > 0 else "-"
        mag = abs(i)
        coef = "" if mag == 1 else str(mag)
        return f"{r} {sign} {coef}i"

    _ALPHAS = {1: (1, 2), 2: (2, -1), 3: (-2, 1), 4: (-1, -2)}

    ec_ctrl, ec_plot = st.columns([1, 2])

    with ec_ctrl:
        a_ec = st.select_slider("curve  $y^2 = x^3 + a\\,x$  (mod 5)",
                                options=[1, 2, 3, 4], value=3, key="ec_a")
        alpha = _ALPHAS[a_ec]
        beta1 = (alpha[0] - 1, alpha[1])
        _, N1 = _fixed_pts(beta1)
        st.latex(rf"\alpha = {_fmt_gauss(alpha)}")
        st.latex(rf"\#E(\mathbb{{F}}_5) = |\alpha-1|^2 = {N1}")
        show25 = st.checkbox("also show $\\mathbb{F}_{25}$ points on the lattice",
                             value=False, key="ec_25")
        st.caption(
            "The classical model can't show the $\\mathbb{F}_{25}$-points at all — "
            "they don't live in the $\\mathbb{F}_5$ grid."
        )

    with ec_plot:
        p = 5

        def _sym(v):                          # symmetric representative in F_5
            return v if 2 * v < p else v - p

        # classical affine points (symmetric coords -2..2)
        aff = [(_sym(x), _sym(y)) for x in range(p) for y in range(p)
               if (y * y) % p == (x * x * x + a_ec * x) % p]
        amb = [(_sym(x), _sym(y)) for x in range(p) for y in range(p)]
        # lattice points
        f5_pts, _ = _fixed_pts(beta1)
        if show25:
            beta2 = (_gpow(alpha, 2)[0] - 1, _gpow(alpha, 2)[1])
            f25_pts, N2 = _fixed_pts(beta2)
            f5_set = set(f5_pts)
            f25_extra = [q for q in f25_pts if q not in f5_set]

        fig_ec, (axC, axL) = plt.subplots(1, 2, figsize=(9, 4.6))

        # ── classical ─────────────────────────────────────────────────────────
        h = p // 2
        axC.scatter([q[0] for q in amb], [q[1] for q in amb],
                    color="gray", alpha=0.2, s=18, zorder=1)
        axC.scatter([q[0] for q in aff], [q[1] for q in aff],
                    color="steelblue", s=55, zorder=3)
        axC.set_xlim(-h - 0.6, h + 0.6); axC.set_ylim(-h - 0.6, h + 0.6)
        axC.set_aspect("equal"); axC.set_frame_on(False)
        axC.set_xticks(range(-h, h + 1)); axC.set_yticks(range(-h, h + 1))
        axC.tick_params(labelsize=7)
        axC.set_title(f"Classical: $y^2=x^3+{a_ec}x$ over $\\mathbb{{F}}_5$\n"
                      f"({len(aff)} affine points $+\\,\\mathcal{{O}}$ at $\\infty$)",
                      fontsize=9)

        # ── lattice ───────────────────────────────────────────────────────────
        axL.add_patch(MplPolygon([(0, 0), (1, 0), (1, 1), (0, 1)], closed=True,
                                 facecolor="gray", alpha=0.12, edgecolor="none",
                                 zorder=0))
        if show25:
            axL.scatter([q[0] for q in f25_extra], [q[1] for q in f25_extra],
                        color="goldenrod", s=40, zorder=2,
                        label=f"$\\mathbb{{F}}_{{25}}$ ({N2} pts)")
        axL.scatter([q[0] for q in f5_pts], [q[1] for q in f5_pts],
                    color="steelblue", s=55, zorder=3,
                    label=f"$\\mathbb{{F}}_5$ ({N1} pts)")
        # the zero point at the origin
        axL.scatter([0], [0], color="crimson", s=90, zorder=5)
        axL.annotate("$0$", (0, 0), xytext=(6, 6), textcoords="offset points",
                     color="crimson", fontsize=12, fontweight="bold")
        axL.set_xlim(-0.12, 1.12); axL.set_ylim(-0.12, 1.12)
        axL.set_aspect("equal"); axL.set_frame_on(False)
        axL.set_xticks([0, 1]); axL.set_yticks([0, 1])
        axL.tick_params(labelsize=7)
        ttl = ("Lattice: fixed points of $\\times\\alpha$ on $\\mathbb{C}/\\mathbb{Z}[i]$"
               + ("\n($\\mathbb{F}_5$ and $\\mathbb{F}_{25}$ together)" if show25
                  else "\n($0$ is visible, at the corner)"))
        axL.set_title(ttl, fontsize=9)
        if show25:
            axL.legend(loc="upper right", fontsize=7, frameon=False)

        fig_ec.tight_layout()
        st.pyplot(fig_ec)
        plt.close(fig_ec)

    st.markdown(
        "Superficially the two pictures are alike — each is a discrete set of "
        "points in a fundamental domain of a torus. But the lattice picture has "
        "real advantages:"
    )
    st.markdown(
        "- **The zero point is visible.** On the algebraic model the identity "
        "$\\mathcal{O}$ is the point at infinity, off the chart; in the lattice "
        "picture it is simply $0$, at the corner of the parallelogram.\n"
        "- **The group structure is clearer.** The points are the group "
        "$\\Lambda/(\\alpha-1)\\Lambda$, laid out regularly in the square.\n"
        "- **Extension points fit in the same picture.** The points over "
        "$\\mathbb{F}_{p^2}, \\mathbb{F}_{p^3}, \\ldots$ are just the fixed points "
        "of $\\alpha^2, \\alpha^3, \\ldots$ — more dots in the *same* fundamental "
        "domain, no need to enlarge the ambient space (toggle $\\mathbb{F}_{25}$ "
        "above). The classical $\\mathbb{F}_5$ grid simply has nowhere to put them."
    )

    st.divider()

    st.markdown("#### Explore the lattice picture")
    st.markdown(
        "The applet above is pinned to the four $j = 1728$ twists over "
        "$\\mathbb{F}_5$. Here is the same construction set loose: pick any prime "
        "$p$, any trace $a$ (so $\\#E(\\mathbb{F}_p) = p + 1 - a$), and any lattice "
        "class of discriminant $a^2 - 4p$. The dots are the "
        "$N = \\#\\,\\mathrm{Fix}([\\alpha]) = \\chi(1)$ solutions of "
        "$(\\alpha - 1)\\,z \\in \\Lambda$ inside the fundamental domain "
        "$\\langle 1, \\tau\\rangle$ — i.e. the group $E(\\mathbb{F}_p)$ drawn on "
        "the CM torus. Every lattice class of a given discriminant carries the "
        "*same* number of points; they are the different curves in one isogeny "
        "class."
    )
    components.html(slide_viz.cm_torus_html(), height=440, scrolling=False)


# ── Tab: Isogenies — kernels and degree ───────────────────────────────────────
with tab_isog:
    st.subheader("Isogenies: Kernels and Degree")

    # ── Definition ────────────────────────────────────────────────────────────
    st.markdown(
        "An **isogeny** is a non-constant morphism $\\varphi : E \\to E'$ between "
        "two elliptic curves that sends the identity to the identity, "
        "$\\varphi(\\mathcal{O}) = \\mathcal{O}'$. A basic theorem says that this "
        "single condition forces $\\varphi$ to be a **group homomorphism** "
        "automatically — isogenies respect the group law for free."
    )
    st.markdown(
        "Every isogeny has a finite **kernel** $\\ker\\varphi = "
        "\\varphi^{-1}(\\mathcal{O}')$, and (for the separable isogenies we care "
        "about) its **degree** is exactly the size of that kernel:"
    )
    st.latex(r"\deg \varphi = \#\ker\varphi.")
    st.markdown(
        "The kernel determines the isogeny: for **every** finite subgroup "
        "$C \\subseteq E$ there is an essentially unique isogeny with that kernel, "
        "the quotient map"
    )
    st.latex(r"\varphi : E \longrightarrow E/C, \qquad \deg\varphi = \#C.")
    st.markdown(
        "An isogeny of degree $\\ell$ is called an **$\\ell$-isogeny**; its kernel "
        "is an order-$\\ell$ subgroup of the $\\ell$-torsion $E[\\ell]$. "
        "Isogeny is an *equivalence relation* — every $\\varphi$ has a **dual** "
        "$\\hat\\varphi$ with $\\hat\\varphi \\circ \\varphi = [\\deg\\varphi]$ — "
        "and isogenous curves over $\\mathbb{F}_p$ share the same trace $a$, hence "
        "the same point count."
    )

    # ── Two descriptions ──────────────────────────────────────────────────────
    st.markdown("#### Two descriptions of the same map")
    st.markdown(
        "Just as with the group law and with endomorphisms, an isogeny looks very "
        "different in the algebraic and analytic pictures:"
    )
    st.markdown(
        "**Algebraically** ($E : y^2 = x^3 + fx + g$), an isogeny is a pair of "
        "rational functions $\\varphi(x, y) = \\bigl(r(x),\\, y\\,s(x)\\bigr)$ "
        "carrying one Weierstrass equation to another. Writing these down from a "
        "kernel is the content of **Vélu's formulas** (next tab). The 2-isogenies "
        "have the smallest such formulas; higher degrees grow quickly."
    )
    st.markdown(
        "**Analytically** ($E = \\mathbb{C}/\\Lambda$), an isogeny is induced by a "
        "complex number $\\beta$ with $\\beta\\Lambda \\subseteq \\Lambda'$ — and "
        "in the cleanest form, by simply taking a **finer lattice**: any sublattice "
        "relationship $\\Lambda \\subseteq \\Lambda'$ of finite index gives"
    )
    st.latex(r"\varphi : \mathbb{C}/\Lambda \longrightarrow \mathbb{C}/\Lambda', "
             r"\qquad \ker\varphi = \Lambda'/\Lambda, \qquad "
             r"\deg\varphi = [\Lambda' : \Lambda].")
    st.markdown(
        "So an $\\ell$-isogeny is just an **index-$\\ell$ overlattice** "
        "$\\Lambda \\subseteq \\Lambda'$. This is dramatically more transparent than "
        "the algebraic version — and it has a vivid geometric picture: the isogeny "
        "**folds** the torus $\\mathbb{C}/\\Lambda$ onto the smaller torus "
        "$\\mathbb{C}/\\Lambda'$, exactly $\\ell$-to-$1$. The next tab is an applet "
        "for that folding."
    )

    # ── Applet: isogenies are determined by their kernels ──────────────────────
    st.markdown("#### Isogenies are determined by their kernels")
    st.markdown(
        "Every finite subgroup $C \\subseteq E$ is the kernel of an isogeny, and "
        "$C$ pins down the codomain $E/C$. On the left is a fundamental domain "
        "of $\\Lambda$ with the $d \\times d$ grid of **$d$-torsion points** "
        "$E[d]$: **drag the gold $\\tau$ corner** to reshape the lattice, pick "
        "the prime degree $d$ with the buttons, and **click a nonzero point** to "
        "choose the generator $P$ of an order-$d$ kernel $C = \\langle P"
        "\\rangle$. The codomain $E/C = \\mathbb{C}/\\Lambda'$ on the right, "
        "with $\\Lambda' = \\Lambda + \\mathbb{Z}P$, follows along live."
    )
    components.html(basics_viz.isogeny_kernel_html(), height=560, scrolling=False)
    st.caption(
        "Left: the $d$-torsion $E[d]$ in a fundamental domain of $\\Lambda$; the "
        "chosen generator $P$ is ringed and its subgroup $\\langle P\\rangle = "
        "\\{O, P, \\ldots, (d{-}1)P\\}$ is red. Right: the codomain lattice "
        "$\\Lambda' = \\Lambda + \\mathbb{Z}P$ — the original lattice $\\Lambda$ "
        "(gray) plus the cosets the kernel adds (gold). The green parallelogram is "
        "a fundamental domain of $\\Lambda'$; its area is $1/d$ of the original, "
        "matching the $d$ points of the kernel. A different kernel gives a "
        "different $\\Lambda'$."
    )


# ── Tab: Analytic pictures — folding the torus ────────────────────────────────
with tab_fold:
    st.subheader("Analytic Pictures: Folding the Torus")

    st.markdown(
        "Here is the payoff of the analytic picture. An $\\ell$-isogeny "
        "$\\varphi : \\mathbb{C}/\\Lambda \\to \\mathbb{C}/\\Lambda'$ comes from an "
        "index-$\\ell$ overlattice $\\Lambda \\subseteq \\Lambda'$, and the map is "
        "just 'reduce mod the finer lattice'. Geometrically it cuts the fundamental "
        "domain of $\\Lambda$ into $\\ell$ congruent slabs and **folds** them onto "
        "the smaller domain of $\\Lambda'$ — so $\\varphi$ is exactly "
        "$\\ell$-to-$1$, and the $\\ell$ preimages of any point differ by the "
        "kernel."
    )
    st.markdown(
        "Pick a degree $\\ell$ and a cyclic kernel direction with the buttons; "
        "**drag $P$** anywhere in the domain (and the gold $\\tau$ corner to "
        "reshape the lattice). The whole fibre of $P$ — its $\\ell$ "
        "kernel-translates, one per slab — moves with it on the left, and all "
        "of them land on the single image $\\varphi(P)$ on the right."
    )
    components.html(basics_viz.torus_folding_html(), height=560, scrolling=False)
    st.caption(
        "Left: the $\\ell$ coloured slabs of $\\mathbb{C}/\\Lambda$, the kernel "
        "points (red), and the fibre of $P$ — its $\\ell$ kernel-translates "
        "(white). Right: all $\\ell$ of them fold onto the single image "
        "$\\varphi(P)$ in the smaller torus $\\mathbb{C}/\\Lambda'$."
    )


# ── Tab: Vélu's formulas over 𝔽ₚ ──────────────────────────────────────────────
with tab_velu:
    st.subheader("Vélu's Formulas over $\\mathbb{F}_p$")

    st.markdown(
        "The analytic 'finer lattice' picture is beautiful but lives over "
        "$\\mathbb{C}$. To actually *compute* an isogeny over $\\mathbb{F}_p$ — "
        "which is what this project needs — we use **Vélu's formulas**: given a "
        "curve $E$ and a finite kernel subgroup $C$, they write down the codomain "
        "$E/C$ and the rational map explicitly, with no lattice in sight."
    )

    # ── The recipe ────────────────────────────────────────────────────────────
    st.markdown("#### The recipe")
    st.markdown(
        "Start from $E : y^2 = x^3 + fx + g$ and a kernel $C$. Pick a set $S$ of "
        "kernel points with one representative from each pair $\\{Q, -Q\\}$ "
        "(2-torsion points are their own inverse). For each $Q = (x_Q, y_Q) \\in S$ set"
    )
    st.latex(r"""
        g^x_Q = 3x_Q^2 + f, \qquad u_Q = (2y_Q)^2, \qquad
        v_Q = \begin{cases} g^x_Q & 2Q = \mathcal{O} \\ 2g^x_Q & \text{otherwise,}\end{cases}
    """)
    st.markdown("and accumulate $v = \\sum_Q v_Q$, $w = \\sum_Q (u_Q + x_Q v_Q)$. Then the **codomain** is")
    st.latex(r"E/C : \; Y^2 = X^3 + (f - 5v)\,X + (g - 7w),")
    st.markdown("and the **isogeny** $\\varphi(x, y) = (X, Y)$ is")
    st.latex(r"""
        X = x + \sum_{Q \in S}\!\left[\frac{v_Q}{x - x_Q} + \frac{u_Q}{(x-x_Q)^2}\right],
        \qquad
        Y = y\,\frac{dX}{dx},
    """)
    st.markdown(
        "the last equation because Vélu's map is *normalised*: it pulls the "
        "invariant differential $dx/y$ back to $dX/Y$."
    )

    # ── The 2-isogeny in closed form ──────────────────────────────────────────
    st.markdown("#### The simplest case: a 2-isogeny")
    st.markdown(
        "A 2-torsion point is $Q = (x_0, 0)$ with $x_0$ a root of $x^3 + fx + g$. "
        "Here $u_Q = 0$ and $v = g^x_Q = 3x_0^2 + f$, and everything collapses to"
    )
    st.latex(r"""
        E/C : \; Y^2 = X^3 + (f - 5v)X + (g - 7x_0 v), \qquad v = 3x_0^2 + f,
    """)
    st.latex(r"""
        X = x + \frac{v}{x - x_0}, \qquad Y = y\left(1 - \frac{v}{(x-x_0)^2}\right).
    """)
    st.markdown("Higher degrees follow the same recipe with more kernel points.")

    st.divider()

    # ── Numeric applet ────────────────────────────────────────────────────────
    st.markdown("#### Compute one over $\\mathbb{F}_p$")
    st.markdown(
        "Choose a curve and a degree. We find a kernel of that order, run Vélu, and "
        "show both curves over $\\mathbb{F}_p$ with a point and its image."
    )

    _VPRIMES = [p for p in range(11, 60)
                if all(p % d for d in range(2, int(p**0.5) + 1))]

    def _v_sym(x, p):
        xr = x % p
        return xr if 2 * xr < p else xr - p

    def _v_ec_pts(f, g, p):
        out = []
        for x in range(p):
            r = (pow(x, 3, p) + f * x + g) % p
            for y in range(p):
                if (y * y) % p == r:
                    out.append((x, y))
        return out

    def _v_kernel(f, g, p, ell):
        """Return (S, two_torsion_index_set) or None."""
        if ell == 2:
            roots = [x for x in range(p) if (pow(x, 3, p) + f * x + g) % p == 0]
            if not roots:
                return None
            return [(roots[0], 0)], {0}
        # ell == 3: order-3 points are roots of the 3-division polynomial
        for x in range(p):
            if (3 * pow(x, 4, p) + 6 * f * x * x + 12 * g * x - f * f) % p == 0:
                r = (pow(x, 3, p) + f * x + g) % p
                for y in range(1, p):
                    if (y * y) % p == r:
                        return [(x, y)], set()
        return None

    def _v_run(f, g, p, S, tt):
        v = w = 0
        data = []
        for i, (xq, yq) in enumerate(S):
            gx = (3 * xq * xq + f) % p
            uq = (4 * yq * yq) % p
            vq = gx % p if i in tt else (2 * gx) % p
            v = (v + vq) % p
            w = (w + uq + xq * vq) % p
            data.append((xq, yq, uq, vq))
        f2, g2 = (f - 5 * v) % p, (g - 7 * w) % p

        def phi(x, y):
            X = x % p
            corr = 0
            for (xq, yq, uq, vq) in data:
                inv = pow((x - xq) % p, -1, p)
                X = (X + vq * inv + uq * inv * inv) % p
                corr = (corr + vq * inv * inv + 2 * uq * pow(inv, 3, p)) % p
            return X, (y * ((1 - corr) % p)) % p

        return f2, g2, phi

    v_ctrl, v_plot = st.columns([1, 2])

    with v_ctrl:
        vp = st.selectbox("p", _VPRIMES, index=_VPRIMES.index(23), key="velu_p")
        vf = int(st.number_input("f", value=1, step=1, key="velu_f"))
        vg = int(st.number_input("g", value=1, step=1, key="velu_g"))
        vell = st.select_slider("degree ℓ", options=[2, 3], value=2, key="velu_ell")

        vdisc = (-16 * (4 * pow(vf, 3) + 27 * pow(vg, 2))) % vp
        velu_ok = False
        if vdisc == 0:
            st.warning("Singular curve mod p — adjust f or g.")
        else:
            ker = _v_kernel(vf, vg, vp, vell)
            if ker is None:
                st.warning(
                    f"No kernel of order {vell} over $\\mathbb{{F}}_{{{vp}}}$ for "
                    "this curve — try another $(f, g, p)$ or degree."
                )
            else:
                S, tt = ker
                f2, g2, phi = _v_run(vf, vg, vp, S, tt)
                velu_ok = True
                st.success("Smooth curve ✓  kernel found ✓")
                st.markdown("**Codomain $E/C$**")
                st.latex(rf"Y^2 = X^3 + {f2}\,X + {g2} \pmod{{{vp}}}")
                xq0 = S[0][0]
                st.caption(
                    f"kernel generator at $x = {_v_sym(xq0, vp)}$ "
                    + ("(2-torsion, $y=0$)" if vell == 2
                       else f", $y = {_v_sym(S[0][1], vp)}$")
                )

    with v_plot:
        if not velu_ok:
            st.info("Pick a curve with a kernel of the chosen order.")
        else:
            E_pts  = _v_ec_pts(vf, vg, vp)
            E2_pts = _v_ec_pts(f2, g2, vp)
            ker_x  = {q[0] % vp for q in S}

            # choose a source point P (not in the kernel)
            xs_av = sorted({x for (x, y) in E_pts if x not in ker_x})
            with v_ctrl:
                st.markdown("**Point $P$ on $E$**")
                if xs_av:
                    px = st.select_slider("x (P)",
                                          options=[_v_sym(x, vp) for x in xs_av],
                                          key="velu_px")
                    pxr = px % vp
                    py_opts = sorted({y for (x, y) in E_pts if x == pxr})
                    py = (st.radio("y (P)", [_v_sym(y, vp) for y in py_opts],
                                   horizontal=True, key="velu_py")
                          if len(py_opts) > 1 else _v_sym(py_opts[0], vp))
                    pyr = py % vp
                    PX, PY = phi(pxr, pyr)
                    st.markdown(
                        f"$\\varphi(P) = ({_v_sym(PX, vp)},\\ {_v_sym(PY, vp)})$"
                    )
                else:
                    px = py = None

            st.caption(
                f"$\\#E(\\mathbb{{F}}_{{{vp}}}) = {len(E_pts)+1}"
                f" = \\#(E/C)(\\mathbb{{F}}_{{{vp}}}) = {len(E2_pts)+1}$ "
                "— isogenous curves have equal point counts."
            )

            h = vp // 2
            amb = [(_v_sym(x, vp), _v_sym(y, vp))
                   for x in range(vp) for y in range(vp)]

            fig_v, (axE, axE2) = plt.subplots(1, 2, figsize=(9, 4.6))

            for ax, pts, ttl in [
                (axE,  E_pts,  rf"$E:\ y^2=x^3+{vf%vp}x+{vg%vp}$"),
                (axE2, E2_pts, rf"$E/C:\ y^2=x^3+{f2}x+{g2}$"),
            ]:
                ax.scatter([q[0] for q in amb], [q[1] for q in amb],
                           color="gray", alpha=0.18, s=6, zorder=1)
                ax.scatter([_v_sym(x, vp) for (x, y) in pts],
                           [_v_sym(y, vp) for (x, y) in pts],
                           color="steelblue", s=22, zorder=3)
                ax.set_xlim(-h - 0.5, h + 0.5)
                ax.set_ylim(-h - 0.5, h + 0.5)
                ax.set_aspect("equal"); ax.set_frame_on(False)
                ax.set_xticks([]); ax.set_yticks([])
                ax.set_title(ttl, fontsize=10)

            # kernel points (red) on E
            for (xq, yq) in S:
                for sy in ({0} if vell == 2 else {yq, (-yq) % vp}):
                    axE.scatter(_v_sym(xq, vp), _v_sym(sy, vp),
                                color="crimson", s=70, zorder=5,
                                edgecolor="white", linewidth=0.7)
            # P on E and its image on E/C
            if px is not None:
                axE.scatter(px, py, color="seagreen", s=90, zorder=6)
                axE.annotate("$P$", (px, py), xytext=(5, 5),
                             textcoords="offset points", color="seagreen",
                             fontsize=12, fontweight="bold")
                axE2.scatter(_v_sym(PX, vp), _v_sym(PY, vp),
                             color="seagreen", s=90, zorder=6)
                axE2.annotate("$\\varphi(P)$", (_v_sym(PX, vp), _v_sym(PY, vp)),
                              xytext=(5, 5), textcoords="offset points",
                              color="seagreen", fontsize=12, fontweight="bold")

            fig_v.tight_layout()
            st.pyplot(fig_v)
            plt.close(fig_v)
            st.caption(
                "Red: the kernel $C$ on $E$. Green: a point $P$ and its image "
                "$\\varphi(P)$ on the codomain $E/C$. Every point of $E$ maps to a "
                "point of $E/C$, exactly $\\ell$-to-$1$."
            )


# ══ §5 — Modular tools ════════════════════════════════════════════════════════

# ── Tab: the j-function and modular curves ────────────────────────────────────
with tab_jfun:
    st.subheader("The $j$-function and Modular Curves")

    st.markdown(
        "So far we have moved between a curve and its lattice by hand. To do this "
        "*in bulk* — and to read off isogenies without ever writing down a kernel "
        "— we need the classical machinery of **modular curves**."
    )
    st.markdown(
        "The **$j$-invariant** is a single number $j(E) \\in \\overline{\\mathbb{F}}_p$ "
        "(or $\\mathbb{C}$) that determines an elliptic curve up to isomorphism over "
        "the algebraic closure. Over $\\mathbb{C}$, writing $E = \\mathbb{C}/\\langle "
        "1, \\tau\\rangle$, the value $j(\\tau)$ depends only on the lattice, and "
        "$j$ identifies the quotient"
    )
    st.latex(r"X(1) \;=\; \mathrm{SL}_2(\mathbb{Z}) \backslash \mathcal{H}"
             r"\;\xrightarrow{\;\sim\;}\; \mathbb{C}, \qquad \tau \longmapsto j(\tau).")
    st.markdown(
        "So a point of the **modular curve $X(1)$** *is* an isomorphism class of "
        "elliptic curves. This is the space our whole gallery lives in."
    )
    st.markdown(
        "To talk about isogenies we go one level up. The modular curve "
        "**$X_0(\\ell)$** parametrizes pairs $(E, C)$ where $C \\subset E$ is a "
        "cyclic subgroup of order $\\ell$ — equivalently, a cyclic $\\ell$-isogeny "
        "$E \\to E/C$. It comes with two maps down to $X(1)$,"
    )
    st.latex(r"(E, C) \;\longmapsto\; j(E) \quad\text{and}\quad (E, C) \;\longmapsto\; j(E/C),")
    st.markdown(
        "the *source* and *target* $j$-invariants. The image of the combined map "
        "$X_0(\\ell) \\to X(1) \\times X(1)$ is a plane curve — and its defining "
        "equation is exactly the **modular polynomial** of the next tab. For small "
        "$\\ell$ these curves have genus $0$ (the *Atkin primes* "
        "$\\ell \\in \\{2, \\dots, 71\\}$ this project leans on), which is what makes "
        "them cheap to compute with."
    )

# ── Tab: modular polynomials Φ_ℓ ──────────────────────────────────────────────
with tab_modpoly:
    st.subheader("Modular Polynomials $\\Phi_\\ell$")

    st.markdown(
        "The **modular polynomial** $\\Phi_\\ell(X, Y) \\in \\mathbb{Z}[X, Y]$ is a "
        "symmetric polynomial with integer coefficients characterized by a single "
        "property:"
    )
    st.latex(r"\Phi_\ell\big(j(E),\, j(E')\big) = 0 \quad\Longleftrightarrow\quad"
             r"\text{there is a cyclic } \ell\text{-isogeny } E \to E'.")
    st.markdown(
        "Because the defining property is about $j$-invariants alone, the *same* "
        "polynomial works over $\\mathbb{C}$ and over every $\\mathbb{F}_p$ (reduce "
        "the coefficients mod $p$). This is the tool that turns isogenies into "
        "**root-finding**:"
    )
    st.markdown(
        "- The neighbours of a curve $E$ in the **$\\ell$-isogeny graph** are exactly "
        "the roots $Y$ of the one-variable polynomial $\\Phi_\\ell\\big(j(E), Y\\big) "
        "\\bmod p$.\n"
        "- So we can walk the isogeny graph **without ever computing a kernel or a "
        "quotient** — just factor $\\Phi_\\ell(j, Y)$ over $\\mathbb{F}_p$. (Vélu's "
        "formulas of §4 are the alternative, and remain a fallback in the code.)"
    )
    st.markdown(
        "The catch is size: $\\Phi_\\ell$ has degree $\\ell + 1$ in each variable and "
        "coefficients that grow very quickly, so only finitely many are practical to "
        "store. `ecfplat` ships the genus-$0$ Atkin polynomials $\\ell \\in \\{2, "
        "\\dots, 71\\}$ and extends the pool with classical $\\Phi_\\ell$ computed on "
        "demand (via a CRT / Hilbert-class-polynomial construction) whenever a "
        "discriminant needs a prime the current pool cannot provide — the "
        "**bootstrapping loop** behind the $97.04\\% \\to 99.86\\%$ jump on the home "
        "page."
    )


# ══ §6 — The CM bijection ═════════════════════════════════════════════════════

# ── Tab: the gallery problem and the equivalence ──────────────────────────────
with tab_setup:
    st.subheader("The Gallery Problem and the Equivalence")

    st.markdown(
        "The lattice pictures of §3 look great and they faithfully capture the "
        "arithmetic, so of course we want a whole **gallery**: one labelled picture "
        "for every curve in an isogeny class. Here is the obstruction."
    )
    st.markdown(
        "Over $\\mathbb{F}_p$ there are always **many curves with the same trace of "
        "Frobenius $a$ and the same endomorphism ring**. We can draw all their "
        "pictures immediately — one per lattice class of discriminant $a^2 - 4p$ — "
        "and we can label *one* of them freely: as long as the endomorphism rings "
        "match, some prime ideal $\\mathfrak{P}$ makes that one labelling correct. "
        "**After that first choice we are stuck.** Labelling the rest of the gallery "
        "requires real mathematics — and that is what §6 is about."
    )

    st.divider()
    st.markdown("#### Setup: two categories")
    st.markdown(
        "Fix $\\chi(x) = x^2 - a x + p$ with $a^2 - 4p < 0$, and let $\\alpha$ be a "
        "root — a well-defined algebraic integer once the class is **ordinary** "
        "($a \\neq 0$). We compare two categories:"
    )
    st.markdown(
        "- $\\mathcal{A}(\\alpha)$ — the category of **analytic curves** "
        "$E_\\Lambda = \\mathbb{C}/\\Lambda$ equipped with a root "
        "$[\\alpha] \\in \\mathrm{End}(E_\\Lambda)$ of $\\chi$.\n"
        "- $\\mathcal{E}_{\\mathbb{F}_p}(a)$ — the category of **elliptic curves over "
        "$\\mathbb{F}_p$** with trace of Frobenius $a$; morphisms are isogenies "
        "defined over $\\mathbb{F}_p$. (The “defined over $\\mathbb{F}_p$” "
        "condition only bites when $a = 0$.)"
    )

    st.divider()
    st.markdown("#### Lifting is an equivalence of categories")
    st.markdown(
        "Every curve downstairs has a **unique** lift of Frobenius "
        "$(E_\\Lambda, [\\alpha])$ upstairs, and every object upstairs is the lift of "
        "a unique curve downstairs. Choosing a prime ideal $\\mathfrak{P}$ turns this "
        "into an **equivalence of categories**"
    )
    st.latex(r"\mathcal{F}_\mathfrak{P} : \mathcal{A}(\alpha) \;\xrightarrow{\;\sim\;}\; \mathcal{E}_{\mathbb{F}_p}(a),")
    st.markdown(
        "and conversely **every** equivalence $\\mathcal{A}(\\alpha) \\to "
        "\\mathcal{E}_{\\mathbb{F}_p}(a)$ arises as $\\mathcal{F}_\\mathfrak{P}$ for "
        "some $\\mathfrak{P}$. The upshot reframes the whole problem:"
    )
    st.info(
        "**Labelling the gallery $\\iff$ computing any equivalence of categories** — "
        "and both categories have only finitely many objects."
    )
    st.markdown(
        "**An aside worth the price of admission.** The category "
        "$\\mathcal{E}_{\\mathbb{F}_p}(a)$ depends only on the discriminant "
        "$a^2 - 4p$. So if $a_1^2 - 4p_1 = a_2^2 - 4p_2$ then"
    )
    st.latex(r"\mathcal{E}_{\mathbb{F}_{p_1}}(a_1) \;\simeq\; \mathcal{E}_{\mathbb{F}_{p_2}}(a_2)")
    st.markdown(
        "— an equivalence between elliptic curves in **different characteristics**, "
        "invisible to methods that work one prime at a time."
    )

# ── Tab: isogeny graphs and rigidity ──────────────────────────────────────────
with tab_rigid:
    st.subheader("Isogeny Graphs and Rigidity")

    st.markdown(
        "How do we pin the equivalence down? Through its shadow on **isogeny "
        "graphs**. Write $X(\\mathcal{C})$ for the set of isomorphism classes in a "
        "category $\\mathcal{C}$; the $\\ell$-isogeny graph puts one edge "
        "$[E_1] \\to [E_2]$ per degree-$\\ell$ morphism. An equivalence of categories "
        "induces **isomorphisms of $\\ell$-isogeny graphs for every $\\ell$** at "
        "once."
    )
    st.markdown(
        "Below is the $2$-isogeny graph for $\\mathrm{disc} = -368$: a **volcano** "
        "(a crater rim of curves with the maximal endomorphism ring, plus trees "
        "descending to the sub-orders). The point of the equivalence is that this is "
        "the **same graph** for the lattice classes and for *every* $(a, p)$ with "
        "$a^2 - 4p = -368$. Use the pills to relabel the identical graph — abstract "
        "lattice classes, or the $j$-invariants mod $p$ for a series of primes:"
    )
    components.html(slide_viz.volcano_html(), height=650, scrolling=False)

    st.divider()
    st.markdown("#### Rigidity: a finite certificate")
    st.markdown(
        "The converse is what makes the labelling *computable*. For each $\\chi$ "
        "there is a **finite** set $S = \\{\\ell_1, \\dots, \\ell_n\\}$ — a **rigid "
        "spanning set** — with the property that any bijection "
        "$X(\\mathcal{A}(\\alpha)) \\to X(\\mathcal{E}_{\\mathbb{F}_p}(a))$ that "
        "induces isomorphisms of $\\ell$-isogeny graphs **for all $\\ell \\in S$** "
        "already arises from an equivalence of categories. Checking finitely many "
        "graphs certifies the whole labelling."
    )
    st.markdown(
        "An explicit small $S$ can be read off the class group "
        "$\\mathrm{Cl}(\\mathbb{Z}[\\alpha])$:"
    )
    st.markdown(
        "- one prime $\\ell_\\mathfrak{b}$ for each generator $(\\mathfrak{b})$ in a "
        "**minimal generating set** $B$ of $\\mathrm{Cl}(\\mathbb{Z}[\\alpha])$ "
        "(spanning $\\iff$ generating);\n"
        "- one **orientation prime** $\\ell^{*}$ when two or more generators have "
        "order $> 2$;\n"
        "- the prime factors of $\\mathrm{cond}(\\mathbb{Z}[\\alpha])$, to handle the "
        "vertical extension up and down the volcano."
    )

# ── Tab: the algorithm and the gallery ────────────────────────────────────────
with tab_algo:
    st.subheader("The Algorithm and the Gallery")

    st.markdown(
        "The equivalences form a **single orbit** under "
        "$\\mathrm{Cl}(\\mathcal{O}) \\times \\{\\pm 1\\}$ — there are far more "
        "curves to label than there are degrees of freedom in the labelling. That "
        "gap is exactly what makes the gallery rigid, and it splits the algorithm "
        "into **a few free choices** followed by **forced propagation**."
    )
    st.markdown("#### The initial values (the only free choices)")
    st.markdown(
        "- Pick $\\mathcal{F}((\\mathcal{O}))$ — where the **trivial class** goes. "
        "This is essentially the choice of prime ideal $\\mathfrak{P}$.\n"
        "- Pick an **orientation**: $\\mathcal{F}((\\mathfrak{b}^{*}))$ among the two "
        "neighbours of $\\mathcal{F}((\\mathcal{O}))$ in the $\\ell^{*}$-graph "
        "(this is the $\\mathfrak{P} \\leftrightarrow \\overline{\\mathfrak{P}}$ "
        "ambiguity)."
    )
    st.markdown("#### Propagation (everything else is forced)")
    st.markdown(
        "Read the rest straight off the graphs for $\\ell \\in S$, since these "
        "graphs are **Cayley graphs of $\\mathrm{Cl}(\\mathcal{O})$**:"
    )
    st.markdown(
        "- **radical words** via the isogeny trees $\\to$ **oriented cycle walks** "
        "around the crater;\n"
        "- **order-$2$ generators** are unambiguous (no orientation needed);\n"
        "- **vertical** steps: match $\\ell$-ancestors up the volcano through "
        "$\\mathrm{cond}(\\mathbb{Z}[\\alpha])$."
    )
    st.markdown(
        "The output is the **unique** bijection extending the chosen initial values "
        "that arises from an equivalence: a complete, correctly labelled gallery."
    )

    st.divider()
    st.markdown("#### The gallery — and all of its labellings")
    st.markdown(
        "Here is the payoff for $(a, p) = (6, 101)$, $\\mathrm{disc} = -368$: the "
        "same graph carrying the $2$-isogeny edges (solid) and $3$-isogeny cycles "
        "(dashed), with a lattice picture beside each floor class. The buttons walk "
        "the orbit of valid labellings:"
    )
    st.markdown(
        "- **rotate** = act by a class of norm $3$; every label steps along its "
        "$3$-isogeny cycle.\n"
        "- **reflect** = complex conjugation, "
        "$\\mathfrak{P} \\leftrightarrow \\overline{\\mathfrak{P}}$."
    )
    components.html(slide_viz.gallery_html(), height=580, scrolling=False)
    st.caption(
        "The cards show the six floor classes with their lattice pictures; hover a "
        "card for its Weierstrass equation. All "
        "$12 = |\\mathrm{Cl}(\\mathbb{Z}[\\alpha])| \\cdot |\\{\\pm 1\\}|$ labellings "
        "are equally correct, and none is canonical. This computation has been "
        "carried out for 117,155 ordinary classes with $4 \\leq p \\leq 8192$ "
        "(99.86%), plus every supersingular class in that range."
    )
