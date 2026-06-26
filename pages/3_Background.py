import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.patches import Arc as MplArc


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

with st.expander("Miscellaneous topics", expanded=False):
    (tab_volc,) = st.tabs([
        "Isogenies over $\\mathbb{F}_p$: volcanoes",
    ])



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
    ex_col, ex_plot = st.columns([1, 2])
    with ex_col:
        f_ex    = st.slider("f", -5.0, 5.0, -1.0, 0.1, key="bg1_f_ex")
        g_ex    = st.slider("g", -5.0, 5.0,  1.0, 0.1, key="bg1_g_ex")
        disc_ex = -16 * (4*f_ex**3 + 27*g_ex**2)
        if abs(disc_ex) < 1e-4:
            st.warning("Singular curve (Δ ≈ 0)")
        elif disc_ex > 0:
            st.info("Two components (Δ > 0)")
        else:
            st.info("One component (Δ < 0)")

    with ex_plot:
        xs_ex = np.linspace(-3.3, 3.3, 3000)
        y2_ex = xs_ex**3 + f_ex*xs_ex + g_ex
        yp_ex = np.where(y2_ex >= 0, np.sqrt(np.clip(y2_ex, 0, None)), np.nan)
        fig_ex2, ax_ex2 = plt.subplots(figsize=(5, 4))
        ax_ex2.plot(xs_ex,  yp_ex, color="steelblue", lw=2)
        ax_ex2.plot(xs_ex, -yp_ex, color="steelblue", lw=2)
        ax_ex2.axhline(0, color="k", lw=0.4)
        ax_ex2.axvline(0, color="k", lw=0.4)
        ax_ex2.set_xlim(-3.3, 3.3)
        ax_ex2.set_ylim(-4.5, 4.5)
        ax_ex2.set_title(f"$y^2 = x^3 {f_ex:+.1f}x {g_ex:+.1f}$", fontsize=11)
        ax_ex2.set_frame_on(False)
        st.pyplot(fig_ex2)
        plt.close(fig_ex2)

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
        "Use the applet below to explore the group law on a torus. "
        "Points $z_1$ and $z_2$ are specified by their fractional coordinates "
        "$(s, t)$ in the parallelogram: $z = s \\cdot 1 + t \\cdot \\tau$, with $s, t \\in [0, 1)$. "
        "When their sum falls outside the fundamental domain, it is translated back in — "
        "this is the 'wrapping around' of the torus."
    )

    st.divider()

    # ── Group law applet: click two points, compute their sum ─────────────────
    ctrl3, plot3 = st.columns([1, 2])

    with ctrl3:
        st.markdown("**Lattice**")
        tau_re3 = st.slider("Re(τ)", -0.5, 0.5, 0.2, 0.01, key="bg3c_re")
        tau_im3 = st.slider("Im(τ)",  0.1, 3.0, 1.2, 0.05, key="bg3c_im")
        _sign3  = "+" if tau_re3 >= 0 else "-"
        st.latex(rf"\tau = {tau_re3:.2f} {_sign3} {abs(tau_im3):.2f}\,i")

    # Clickable grid of fractional coordinates (s, t) ∈ [0,1)².
    _STEP3 = 0.05
    _grid3 = [(round(s, 3), round(t, 3))
              for s in np.arange(0.0, 1.0, _STEP3)
              for t in np.arange(0.0, 1.0, _STEP3)]
    _grid3_set = set(_grid3)

    sel3 = [pt for pt in st.session_state.get("bg3_sel", []) if pt in _grid3_set]
    st.session_state["bg3_sel"] = sel3
    show3 = (st.session_state.get("bg3_mode") == "sum") and len(sel3) == 2

    s3 = t3 = s3_raw = t3_raw = None
    wrapped = False
    if show3:
        (s1, t1), (s2, t2) = sel3[0], sel3[1]
        s3_raw, t3_raw = s1 + s2, t1 + t2
        s3, t3 = s3_raw % 1.0, t3_raw % 1.0
        wrapped = (s3_raw >= 1.0) or (t3_raw >= 1.0)

    with ctrl3:
        st.markdown("**Pick two points**")
        st.caption("Click inside the parallelogram to choose $z_1$ then $z_2$ "
                   "(a third click clears). The identity $0$ is the olive corner.")
        for nm, pt in zip(("z₁", "z₂"), sel3):
            st.markdown(f"- **{nm}**: $(s, t) = ({pt[0]:.2f}, {pt[1]:.2f})$")
        c3a, c3b = st.columns(2)
        if c3a.button("Clear", key="bg3_clear", width="stretch"):
            st.session_state["bg3_sel"] = []
            st.session_state["bg3_mode"] = None
            st.rerun()
        if c3b.button("Compute sum", key="bg3_compute",
                      width="stretch", disabled=(len(sel3) != 2)):
            st.session_state["bg3_mode"] = "sum"
            st.rerun()
        if show3:
            st.success(f"z₁ + z₂:  (s, t) = ({s3:.2f}, {t3:.2f})")
            if wrapped:
                st.caption("(wrapped back into the fundamental domain)")

    with plot3:
        tau3 = np.array([tau_re3, tau_im3])
        one3 = np.array([1.0, 0.0])

        def _xy(s, t):
            return s * one3 + t * tau3

        verts3 = [np.zeros(2), one3, one3 + tau3, tau3]

        fig = go.Figure()
        # Adjacent copies (torus context)
        for dm in (-1, 0, 1):
            for dn in (-1, 0, 1):
                if dm == 0 and dn == 0:
                    continue
                shift = dm * one3 + dn * tau3
                vv = [v + shift for v in verts3] + [verts3[0] + shift]
                fig.add_trace(go.Scatter(
                    x=[p[0] for p in vv], y=[p[1] for p in vv],
                    mode="lines", fill="toself",
                    fillcolor="rgba(220,220,220,0.18)",
                    line=dict(color="lightgray", width=0.8),
                    hoverinfo="skip", showlegend=False))
        # Fundamental domain
        vv0 = verts3 + [verts3[0]]
        fig.add_trace(go.Scatter(
            x=[p[0] for p in vv0], y=[p[1] for p in vv0],
            mode="lines", fill="toself", fillcolor="rgba(178,178,210,0.45)",
            line=dict(color="steelblue", width=2),
            hoverinfo="skip", showlegend=False))

        # Sum as vector addition: the parallelogram 0–z₁–(z₁+z₂)–z₂ and the
        # two vectors out of the origin.
        if show3:
            z1xy = _xy(*sel3[0])
            z2xy = _xy(*sel3[1])
            raw_xy = _xy(s3_raw, t3_raw)
            para = [(0.0, 0.0), tuple(z1xy), tuple(raw_xy), tuple(z2xy), (0.0, 0.0)]
            fig.add_trace(go.Scatter(
                x=[p[0] for p in para], y=[p[1] for p in para], mode="lines",
                line=dict(color="gray", width=1, dash="dot"),
                hoverinfo="skip", showlegend=False))
            for tip, col in ((z1xy, "red"), (z2xy, "green")):
                fig.add_annotation(x=tip[0], y=tip[1], ax=0.0, ay=0.0,
                                   xref="x", yref="y", axref="x", ayref="y",
                                   showarrow=True, arrowcolor=col,
                                   arrowwidth=2, arrowhead=2)
            # Wrapping: unreduced point + arrow back into the domain
            if wrapped:
                red_xy = _xy(s3, t3)
                fig.add_trace(go.Scatter(
                    x=[raw_xy[0]], y=[raw_xy[1]], mode="markers",
                    marker=dict(color="orange", size=12, opacity=0.35),
                    hoverinfo="skip", showlegend=False))
                fig.add_annotation(x=red_xy[0], y=red_xy[1], ax=raw_xy[0], ay=raw_xy[1],
                                   xref="x", yref="y", axref="x", ayref="y",
                                   showarrow=True, arrowcolor="orange",
                                   arrowwidth=1.5, arrowhead=2)

        # "0" label at the identity (origin)
        fig.add_annotation(x=0.0, y=0.0, text="0", showarrow=False,
                           xshift=-8, yshift=-8, font=dict(size=13, color="olive"))

        # Clickable grid, with the origin, z1/z2 and sum coloured
        _role3  = {pt: nm for nm, pt in zip(("z1", "z2"), sel3)}
        sum_pt  = (round(s3, 3), round(t3, 3)) if show3 else None
        gcols, gsizes = [], []
        for pt in _grid3:
            if show3 and sum_pt is not None and pt == sum_pt:
                gcols.append("orange"); gsizes.append(15)
            elif pt in _role3:
                gcols.append("red" if _role3[pt] == "z1" else "green")
                gsizes.append(15)
            elif pt == (0.0, 0.0):
                gcols.append("olive"); gsizes.append(13)
            else:
                gcols.append("rgba(70,90,160,0.55)"); gsizes.append(6)
        gx = [_xy(s, t)[0] for (s, t) in _grid3]
        gy = [_xy(s, t)[1] for (s, t) in _grid3]
        fig.add_trace(go.Scatter(
            x=gx, y=gy, mode="markers",
            marker=dict(color=gcols, size=gsizes),
            customdata=list(range(len(_grid3))),
            hovertext=[f"(s, t) = ({s:.2f}, {t:.2f})" for (s, t) in _grid3],
            hoverinfo="text", showlegend=False,
            selected=dict(marker=dict(opacity=1.0)),
            unselected=dict(marker=dict(opacity=1.0))))

        fig.update_layout(
            margin=dict(l=10, r=10, t=38, b=10), height=520,
            plot_bgcolor="white", showlegend=False,
            title=dict(text="Group law on ℂ/Λ", x=0.5, xanchor="center",
                       font=dict(size=14)))
        fig.update_xaxes(visible=False)
        fig.update_yaxes(visible=False, scaleanchor="x", scaleratio=1)
        ev = st.plotly_chart(fig, width="stretch", on_select="rerun",
                             selection_mode="points", key="bg3_chart")

        idx = _click_index(ev)
        if idx is not None and 0 <= idx < len(_grid3):
            _two_point_pick("bg3", _grid3[idx])
            st.rerun()


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
        "Pick a prime $p$ and a curve $y^2 = x^3 + fx + g$. We count "
        "$\\#E(\\mathbb{F}_p)$ by brute force, read off $a = p + 1 - \\#E$, and "
        "place it inside the Hasse interval."
    )

    _FPRIMES = [n for n in range(5, 200)
                if all(n % dvr != 0 for dvr in range(2, int(n**0.5) + 1))]

    frob_ctrl, frob_plot = st.columns([1, 2])

    with frob_ctrl:
        fp = st.selectbox("p", _FPRIMES, index=_FPRIMES.index(23), key="frob_p")
        ff = int(st.number_input("f", value=1, step=1, key="frob_f"))
        fg = int(st.number_input("g", value=1, step=1, key="frob_g"))

        f_disc = (-16 * (4 * pow(ff, 3) + 27 * pow(fg, 2))) % fp
        if f_disc == 0:
            st.warning("Singular curve mod p — adjust f or g.")
            frob_ok = False
        else:
            frob_ok = True
            # Brute-force affine point count + point at infinity.
            n_affine = 0
            for x in range(fp):
                rhs = (pow(x, 3, fp) + ff * x + fg) % fp
                for y in range(fp):
                    if pow(y, 2, fp) == rhs:
                        n_affine += 1
            n_pts = n_affine + 1          # + O
            a_tr  = fp + 1 - n_pts
            d_fr  = a_tr * a_tr - 4 * fp
            st.success("Smooth curve ✓")
            st.markdown(
                f"$\\#E(\\mathbb{{F}}_{{{fp}}}) = {n_pts}$  \n"
                f"$a = p + 1 - \\#E = {a_tr}$  \n"
                f"$d = a^2 - 4p = {d_fr}$"
            )
            if a_tr % fp == 0:
                st.info("Supersingular ($a \\equiv 0 \\bmod p$).")
            else:
                st.info("Ordinary ($a \\not\\equiv 0 \\bmod p$).")

    with frob_plot:
        if frob_ok:
            bound = 2 * np.sqrt(fp)
            a_lo, a_hi = int(np.ceil(-bound)), int(np.floor(bound))
            fig_fr, ax_fr = plt.subplots(figsize=(7, 1.8))
            # all admissible traces
            ax_fr.scatter(range(a_lo, a_hi + 1), [0] * (a_hi - a_lo + 1),
                          color="lightgray", s=40, zorder=2)
            # Hasse bounds
            for s in (-bound, bound):
                ax_fr.axvline(s, color="steelblue", ls="--", lw=1.2, alpha=0.7)
            ax_fr.axvline(0, color="k", lw=0.4)
            # current trace
            ax_fr.scatter([a_tr], [0], color="red", s=130, zorder=4)
            ax_fr.annotate(f"$a = {a_tr}$", (a_tr, 0), xytext=(0, 12),
                           textcoords="offset points", ha="center",
                           color="red", fontsize=12, fontweight="bold")
            ax_fr.annotate(f"$-2\\sqrt{{p}} \\approx {-bound:.1f}$", (-bound, 0),
                           xytext=(0, -20), textcoords="offset points",
                           ha="center", color="steelblue", fontsize=9)
            ax_fr.annotate(f"$+2\\sqrt{{p}} \\approx {bound:.1f}$", (bound, 0),
                           xytext=(0, -20), textcoords="offset points",
                           ha="center", color="steelblue", fontsize=9)
            ax_fr.set_xlim(-bound - 2, bound + 2)
            ax_fr.set_ylim(-0.6, 0.6)
            ax_fr.set_yticks([])
            ax_fr.set_frame_on(False)
            ax_fr.set_title("Admissible traces in the Hasse interval",
                            fontsize=10)
            st.pyplot(fig_fr)
            plt.close(fig_fr)
            st.caption(
                "Every gray dot is an integer $a$ with $|a| \\le 2\\sqrt p$; each "
                "corresponds to an isogeny class over $\\mathbb{F}_p$. The red dot "
                "is the class of the curve above."
            )
        else:
            st.info("Adjust f and g to get a smooth curve.")


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
        "Choose $(a, p)$ with $a^2 < 4p$. We take the lattice "
        "$\\Lambda = \\mathbb{Z} + \\alpha\\,\\mathbb{Z}$ (the order $\\mathcal{O} = "
        "\\mathbb{Z}[\\alpha]$ itself) and watch multiplication by $\\alpha$ send its "
        "generators to other lattice points — the analytic shadow of Frobenius."
    )

    lift_ctrl, lift_plot = st.columns([1, 2])

    # Small primes keep the picture compact: multiplication by π scales by √p,
    # so π·π flies off-screen for large p.
    _LIFTPRIMES = [5, 7, 11, 13, 17, 19, 23]

    with lift_ctrl:
        lp = st.selectbox("p", _LIFTPRIMES, index=_LIFTPRIMES.index(7), key="lift_p")
        # admissible traces for this p (exclude a=0 so the lattice is non-real)
        a_max = int(np.floor(2 * np.sqrt(lp)))
        a_opts = [v for v in range(-a_max, a_max + 1) if v != 0]
        la = st.select_slider("a (trace of Frobenius)", options=a_opts,
                              value=1 if 1 in a_opts else a_opts[0], key="lift_a")
        d_lift = la * la - 4 * lp
        st.latex(rf"d = a^2 - 4p = {d_lift}")
        # alpha = a/2 + i sqrt(4p - a^2)/2
        pi_re = la / 2.0
        pi_im = np.sqrt(4 * lp - la * la) / 2.0
        st.latex(rf"\alpha = {pi_re:.2f} + {pi_im:.2f}\,i")
        st.latex(rf"|\alpha| = \sqrt{{{lp}}} \approx {np.sqrt(lp):.2f}")
        st.markdown(
            "On the basis $(1, \\alpha)$, multiplication by $\\alpha$ has the integer "
            "matrix (from $\\alpha\\cdot 1 = \\alpha$ and $\\alpha^2 = a\\alpha - p$):"
        )
        st.latex(
            rf"[\alpha] = \begin{{pmatrix}} 0 & -{lp} \\ 1 & {la} \end{{pmatrix}}"
        )
        st.caption(
            "Integer entries ⇒ $\\alpha\\Lambda \\subseteq \\Lambda$: Frobenius "
            "really is an endomorphism of the lattice."
        )

    with lift_plot:
        pi  = np.array([pi_re, pi_im])
        one = np.array([1.0, 0.0])

        # images of the two generators under mult-by-alpha:
        #   alpha·1     = alpha          -> (0, 1)  in (1, alpha) coords
        #   alpha·alpha = a·alpha - p    -> (-p, a) in (1, alpha) coords
        img_1  = pi                       # = 0·one + 1·alpha
        img_pi = la * pi - lp * one       # = -p·one + a·alpha

        # Window must contain every vertex of both fundamental cells:
        #   blue cell   : 0, 1, 1+alpha, alpha
        #   orange cell : 0, alpha, alpha+(alpha·alpha), alpha·alpha
        key_pts = np.array([
            np.zeros(2), one, pi, one + pi,           # blue cell vertices
            img_1, img_pi, img_1 + img_pi,            # orange cell vertices
        ])
        pad = 1.0
        xlo, xhi = key_pts[:, 0].min() - pad, key_pts[:, 0].max() + pad
        ylo, yhi = key_pts[:, 1].min() - pad, key_pts[:, 1].max() + pad

        # Fill the window with the actual lattice m·1 + n·pi.
        lat = []
        n_lo = int(np.floor(ylo / pi_im)) - 1
        n_hi = int(np.ceil(yhi / pi_im)) + 1
        for n in range(n_lo, n_hi + 1):
            base_x = n * pi_re
            m_lo = int(np.floor(xlo - base_x)) - 1
            m_hi = int(np.ceil(xhi - base_x)) + 1
            for m in range(m_lo, m_hi + 1):
                lat.append((m + base_x, n * pi_im))

        fig_lf, ax_lf = plt.subplots(figsize=(6, 6))
        ax_lf.scatter([q[0] for q in lat], [q[1] for q in lat],
                      color="gray", alpha=0.35, s=14, zorder=1)

        # original fundamental cell (spanned by 1 and alpha) ...
        ax_lf.add_patch(MplPolygon(
            [np.zeros(2), one, one + pi, pi],
            facecolor=[0.85, 0.85, 0.95, 0.5],
            edgecolor="steelblue", lw=1.8, zorder=2))
        # ... and its image under mult-by-alpha (spanned by alpha and alpha^2)
        ax_lf.add_patch(MplPolygon(
            [np.zeros(2), img_1, img_1 + img_pi, img_pi],
            facecolor=[0.95, 0.85, 0.70, 0.30],
            edgecolor="darkorange", lw=1.5, ls="--", zorder=2))

        # generator arrows
        for vec, lbl, col in [(one, "$1$", "black"), (pi, "$\\alpha$", "red")]:
            ax_lf.annotate("", xy=vec, xytext=(0, 0),
                           arrowprops=dict(arrowstyle="->", color=col, lw=2),
                           zorder=4)
            ax_lf.annotate(lbl, vec, xytext=(6, 6),
                           textcoords="offset points",
                           color=col, fontsize=13, fontweight="bold")

        # image of alpha^2 (far generator image); alpha·1 = alpha is the red arrow
        ax_lf.annotate("", xy=img_pi, xytext=(0, 0),
                       arrowprops=dict(arrowstyle="->", color="darkorange",
                                       lw=1.6, linestyle="dashed", alpha=0.9),
                       zorder=3)
        ax_lf.scatter(*img_pi, color="darkorange", s=80, zorder=5)
        ax_lf.annotate("$\\alpha\\cdot\\alpha = a\\alpha - p$", img_pi,
                       xytext=(6, 6), textcoords="offset points",
                       color="darkorange", fontsize=10)
        ax_lf.scatter(*img_1, color="red", s=60, zorder=5)
        ax_lf.annotate("$\\alpha\\cdot 1 = \\alpha$", img_1, xytext=(6, -14),
                       textcoords="offset points", color="red", fontsize=9)

        ax_lf.set_xlim(xlo, xhi)
        ax_lf.set_ylim(ylo, yhi)
        ax_lf.set_aspect("equal")
        ax_lf.axhline(0, color="k", lw=0.4)
        ax_lf.axvline(0, color="k", lw=0.4)
        ax_lf.set_frame_on(False)
        ax_lf.set_title("$\\Lambda = \\mathbb{Z} + \\alpha\\mathbb{Z}$ "
                        "and multiplication by $\\alpha$", fontsize=11)
        st.pyplot(fig_lf)
        plt.close(fig_lf)
        st.caption(
            "Blue cell: the fundamental domain spanned by $1, \\alpha$. Orange "
            "cell: its image under multiplication by $\\alpha$ — a rotation by "
            "$\\arg\\alpha$ and scaling by $|\\alpha| = \\sqrt p$. Both image "
            "generators $\\alpha\\cdot 1$ and $\\alpha\\cdot\\alpha$ land on lattice "
            "points, so $\\alpha\\Lambda \\subseteq \\Lambda$."
        )


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
        "$C$ pins down the codomain $E/C$. Pick a lattice $\\Lambda$ and a prime "
        "degree $d$; the $d \\times d$ grid of **$d$-torsion points** $E[d]$ is "
        "shown in the fundamental domain on the left. **Click a nonzero point** to "
        "choose a generator $P$ of an order-$d$ kernel $C = \\langle P\\rangle$ — "
        "the codomain $E/C = \\mathbb{C}/\\Lambda'$ then appears on the right, with "
        "$\\Lambda' = \\Lambda + \\mathbb{Z}P$."
    )

    ik_ctrl, ik_left, ik_right = st.columns([1, 2, 2])

    with ik_ctrl:
        ik_re = st.slider("Re(τ)", -0.5, 0.5, 0.30, 0.01, key="ik_re")
        ik_im = st.slider("Im(τ)", 0.4, 2.0, 1.00, 0.05, key="ik_im")
        ik_d = st.select_slider("degree d (prime)", options=[2, 3, 5, 7],
                                value=3, key="ik_d")

    d = int(ik_d)
    u = np.array([1.0, 0.0])
    v = np.array([float(ik_re), float(ik_im)])

    def _pos(cu, cv):                       # (1, τ)-coordinates → (x, y)
        return cu * u + cv * v

    grid = [(a, b) for a in range(d) for b in range(d)]   # index ↔ this order
    kkey = f"ik_ker_{d}"
    ker = st.session_state.get(kkey)        # chosen generator (a, b) or None
    if ker is not None and (ker[0] >= d or ker[1] >= d):
        ker = None

    # kernel subgroup ⟨P⟩ as a set of torsion coordinates
    ker_set = set()
    if ker is not None:
        ker_set = {((ker[0] * k) % d, (ker[1] * k) % d) for k in range(d)}

    # shared square view box (from the fundamental-domain corners)
    corners = [_pos(0, 0), _pos(1, 0), _pos(1, 1), _pos(0, 1)]
    cxs = [c[0] for c in corners]; cys = [c[1] for c in corners]
    cx, cy = (min(cxs) + max(cxs)) / 2, (min(cys) + max(cys)) / 2
    hs = max(max(cxs) - min(cxs), max(cys) - min(cys)) / 2 + 0.2
    xrng, yrng = [cx - hs, cx + hs], [cy - hs, cy + hs]

    def _fd_path():
        c = corners
        return (f"M {c[0][0]},{c[0][1]} L {c[1][0]},{c[1][1]} "
                f"L {c[2][0]},{c[2][1]} L {c[3][0]},{c[3][1]} Z")

    def _layout(fig, title):
        fig.add_shape(type="path", path=_fd_path(),
                      fillcolor="rgba(140,140,140,0.10)",
                      line=dict(color="rgba(130,130,130,0.6)", width=1),
                      layer="below")
        fig.update_layout(
            title=dict(text=title, font=dict(size=13), x=0.5, xanchor="center"),
            xaxis=dict(visible=False, range=xrng),
            yaxis=dict(visible=False, range=yrng, scaleanchor="x", scaleratio=1),
            margin=dict(l=0, r=0, t=34, b=0), height=400, showlegend=False,
            plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)",
        )

    # ── Left figure: domain with clickable d-torsion grid ─────────────────────
    gx, gy, colors, sizes, symbols = [], [], [], [], []
    for (a, b) in grid:
        xy = _pos(a / d, b / d)
        gx.append(xy[0]); gy.append(xy[1])
        if (a, b) == (0, 0):
            colors.append("#888"); sizes.append(10); symbols.append("circle")
        elif (a, b) == ker:
            colors.append("crimson"); sizes.append(16); symbols.append("star")
        elif (a, b) in ker_set:
            colors.append("crimson"); sizes.append(11); symbols.append("circle")
        else:
            colors.append("steelblue"); sizes.append(10); symbols.append("circle")

    figL = go.Figure()
    figL.add_trace(go.Scatter(
        x=gx, y=gy, mode="markers",
        marker=dict(color=colors, size=sizes, symbol=symbols,
                    line=dict(width=0.6, color="white")),
        customdata=grid,
        hovertemplate="P = (%{customdata[0]}, %{customdata[1]}) / "
                      + str(d) + "<extra></extra>",
    ))
    _layout(figL, "Domain  E = ℂ/Λ   (grid = E[d])")

    with ik_left:
        ev = st.plotly_chart(figL, on_select="rerun", selection_mode="points",
                             key=f"ik_left_{d}", config={"displayModeBar": False})

    # read the click → new generator
    new_ab = None
    try:
        pts = ev["selection"]["points"]
    except (TypeError, KeyError):
        pts = None
    if pts:
        pt = pts[-1]
        idx = pt.get("point_index", pt.get("point_number"))
        if idx is not None and 0 <= idx < len(grid):
            cand = grid[idx]
            if cand != (0, 0):
                new_ab = cand
    if new_ab is not None and new_ab != ker:
        st.session_state[kkey] = new_ab
        st.rerun()

    # ── Right figure: codomain lattice Λ' = Λ + ZP ────────────────────────────
    with ik_right:
        if ker is None:
            st.info("Click a nonzero point on the left to choose a kernel "
                    "generator $P$.")
        else:
            ox, oy, nx, ny = [], [], [], []
            for m in range(-1, 3):
                for n in range(-1, 3):
                    for k in range(d):
                        xy = _pos(m + k * ker[0] / d, n + k * ker[1] / d)
                        if k == 0:
                            ox.append(xy[0]); oy.append(xy[1])
                        else:
                            nx.append(xy[0]); ny.append(xy[1])
            figR = go.Figure()
            figR.add_trace(go.Scatter(
                x=nx, y=ny, mode="markers", name="added by C",
                marker=dict(color="goldenrod", size=9,
                            line=dict(width=0.5, color="white")),
                hoverinfo="skip"))
            figR.add_trace(go.Scatter(
                x=ox, y=oy, mode="markers", name="Λ",
                marker=dict(color="#888", size=9,
                            line=dict(width=0.5, color="white")),
                hoverinfo="skip"))
            _layout(figR,
                    f"Codomain  E/C = ℂ/Λ′   (P = ({ker[0]},{ker[1]})/{d})")
            st.plotly_chart(figR, key=f"ik_right_{d}",
                            config={"displayModeBar": False})

    st.caption(
        "Left: the $d$-torsion $E[d]$ in a fundamental domain of $\\Lambda$; the "
        "chosen generator $P$ is the star and its subgroup $\\langle P\\rangle$ is "
        "red. Right: the codomain lattice $\\Lambda' = \\Lambda + \\mathbb{Z}P$ — "
        "the original lattice $\\Lambda$ (gray) together with the cosets the kernel "
        "adds (gold), inside the same fundamental domain. A different kernel gives "
        "a different $\\Lambda'$."
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
        "Pick a lattice $\\Lambda = \\mathbb{Z} + \\tau\\mathbb{Z}$, a degree "
        "$\\ell$, and a cyclic kernel direction. The point $P$ (and its $\\ell$ "
        "kernel-translates, the whole fibre) is shown on the left; its single image "
        "$\\varphi(P)$ is shown on the right."
    )

    st.divider()

    fold_ctrl, fold_plot = st.columns([1, 2])

    with fold_ctrl:
        st.markdown("**Lattice**")
        ftau_re = st.slider("Re(τ)", -0.5, 0.5, 0.25, 0.01, key="fold_re")
        ftau_im = st.slider("Im(τ)",  0.4, 2.0, 1.0, 0.05, key="fold_im")
        st.markdown("**Isogeny**")
        fell = st.select_slider("degree ℓ", options=[2, 3, 4, 5], value=3,
                                key="fold_ell")
        fdir = st.radio(
            "cyclic kernel",
            ["⟨1/ℓ⟩  (vertical slabs)", "⟨τ/ℓ⟩  (horizontal slabs)"],
            key="fold_dir")
        vertical = fdir.startswith("⟨1/ℓ⟩")
        st.markdown("**Point $P = s\\cdot 1 + t\\cdot\\tau$**")
        fps = st.slider("s", 0.0, 1.0, 0.30, 0.02, key="fold_s")
        fpt = st.slider("t", 0.0, 1.0, 0.55, 0.02, key="fold_t")
        st.caption(
            f"$\\Lambda' = "
            + (r"\tfrac{1}{%d}\mathbb{Z} + \tau\mathbb{Z}$" % fell if vertical
               else r"\mathbb{Z} + \tfrac{1}{%d}\tau\mathbb{Z}$" % fell)
            + f"  — kernel of order {fell}."
        )

    with fold_plot:
        one_f = np.array([1.0, 0.0])
        tau_f = np.array([ftau_re, ftau_im])

        def _xy(s, t):
            return s * one_f + t * tau_f

        # kernel generator fractions and the target cell basis
        if vertical:
            ker_pts = [_xy(k / fell, 0.0) for k in range(fell)]
            tgt_u, tgt_v = _xy(1.0 / fell, 0.0), tau_f      # small cell basis
            s_img = (fps % (1.0 / fell)); t_img = fpt
        else:
            ker_pts = [_xy(0.0, k / fell) for k in range(fell)]
            tgt_u, tgt_v = one_f, _xy(0.0, 1.0 / fell)
            s_img = fps; t_img = (fpt % (1.0 / fell))

        # fibre of P = the ell kernel-translates
        if vertical:
            fibre = [_xy((fps + k / fell) % 1.0, fpt) for k in range(fell)]
        else:
            fibre = [_xy(fps, (fpt + k / fell) % 1.0) for k in range(fell)]
        img_pt = s_img * one_f + t_img * tau_f

        strip_cols = plt.cm.viridis(np.linspace(0.15, 0.85, fell))

        fig_fd, (axL, axR) = plt.subplots(1, 2, figsize=(9, 4.2))

        # ── Left: E = C/Lambda, sliced into ell slabs ─────────────────────────
        for k in range(fell):
            if vertical:
                corners = [_xy(k / fell, 0), _xy((k + 1) / fell, 0),
                           _xy((k + 1) / fell, 1), _xy(k / fell, 1)]
            else:
                corners = [_xy(0, k / fell), _xy(1, k / fell),
                           _xy(1, (k + 1) / fell), _xy(0, (k + 1) / fell)]
            axL.add_patch(MplPolygon(corners, facecolor=strip_cols[k],
                                     alpha=0.35, edgecolor="gray", lw=0.7,
                                     zorder=1))
        # kernel points
        for K in ker_pts:
            axL.scatter(*K, color="crimson", s=55, zorder=4,
                        edgecolor="white", linewidth=0.7)
        # fibre of P
        for Q in fibre:
            axL.scatter(*Q, color="black", s=50, zorder=5)
        axL.scatter(*fibre[0], color="black", s=90, zorder=6)
        axL.annotate("$P$", fibre[0], xytext=(6, 6),
                     textcoords="offset points", fontsize=12, fontweight="bold")
        axL.set_title("$E = \\mathbb{C}/\\Lambda$  (sliced into $\\ell$)",
                      fontsize=10)

        # ── Right: E' = C/Lambda' ─────────────────────────────────────────────
        axR.add_patch(MplPolygon(
            [np.zeros(2), tgt_u, tgt_u + tgt_v, tgt_v],
            facecolor=[0.85, 0.9, 0.85, 0.5], edgecolor="seagreen",
            lw=1.8, zorder=1))
        axR.scatter(0, 0, color="crimson", s=55, zorder=4,
                    edgecolor="white", linewidth=0.7)
        axR.scatter(*img_pt, color="black", s=90, zorder=5)
        axR.annotate("$\\varphi(P)$", img_pt, xytext=(6, 6),
                     textcoords="offset points", fontsize=12, fontweight="bold")
        axR.set_title("$E' = \\mathbb{C}/\\Lambda'$  (folded)", fontsize=10)

        for ax in (axL, axR):
            ax.axhline(0, color="k", lw=0.4)
            ax.axvline(0, color="k", lw=0.4)
            ax.set_aspect("equal")
            ax.set_frame_on(False)
            ax.set_xticks([]); ax.set_yticks([])
        # shared limits from the (larger) domain
        allx = [0, 1, 1 + tau_f[0], tau_f[0]]
        ally = [0, 0, tau_f[1], tau_f[1]]
        padf = 0.25
        for ax in (axL, axR):
            ax.set_xlim(min(allx) - padf, max(allx) + padf)
            ax.set_ylim(min(ally) - padf, max(ally) + padf)

        fig_fd.tight_layout()
        st.pyplot(fig_fd)
        plt.close(fig_fd)
        st.caption(
            "Left: the $\\ell$ coloured slabs of $\\mathbb{C}/\\Lambda$, the kernel "
            "points (red), and the fibre of $P$ — its $\\ell$ kernel-translates "
            "(black). Right: all $\\ell$ of them fold onto the single image "
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


# ── Tab: Isogenies over 𝔽ₚ — volcanoes ────────────────────────────────────────
with tab_volc:
    st.subheader("Isogenies over $\\mathbb{F}_p$: Volcanoes")

    st.markdown(
        "Fix a prime $\\ell$ and draw the **$\\ell$-isogeny graph**: one vertex per "
        "curve in the isogeny class (up to isomorphism), with an edge for each "
        "$\\ell$-isogeny. Every vertex has exactly $\\ell + 1$ edges. For ordinary "
        "curves over $\\mathbb{F}_p$, Kohel discovered that this graph always has "
        "the same striking shape — an **isogeny volcano**."
    )

    st.markdown("#### Anatomy of a volcano")
    st.markdown(
        "- **The crater (surface).** At the top is a single cycle of vertices — the "
        "rim. These are the curves whose endomorphism ring is *maximal at $\\ell$*. "
        "The edges around the rim are the **horizontal** isogenies.\n"
        "- **The descending trees.** From each rim vertex hang trees that descend "
        "level by level to the **floor**. Going down one level multiplies the "
        "conductor of $\\mathrm{End}$ by $\\ell$ — the curve gets a 'smaller' "
        "endomorphism ring. Interior vertices have one edge up and the rest down; "
        "the descent stops at the floor.\n"
        "- **The depth** equals the $\\ell$-adic valuation of the conductor "
        "$[\\mathcal{O}_K : \\mathbb{Z}[\\pi]]$, and the **rim length** equals the "
        "order of a prime $\\mathfrak{l} \\mid \\ell$ in the class group of the "
        "surface order."
    )
    st.info(
        "**Why this matters here.** The horizontal isogenies around the rim "
        "*realise the class-group action* — the very same action that permutes the "
        "lattice classes in the analytic picture. Matching the rim to that action "
        "is the heart of the CM bijection (§5). Vertical (descending) isogenies "
        "change the order, so the project's rigid $\\ell$-sets are chosen from "
        "primes that give **horizontal** steps."
    )

    st.divider()

    st.markdown("#### A schematic volcano")
    st.markdown(
        "Adjust the degree, the rim length, and the depth. Colours mark the level: "
        "the crater rim in red, each descending level cooler."
    )

    vol_ctrl, vol_plot = st.columns([1, 2])

    with vol_ctrl:
        vol_ell   = st.select_slider("degree ℓ", options=[2, 3, 5], value=2,
                                     key="vol_ell")
        vol_rim   = st.slider("rim length", 1, 6, 3, key="vol_rim")
        vol_depth = st.slider("depth", 0, 3, 2, key="vol_depth")
        n_surf = vol_rim
        st.caption(
            f"Each vertex has $\\ell + 1 = {vol_ell + 1}$ edges. Rim of "
            f"{vol_rim} vert{'ex' if vol_rim == 1 else 'ices'}, "
            f"depth {vol_depth}."
        )

    with vol_plot:
        R0, dR = 1.0, 1.15
        nodes = []   # (level, angle, radius)
        edges = []   # (i, j)

        def _add(level, angle, radius):
            nodes.append((level, angle, radius))
            return len(nodes) - 1

        crater = [_add(0, 2 * np.pi * k / vol_rim, R0) for k in range(vol_rim)]
        if vol_rim >= 3:
            for k in range(vol_rim):
                edges.append((crater[k], crater[(k + 1) % vol_rim]))
        elif vol_rim == 2:
            edges.append((crater[0], crater[1]))

        def _descend(parent, lo, hi, level, n_children):
            if level > vol_depth or n_children <= 0:
                return
            for j in range(n_children):
                ang = lo + (j + 0.5) / n_children * (hi - lo)
                idx = _add(level, ang, R0 + level * dR)
                edges.append((parent, idx))
                _descend(idx, lo + j / n_children * (hi - lo),
                         lo + (j + 1) / n_children * (hi - lo),
                         level + 1, vol_ell)

        sector = 2 * np.pi / vol_rim
        for k in range(vol_rim):
            a = 2 * np.pi * k / vol_rim
            _descend(crater[k], a - sector / 2, a + sector / 2, 1, vol_ell - 1)

        xy = [(r * np.cos(t), r * np.sin(t)) for (lv, t, r) in nodes]
        lvl_cols = plt.cm.coolwarm_r(
            np.linspace(0.0, 0.85, max(2, vol_depth + 1)))

        fig_vol, ax_vol = plt.subplots(figsize=(6, 6))
        for (i, j) in edges:
            ax_vol.plot([xy[i][0], xy[j][0]], [xy[i][1], xy[j][1]],
                        color="gray", lw=1.0, zorder=1, alpha=0.7)
        for idx, (lv, t, r) in enumerate(nodes):
            ax_vol.scatter(*xy[idx],
                           color=("crimson" if lv == 0 else lvl_cols[lv]),
                           s=(110 if lv == 0 else 70),
                           edgecolor="white", linewidth=0.6, zorder=3)

        lim = R0 + vol_depth * dR + 0.4
        ax_vol.set_xlim(-lim, lim); ax_vol.set_ylim(-lim, lim)
        ax_vol.set_aspect("equal"); ax_vol.set_frame_on(False)
        ax_vol.set_xticks([]); ax_vol.set_yticks([])
        ax_vol.set_title(
            f"$\\ell = {vol_ell}$ volcano — rim {vol_rim}, depth {vol_depth}",
            fontsize=11)
        st.pyplot(fig_vol)
        plt.close(fig_vol)
        st.caption(
            "Red rim = the crater (horizontal isogenies / class-group action). "
            "Cooler colours descend to the floor; each step multiplies the "
            "conductor by $\\ell$. This is a schematic — the live, computed "
            "$\\ell$-isogeny graphs for real classes are on the **Isogeny Class** "
            "page."
        )
