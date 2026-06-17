import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon

st.header("Background")
st.markdown(
    "A crash course on the mathematics underlying this project. "
    "Each section pairs a short conceptual overview with an interactive applet."
)

with st.expander("§1 — Elliptic curves: the basics", expanded=False):
    (tab0, tab1, tab2, tab3) = st.tabs([
        "Algebraic curves over ℝ",
        "Real Elliptic Curves",
        "Elliptic curves over $\\mathbb{F}_p$",
        "Elliptic Curves over ℂ",
    ])

with st.expander("§2 — Frobenius and the lattice pictures", expanded=False):
    (tab4, tab_frob, tab_lift) = st.tabs([
        "Endomorphisms and Complex Multiplication",
        "The Frobenius endomorphism over $\\mathbb{F}_p$",
        "Lifting Frobenius (Deuring)",
    ])

with st.expander("§3 — Isogenies", expanded=False):
    (tab_isog, tab_fold, tab_velu, tab_volc) = st.tabs([
        "Isogenies: kernels and degree",
        "Analytic pictures: folding the torus",
        "Vélu's formulas over $\\mathbb{F}_p$",
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

    # ── Chord-tangent applet ──────────────────────────────────────────────────
    ctrl1, plot1 = st.columns([1, 2])

    with ctrl1:
        st.markdown("**Curve**")
        f_val = st.slider("f", -5.0, 5.0, -1.0, 0.1, key="bg1_f")
        g_val = st.slider("g", -5.0, 5.0,  1.0, 0.1, key="bg1_g")
        disc  = -16 * (4*f_val**3 + 27*g_val**2)
        if abs(disc) < 1e-6:
            st.warning("Singular curve (Δ ≈ 0) — adjust f or g.")

        def point_ui(label, default_x, x_key, s_key):
            st.markdown(f"**Point {label}**")
            x = st.slider(f"x ({label})", -3.0, 3.0, default_x, 0.05, key=x_key)
            y2 = x**3 + f_val*x + g_val
            if y2 < 0:
                st.caption(f"x = {x:.2f} is not on the curve.")
                return x, 0.0, False
            if abs(y2) < 1e-8:
                st.caption(f"{label} is a 2-torsion point (y = 0).")
                return x, 0.0, True
            sgn = st.radio("Branch", ["+", "−"], horizontal=True, key=s_key)
            y = np.sqrt(y2) * (1 if sgn == "+" else -1)
            return x, y, True

        x1, y1, p_ok = point_ui("P", -1.0, "bg1_x1", "bg1_s1")
        x2, y2, q_ok = point_ui("Q",  1.0, "bg1_x2", "bg1_s2")

    with plot1:
        xs  = np.linspace(-3.3, 3.3, 3000)
        ys2 = xs**3 + f_val*xs + g_val
        yp  = np.where(ys2 >= 0, np.sqrt(np.clip(ys2, 0, None)), np.nan)
        yn  = -yp

        fig1, ax1 = plt.subplots(figsize=(6, 5))
        ax1.plot(xs, yp, color="steelblue", lw=2)
        ax1.plot(xs, yn, color="steelblue", lw=2)
        ax1.axhline(0, color="k", lw=0.4)
        ax1.axvline(0, color="k", lw=0.4)
        ax1.set_xlim(-3.3, 3.3)
        ax1.set_ylim(-4.5, 4.5)
        fs = f"{f_val:+.1f}"; gs = f"{g_val:+.1f}"
        ax1.set_title(f"$y^2 = x^3 {fs}x {gs}$", fontsize=11)
        ax1.set_frame_on(False)

        if p_ok and q_ok:
            ax1.scatter([x1], [y1], color="red",   s=70, zorder=6)
            ax1.scatter([x2], [y2], color="green", s=70, zorder=6)
            for x, y, lbl, col in [(x1, y1, "P", "red"), (x2, y2, "Q", "green")]:
                ax1.annotate(lbl, (x, y), xytext=(6, 4),
                             textcoords="offset points",
                             color=col, fontsize=11, fontweight="bold")

            is_double = abs(x1-x2) < 1e-5 and abs(y1-y2) < 1e-5

            if is_double and abs(y1) < 1e-8:
                st.info("P = Q is a 2-torsion point, so 2P = 𝒪.")
            elif not is_double and abs(x1-x2) < 1e-5:
                st.info("P and Q are inverses, so P + Q = 𝒪.")
            else:
                m = (3*x1**2 + f_val)/(2*y1) if is_double else (y2-y1)/(x2-x1)
                lbl_sum = "2P" if is_double else "P+Q"

                x3       = m**2 - x1 - x2
                y3_R     = m*(x3-x1) + y1   # y of R (on line, before reflection)
                y3_sum   = -y3_R

                # chord / tangent line
                xl = np.linspace(-3.3, 3.3, 500)
                yl = m*(xl-x1) + y1
                ax1.plot(xl[np.abs(yl) <= 4.5], yl[np.abs(yl) <= 4.5],
                         color="orange", lw=1.3, ls="--", alpha=0.8, zorder=2)

                if -3.3 <= x3 <= 3.3:
                    ax1.scatter([x3], [y3_R],   color="mediumpurple", s=70, zorder=6)
                    ax1.scatter([x3], [y3_sum],  color="orange",       s=90, zorder=6)
                    ax1.plot([x3, x3], [y3_R, y3_sum],
                             color="gray", ls=":", lw=1.2, zorder=3)
                    ax1.annotate("R", (x3, y3_R), xytext=(6, 4),
                                 textcoords="offset points",
                                 color="mediumpurple", fontsize=11, fontweight="bold")
                    ax1.annotate(lbl_sum, (x3, y3_sum), xytext=(6, -14),
                                 textcoords="offset points",
                                 color="orange", fontsize=11, fontweight="bold")
                    st.markdown(
                        f"**P** = ({x1:.3f}, {y1:.3f})  \n"
                        f"**Q** = ({x2:.3f}, {y2:.3f})  \n"
                        f"**{lbl_sum}** = ({x3:.3f}, {y3_sum:.3f})"
                    )
                else:
                    st.info(f"{lbl_sum} lies outside the visible range (x ≈ {x3:.2f}).")

        st.pyplot(fig1)
        plt.close(fig1)


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

    # ── Applet controls ───────────────────────────────────────────────────────
    _PRIMES = [n for n in range(5, 72)
               if n > 1 and all(n % d != 0 for d in range(2, n))]

    # initialise so variables are always defined
    _p_valid = _q_valid = False
    _x1 = _y1 = _x2 = _y2 = 0
    _pts   = []
    _curve_ok = False

    app_left, app_right = st.columns([1, 2])

    with app_left:
        st.markdown("**Field**")
        _p = st.selectbox("p", _PRIMES,
                          index=_PRIMES.index(17), key="bg2_p")

        st.markdown("**Curve**  $y^2 = x^3 + fx + g$")
        _f = int(st.number_input("f", value=0, step=1, key="bg2_f"))
        _g = int(st.number_input("g", value=1, step=1, key="bg2_g"))

        _disc = (-16 * (4 * pow(_f, 3) + 27 * pow(_g, 2))) % _p
        if _disc == 0:
            st.warning("Singular curve mod p — adjust f or g.")
        else:
            st.success("Smooth curve ✓")
            _curve_ok = True
            _pts      = _ec_pts(_f, _g, _p)
            _xs       = sorted(set(pt[0] for pt in _pts))

            if _xs:
                st.markdown("**Point P**")
                _x1      = st.select_slider("x (P)", options=_xs, key="bg2_x1")
                _y1_opts = sorted(pt[1] for pt in _pts if pt[0] == _x1)
                if len(_y1_opts) == 1:
                    _y1 = _y1_opts[0]; _p_valid = True
                    st.caption(f"y = {_y1}")
                elif len(_y1_opts) >= 2:
                    _y1 = st.radio("y (P)", _y1_opts,
                                   horizontal=True, key="bg2_y1")
                    _p_valid = True

                st.markdown("**Point Q**")
                _x2      = st.select_slider("x (Q)", options=_xs, key="bg2_x2")
                _y2_opts = sorted(pt[1] for pt in _pts if pt[0] == _x2)
                if len(_y2_opts) == 1:
                    _y2 = _y2_opts[0]; _q_valid = True
                    st.caption(f"y = {_y2}")
                elif len(_y2_opts) >= 2:
                    _y2 = st.radio("y (Q)", _y2_opts,
                                   horizontal=True, key="bg2_y2")
                    _q_valid = True
            else:
                st.caption("No affine points on this curve mod p.")

    # ── Plot ──────────────────────────────────────────────────────────────────
    with app_right:
        if not _curve_ok:
            st.info("Adjust f and g to get a smooth curve.")
        else:
            _h    = _p // 2
            _amb  = [(_fp_sym(x, _p), _fp_sym(y, _p))
                     for x in range(_p) for y in range(_p)]

            # Compute line/tangent
            _line = []
            if _p_valid and _q_valid:
                if _x1 == _x2 and _y1 == _y2:
                    if _y1 == 0:
                        st.info("P is a 2-torsion point: 2P = 𝒪.")
                    else:
                        _line = _tangent_pts(_x1, _y1, _f, _p)
                else:
                    _line = _chord_pts(_x1, _y1, _x2, _y2, _p)

            fig_app, ax_app = plt.subplots(figsize=(5, 5))
            # Ambient
            ax_app.scatter([q[0] for q in _amb], [q[1] for q in _amb],
                           color="gray", alpha=0.25, s=8, zorder=1)
            # Fp line (slightly darker gray)
            if _line:
                ax_app.scatter([q[0] for q in _line], [q[1] for q in _line],
                               color="gray", alpha=0.6, s=10, zorder=2)
            # Curve points
            ax_app.scatter([q[0] for q in _pts], [q[1] for q in _pts],
                           color="steelblue", s=18, zorder=3)
            # P and Q
            if _p_valid:
                ax_app.scatter([_x1], [_y1], color="red", s=60, zorder=5)
                ax_app.annotate("P", (_x1, _y1), xytext=(5, 5),
                                textcoords="offset points",
                                color="red", fontsize=10, fontweight="bold")
            if _q_valid:
                ax_app.scatter([_x2], [_y2], color="green", s=60, zorder=5)
                ax_app.annotate("Q", (_x2, _y2), xytext=(5, 5),
                                textcoords="offset points",
                                color="green", fontsize=10, fontweight="bold")

            ax_app.set_xlim(-_h - 0.5, _h + 0.5)
            ax_app.set_ylim(-_h - 0.5, _h + 0.5)
            ax_app.set_aspect("equal")
            ax_app.set_frame_on(False)
            ax_app.set_title(
                f"$y^2 = x^3 + {_f % _p}x + {_g % _p}$  (mod {_p})",
                fontsize=10,
            )
            st.pyplot(fig_app)
            plt.close(fig_app)


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

    # ── Group law applet ──────────────────────────────────────────────────────
    ctrl3, plot3 = st.columns([1, 2])

    with ctrl3:
        st.markdown("**Lattice**")
        tau_re3 = st.slider("Re(τ)", -0.5, 0.5, 0.2, 0.01, key="bg3c_re")
        tau_im3 = st.slider("Im(τ)",  0.1, 3.0, 1.2, 0.05, key="bg3c_im")
        _sign3  = "+" if tau_re3 >= 0 else "-"
        st.latex(rf"\tau = {tau_re3:.2f} {_sign3} {abs(tau_im3):.2f}\,i")

        st.markdown("**Point $z_1 = s_1 + t_1\\,\\tau$**")
        s1 = st.slider("s₁", 0.0, 1.0, 0.2, 0.02, key="bg3c_s1")
        t1 = st.slider("t₁", 0.0, 1.0, 0.3, 0.02, key="bg3c_t1")

        st.markdown("**Point $z_2 = s_2 + t_2\\,\\tau$**")
        s2 = st.slider("s₂", 0.0, 1.0, 0.5, 0.02, key="bg3c_s2")
        t2 = st.slider("t₂", 0.0, 1.0, 0.4, 0.02, key="bg3c_t2")

        # Compute sum
        s3_raw = s1 + s2
        t3_raw = t1 + t2
        s3 = s3_raw % 1.0
        t3 = t3_raw % 1.0
        wrapped = (s3_raw >= 1.0) or (t3_raw >= 1.0)

        st.markdown("**Result $z_1 + z_2$**")
        st.latex(
            rf"s_3 = {s3:.2f},\quad t_3 = {t3:.2f}"
        )
        if wrapped:
            st.caption("(wrapped back into the fundamental domain)")

    with plot3:
        tau3 = np.array([tau_re3, tau_im3])
        one3 = np.array([1.0, 0.0])

        def _frac_to_xy(s, t, tau):
            return s * one3 + t * tau

        z1_xy     = _frac_to_xy(s1, t1, tau3)
        z2_xy     = _frac_to_xy(s2, t2, tau3)
        z3_xy     = _frac_to_xy(s3, t3, tau3)        # reduced
        z3_raw_xy = _frac_to_xy(s3_raw, t3_raw, tau3)  # unreduced

        # Vertices of fundamental domain
        verts3 = [np.zeros(2), one3, one3 + tau3, tau3]

        fig3c, ax3c = plt.subplots(figsize=(6, 6))

        # Draw adjacent copies of the parallelogram (torus context)
        for dm in range(-1, 3):
            for dn in range(-1, 3):
                if dm == 0 and dn == 0:
                    continue
                shift = dm * one3 + dn * tau3
                adj_verts = [v + shift for v in verts3]
                ax3c.add_patch(MplPolygon(
                    adj_verts,
                    facecolor=[0.9, 0.9, 0.9, 0.25],
                    edgecolor="gray", lw=0.8, zorder=1,
                ))

        # Fundamental domain
        ax3c.add_patch(MplPolygon(
            verts3,
            facecolor=[0.85, 0.85, 0.95, 0.5],
            edgecolor="steelblue", lw=2, zorder=2,
        ))

        # If wrapped: show unreduced position and translation arrow
        if wrapped:
            ax3c.scatter(*z3_raw_xy, color="orange", s=80, alpha=0.35,
                         zorder=4, marker="o")
            ax3c.annotate(
                "", xy=z3_xy, xytext=z3_raw_xy,
                arrowprops=dict(arrowstyle="->", color="orange",
                                lw=1.5, linestyle="dashed"),
                zorder=5,
            )

        # Points
        for xy, lbl, col in [
            (z1_xy, "$z_1$", "red"),
            (z2_xy, "$z_2$", "green"),
            (z3_xy, "$z_1+z_2$", "orange"),
        ]:
            ax3c.scatter(*xy, color=col, s=90, zorder=6)
            ax3c.annotate(lbl, xy, xytext=(6, 5),
                          textcoords="offset points",
                          color=col, fontsize=11, fontweight="bold")

        # Axis limits: cover fundamental domain plus context copies
        pad3 = 0.3
        all_x = [v[0] for v in verts3] + [z3_raw_xy[0]]
        all_y = [v[1] for v in verts3] + [z3_raw_xy[1]]
        ax3c.set_xlim(min(all_x) - pad3, max(all_x) + pad3 + 1)
        ax3c.set_ylim(-pad3, max(all_y) + pad3 + tau3[1] * 0.5)
        ax3c.set_aspect("equal")
        ax3c.set_frame_on(False)
        ax3c.axhline(0, color="k", lw=0.4)
        ax3c.axvline(0, color="k", lw=0.4)
        ax3c.set_title("Group law on $\\mathbb{C}/\\Lambda$", fontsize=11)

        st.pyplot(fig3c)
        plt.close(fig3c)


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
        "We now have two pictures of an elliptic curve with extra endomorphisms:\n\n"
        "- Over $\\mathbb{C}$, a CM curve is a **lattice** $\\Lambda$, and each "
        "endomorphism is multiplication by a complex number $\\alpha$ with "
        "$\\alpha\\Lambda \\subseteq \\Lambda$.\n"
        "- Over $\\mathbb{F}_p$, an ordinary curve carries the **Frobenius** "
        "endomorphism $\\phi$, a root of $\\phi^2 - a\\phi + p = 0$ living in an "
        "imaginary quadratic order.\n\n"
        "These two pictures are connected by a theorem of **Deuring**."
    )

    st.markdown("#### Deuring's lifting theorem")
    st.info(
        "Let $E/\\mathbb{F}_p$ be an **ordinary** elliptic curve with "
        "$\\mathrm{End}(E) = \\mathcal{O}$, an order in the imaginary quadratic "
        "field $K = \\mathbb{Q}(\\sqrt{a^2 - 4p})$. Then there is an elliptic "
        "curve $\\tilde E$ defined over a number field, with "
        "$\\mathrm{End}(\\tilde E) = \\mathcal{O}$, that **reduces to $E$ mod $p$** "
        "— and the Frobenius $\\phi \\in \\mathcal{O}$ is the reduction of an honest "
        "complex-analytic endomorphism of $\\tilde E$."
    )
    st.markdown(
        "Concretely: the curve $E$ over $\\mathbb{F}_p$ corresponds to a lattice "
        "$\\Lambda$ in $\\mathbb{C}$ with CM by $\\mathcal{O}$, and Frobenius lifts "
        "to **multiplication by the complex number**"
    )
    st.latex(r"\alpha = \frac{a + \sqrt{a^2 - 4p}}{2}, \qquad |\alpha|^2 = \alpha\bar\alpha = p.")
    st.markdown(
        "This $\\alpha$ is just $\\phi$ viewed inside $\\mathbb{C}$ via "
        "$\\mathcal{O} \\hookrightarrow \\mathbb{C}$ — it generates "
        "$\\mathcal{O} = \\mathbb{Z}[\\alpha]$ and satisfies the very same equation "
        "$\\alpha^2 - a\\alpha + p = 0$. We **take the existence of $\\Lambda$ for "
        "granted here** — producing it explicitly, and matching every curve in the "
        "isogeny class to its lattice, is exactly the *CM bijection* worked out "
        "later. For now the payoff is the picture: a curve over $\\mathbb{F}_p$ "
        "becomes a lattice in the plane, carrying a distinguished "
        "rotation-and-scaling $\\alpha$."
    )

    st.markdown("#### This *is* the lattice picture")
    st.markdown(
        "The artwork this project feeds is built from exactly this data. Each curve "
        "in an isogeny class lifts to a lattice $\\Lambda$ with CM by $\\mathcal{O}$, "
        "and all of them share the same marked endomorphism — multiplication by "
        "$\\alpha$ (with $|\\alpha| = \\sqrt p$). Drawing those lattices — and the "
        "action of $\\alpha$ on them — is how we 'draw an elliptic curve over "
        "$\\mathbb{F}_p$.' The applet below shows a single lattice with its "
        "Frobenius action."
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
        "the same point count. The whole project is the study of one such "
        "**isogeny class** at a time."
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

    # ── Small illustration: a sublattice ──────────────────────────────────────
    st.markdown("#### A finer lattice")
    st.markdown(
        "Below, $\\Lambda$ (large blue dots) sits inside an index-3 overlattice "
        "$\\Lambda'$ (all dots). The kernel $\\Lambda'/\\Lambda$ has the 3 cosets "
        "highlighted inside one fundamental cell of $\\Lambda$ — these are the 3 "
        "points killed by the 3-isogeny $\\mathbb{C}/\\Lambda \\to \\mathbb{C}/\\Lambda'$."
    )

    _tau_ill = np.array([0.35, 1.0])
    _one_ill = np.array([1.0, 0.0])
    _ell_ill = 3
    # Lambda' = (1/ell) Z + tau Z  ⊃  Lambda = Z + tau Z
    fig_il, ax_il = plt.subplots(figsize=(6, 4))
    RNG_I = 3
    for m in range(-RNG_I, RNG_I + 1):
        for n in range(-2, 4):
            # fine lattice point: (m/ell)·1 + n·tau
            P_fine = (m / _ell_ill) * _one_ill + n * _tau_ill
            on_coarse = (m % _ell_ill == 0)
            ax_il.scatter(*P_fine,
                          color="steelblue" if on_coarse else "lightgray",
                          s=70 if on_coarse else 22,
                          zorder=3 if on_coarse else 2)
    # one fundamental cell of Lambda, with the 3 kernel cosets
    ax_il.add_patch(MplPolygon(
        [np.zeros(2), _one_ill, _one_ill + _tau_ill, _tau_ill],
        facecolor=[0.85, 0.85, 0.95, 0.30], edgecolor="steelblue",
        lw=1.5, zorder=1))
    for k in range(_ell_ill):
        P_k = (k / _ell_ill) * _one_ill
        ax_il.scatter(*P_k, color="crimson", s=90, zorder=5,
                      edgecolor="white", linewidth=0.8)
    ax_il.annotate("kernel\n$\\Lambda'/\\Lambda$", (1.0 / _ell_ill, 0),
                   xytext=(10, 22), textcoords="offset points",
                   color="crimson", fontsize=10, fontweight="bold",
                   ha="left")
    ax_il.set_xlim(-1.2, 2.0)
    ax_il.set_ylim(-1.3, 2.2)
    ax_il.set_aspect("equal")
    ax_il.axhline(0, color="k", lw=0.4)
    ax_il.axvline(0, color="k", lw=0.4)
    ax_il.set_frame_on(False)
    ax_il.set_title("$\\Lambda \\subseteq \\Lambda'$  (index 3)", fontsize=11)
    st.pyplot(fig_il)
    plt.close(fig_il)


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
