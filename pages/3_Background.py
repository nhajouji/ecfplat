import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon

st.header("Background")
st.markdown(
    "A crash course on the mathematics underlying this project. "
    "Each section pairs a short conceptual overview with an interactive applet."
)

with st.expander("### Algebraic Curves", expanded=False):
    (tab0, tab1, tab2) = st.tabs([
        "Algebraic curves over ℝ",
        "Real Elliptic Curves",
        "Elliptic curves over $\\mathbb{F}_p$",
    ])

with st.expander("### Analytic Methods", expanded=False):
    (tab3, tab4) = st.tabs([
        "Elliptic Curves over ℂ",
        "Endomorphisms and Complex Multiplication",
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
