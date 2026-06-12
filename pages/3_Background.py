import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon

st.header("Background")
st.markdown(
    "A crash course on the mathematics underlying this project. "
    "Each section pairs a short conceptual overview with an interactive applet."
)

(tab1, tab2, tab3, tab4, tab5, tab6) = st.tabs([
    "EC over ℝ (algebraic)",
    "EC over 𝔽ₚ (algebraic)",
    "EC over ℂ (analytic)",
    "EC over ℝ (analytic)",
    "Isogenies, Endomorphisms & CM",
    "EC over 𝔽ₚ via Frobenius",
])


# ── Tab 1: Elliptic curves over ℝ (algebraically) ────────────────────────────
with tab1:
    st.subheader("Elliptic Curves over ℝ — Algebraic Viewpoint")
    st.markdown("An **elliptic curve** over ℝ is the set of points satisfying")
    st.latex(r"y^2 = x^3 + fx + g, \qquad f, g \in \mathbb{R},")
    st.markdown(
        "together with a distinguished 'point at infinity' 𝒪. "
        "The curve is **smooth** when the discriminant Δ = −16(4f³ + 27g²) ≠ 0.\n\n"
        "The set of points forms an **abelian group** with 𝒪 as the identity. "
        "Given two points P and Q, their sum P + Q is constructed by the "
        "**chord-tangent law**:\n"
        "1. Draw the line through P and Q (or the tangent at P if P = Q).\n"
        "2. This line meets the curve in a third point R.\n"
        "3. Reflect R over the x-axis to obtain P + Q."
    )
    st.divider()

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


# ── Tab 2: Elliptic curves over 𝔽ₚ (algebraically) ──────────────────────────
with tab2:
    st.subheader("Elliptic Curves over 𝔽ₚ — Algebraic Viewpoint")
    st.markdown(
        "Over the finite field 𝔽ₚ the same equation y² = x³ + fx + g defines a "
        "**finite** set of points, and the chord-tangent law still makes the set into "
        "an abelian group — but now all arithmetic is done mod p. "
        "The EC Search page already lets you plot the affine points of any such curve."
    )
    st.info("An interactive group-law applet for 𝔽ₚ will be added here.")
    st.page_link("pages/2_EC_Search.py", label="Go to EC Search →")


# ── Tab 3: Elliptic curves over ℂ (analytically) ─────────────────────────────
with tab3:
    st.subheader("Elliptic Curves over ℂ — Analytic Viewpoint")
    st.markdown("Every elliptic curve over ℂ is isomorphic, as a complex manifold, to a **complex torus**")
    st.latex(r"\mathbb{C}/\Lambda, \qquad \Lambda = \mathbb{Z} + \tau\,\mathbb{Z},")
    st.markdown(
        "where τ lies in the **upper half-plane** (Im τ > 0). "
        "Two values of τ give isomorphic tori if and only if they are related by a "
        "Möbius transformation in SL(2,ℤ).\n\n"
        "Use the sliders below to move τ and watch the fundamental parallelogram — "
        "the basic building block of the lattice — change shape."
    )
    st.divider()

    ctrl3, plot3 = st.columns([1, 2])
    with ctrl3:
        tau_re = st.slider("Re(τ)", -0.5, 0.5,  0.2, 0.01, key="bg3_re")
        tau_im = st.slider("Im(τ)",  0.1, 3.0,   1.2, 0.05, key="bg3_im")
        sign   = "+" if tau_re >= 0 else "-"
        st.latex(rf"\tau = {tau_re:.2f} {sign} {abs(tau_im):.2f}\,i")

    with plot3:
        tau  = np.array([tau_re, tau_im])
        one  = np.array([1.0, 0.0])
        verts = [np.zeros(2), one, one + tau, tau]

        N  = 5
        lx = [m + n*tau_re for m in range(-N, N+1) for n in range(-N, N+1)]
        ly = [n*tau_im      for m in range(-N, N+1) for n in range(-N, N+1)]

        fig3, ax3 = plt.subplots(figsize=(6, 6))
        ax3.scatter(lx, ly, s=14, color="steelblue", alpha=0.4, zorder=2)

        poly3 = MplPolygon(verts, facecolor=[0.85, 0.85, 0.95, 0.5],
                           edgecolor="steelblue", lw=2, zorder=3)
        ax3.add_patch(poly3)

        for pt, lbl, col in [
            (np.zeros(2), "0",   "red"),
            (one,         "1",   "black"),
            (tau,         "τ",   "orange"),
            (one + tau,   "1+τ", "black"),
        ]:
            ax3.scatter(*pt, color=col, s=55, zorder=5)
            ax3.annotate(lbl, pt, xytext=(5, 5),
                         textcoords="offset points", fontsize=11, color=col)

        pad = 0.4
        ax3.set_xlim(min(lx) - pad, max(lx) + pad)
        ax3.set_ylim(-pad,           max(ly) + pad)
        ax3.set_aspect("equal")
        ax3.set_frame_on(False)
        ax3.axhline(0, color="k", lw=0.4)
        ax3.axvline(0, color="k", lw=0.4)
        ax3.set_title("Lattice  Λ = ℤ + τℤ", fontsize=11)
        st.pyplot(fig3)
        plt.close(fig3)


# ── Tab 4: Elliptic curves over ℝ (analytically) ─────────────────────────────
with tab4:
    st.subheader("Elliptic Curves over ℝ — Analytic Viewpoint")
    st.info(
        "This section will explain how the real locus of an elliptic curve sits inside "
        "the complex torus ℂ/Λ, and how the real and imaginary periods determine the "
        "real geometry of the curve."
    )


# ── Tab 5: Isogenies, Endomorphisms & CM ─────────────────────────────────────
with tab5:
    st.subheader("Isogenies, Endomorphisms and Complex Multiplication")
    st.info(
        "This section will cover: isogenies as group homomorphisms between elliptic "
        "curves, the endomorphism ring End(E), and the special case of complex "
        "multiplication where End(E) is strictly larger than ℤ."
    )


# ── Tab 6: EC over 𝔽ₚ via Frobenius ─────────────────────────────────────────
with tab6:
    st.subheader("Elliptic Curves over 𝔽ₚ via Lifts of Frobenius")
    st.info(
        "This section will explain the central bijection behind this project: how the "
        "Frobenius endomorphism of E/𝔽ₚ determines a CM lattice class, giving an "
        "explicit correspondence between isogeny classes and classes of binary "
        "quadratic forms. The tools in the main app are built on this bijection."
    )
    st.markdown("**Explore the bijection in the main app:**")
    st.page_link("pages/1_Isogeny_Class.py", label="Go to Isogeny Class →")
    st.page_link("pages/2_EC_Search.py",     label="Go to EC Search →")
