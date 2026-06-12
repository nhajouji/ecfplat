import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "pycode"))

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

from ecqf_tools import ec_look_up, abc_to_tau, frob_to_mw_gens, mw_arr_from_gens

# ── Helpers ───────────────────────────────────────────────────────────────────

def fp_to_sym_fd(x: int, p: int) -> int:
    xr = x % p
    return xr if 2 * xr < p else xr - p


def ecfp_aff_pts(fg: tuple, p: int) -> list:
    square_root_dic = {y2: [] for y2 in range(p)}
    for x in range(p):
        square_root_dic[pow(x, 2, p)].append(x)
    pts = []
    f, g = fg
    for x in range(p):
        x = fp_to_sym_fd(x, p)
        y2 = (pow(x, 3) + f * x + g) % p
        pts += [(x, fp_to_sym_fd(y, p)) for y in square_root_dic[y2]]
    return pts


def classical_plot(ecdata: dict):
    fg = ecdata["coefs"]
    p = ecdata["ap"][1]
    xys = ecfp_aff_pts(fg, p)
    xs = [xy[0] for xy in xys]
    ys = [xy[1] for xy in xys]
    ambient = [(x, y) for x in range(-(p // 2), (p // 2) + 1)
                       for y in range(-(p // 2), (p // 2) + 1)]
    ax_s, ay_s = zip(*ambient)
    fig, ax = plt.subplots()
    ax.scatter(ax_s, ay_s, color="gray", alpha=0.3, s=8)
    ax.scatter(xs, ys, color="purple", s=18, zorder=3)
    ax.scatter([(p + 1) // 2], [(p + 1) // 2], color="olive", s=30, zorder=4, label="∞")
    if p <15:
        ax.set_xticks([x for x in range(-(p // 2), (p + 1) // 2)])
        ax.set_yticks([y for y in range(-(p // 2), (p + 1) // 2)])
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    ax.set_xlim(-(p + 3) // 2, (p + 3) // 2)
    ax.set_ylim(-(p + 3) // 2, (p + 3) // 2)
    ax.set_aspect("equal")
    ax.set_frame_on(False)
    ax.set_title(ecdata["ec_eq"])
    return fig


def pts_to_export_str(pts: np.ndarray) -> str:
    lines = []
    for i in range((len(pts) // 3) + 1):
        chunk = pts[3 * i : 3 * (i + 1)]
        if len(chunk) == 0:
            break
        lines.append(",".join(f"[{str(v[0])[:10]},{str(v[1])[:10]}]" for v in chunk))
    return "\n".join(lines)


def lattice_plot(ecdata: dict, k: int):
    """Returns (fig, pts_array)."""
    qf = ecdata["qf"]
    frm = ecdata["FrobMat"]
    mw_gens = frob_to_mw_gens(frm, k)
    pts = mw_arr_from_gens(qf, mw_gens)
    pts_arr = np.array(pts) if len(pts) > 0 else np.empty((0, 2))
    tau = abc_to_tau(qf)
    one = np.array([1.0, 0.0])
    verts = [np.array([0.0, 0.0]), one, one + tau, tau]
    xs = [v[0] for v in verts]
    ys = [v[1] for v in verts]
    fig, ax = plt.subplots(figsize=(5, 5 * tau[1]))
    poly = plt.Polygon(verts, facecolor=[0.85, 0.85, 0.95, 0.4],
                       edgecolor="steelblue", linewidth=1.5)
    ax.add_patch(poly)
    if len(pts_arr) > 0:
        ax.scatter(pts_arr[:, 0], pts_arr[:, 1], s=18, color="steelblue", zorder=3)
    ax.scatter([0], [0], color="olive", s=30, zorder=4)
    ax.set_xlim(min(xs) - 0.1, max(xs) + 0.1)
    ax.set_ylim(-0.1, max(ys) + 0.1)
    ax.set_aspect("equal")
    ax.set_frame_on(False)
    ax.set_title(f"τ = {ecdata['tau_str']},  k = {k},  {len(pts_arr)} point(s)")
    return fig, pts_arr


# ── Session state ─────────────────────────────────────────────────────────────
if "ecdata" not in st.session_state:
    st.session_state.ecdata = None

# ── Pre-fill from Isogeny Class navigation ────────────────────────────────────
prefill = st.session_state.pop("ec_prefill", None)

# ── Sidebar: curve input ──────────────────────────────────────────────────────
with st.sidebar:
    st.title("EC Search")
    st.markdown(
        "Enter coefficients *(f, g)* and prime *p* for the curve\n\n"
        "**y² = x³ + fx + g  (mod p)**"
    )
    f_input = st.number_input("f", value=prefill["f"] if prefill else 3, step=1)
    g_input = st.number_input("g", value=prefill["g"] if prefill else 0, step=1)
    p_input = st.number_input("p", value=prefill["p"] if prefill else 5, step=1, min_value=5)
    search = st.button("Look up curve", use_container_width=True)

    # Auto-lookup when arriving via navigation link
    if prefill:
        search = True

    if search:
        f, g, p = int(f_input), int(g_input), int(p_input)
        try:
            with st.spinner("Looking up…"):
                st.session_state.ecdata = ec_look_up((f, g), p)
            st.success("Curve found.")
        except ValueError as e:
            st.error(str(e))
            st.session_state.ecdata = None

# ── Main area ─────────────────────────────────────────────────────────────────
ecdata = st.session_state.ecdata

if ecdata is None:
    st.info("Enter curve coefficients in the sidebar and click **Look up curve**.")
    st.stop()

st.header(ecdata["ec_eq"])

# ── Summary ───────────────────────────────────────────────────────────────────
col_l, col_r = st.columns(2)
with col_l:
    st.markdown(f"**j-invariant:** {ecdata['j']}")
    st.markdown(f"**Signature (s):** {ecdata['s']}")
    st.markdown(f"**Trace of Frobenius (a):** {ecdata['frob_tr']}")
    st.markdown(f"**Frobenius discriminant (a²−4p):** {ecdata['frob_disc']}")
    st.markdown(f"**Frobenius conductor:** {ecdata['frob_cond']}")
    st.markdown(f"**Fundamental endomorphism discriminant:** {ecdata['endo_fun_disc']}")

with col_r:
    if ecdata["has_pcqf"]:
        st.markdown(f"**Quadratic form (a, b, c):** {ecdata['qf']}")
        st.markdown(f"**τ:** {ecdata['tau_str']}")
        st.markdown(f"**τ coordinates:** {ecdata['tau_xy']}")
        st.markdown(f"**Frobenius matrix:** {ecdata['FrobMat'].vec}")
        a_frob, p_char = ecdata["ap"]
        if st.button("View isogeny class →", use_container_width=True):
            st.session_state.ic_prefill = {"a": a_frob, "p": p_char}
            st.switch_page("pages/1_Isogeny_Class.py")
    else:
        st.warning(
            "No precomputed quadratic form found for this curve. "
            "Only basic data is available (p must be ≤ 1024 and in the precomputed range)."
        )

st.divider()

# ── Picture tabs ──────────────────────────────────────────────────────────────
if ecdata["has_pcqf"]:
    tab_classic, tab_lattice_tab = st.tabs(["Classical picture", "Lattice picture"])
else:
    tab_classic = st.container()
    tab_lattice_tab = None

with tab_classic:
    st.subheader("Classical picture  (affine points in **F**_p²)")
    p = ecdata["ap"][1]
    if p > 100:
        st.warning(
            f"p = {p} is large — the classical plot has {p}² = {p*p} ambient points "
            "and may be slow to render. Showing it anyway."
        )
    fig = classical_plot(ecdata)
    st.pyplot(fig)
    plt.close(fig)

if tab_lattice_tab is not None:
    with tab_lattice_tab:
        st.subheader("Lattice picture")
        k = st.number_input(
            "Frobenius power k (computes **F**_{p^k}–rational points)",
            min_value=1, value=1, step=1,
            key="ec_search_k"
        )
        fig, pts_arr = lattice_plot(ecdata, int(k))
        st.pyplot(fig)
        plt.close(fig)

        f_val, g_val = ecdata["coefs"]
        p_val = ecdata["ap"][1]
        qf = ecdata["qf"]
        st.markdown(f"**{len(pts_arr)} point(s)** in the fundamental parallelogram.")
        st.download_button(
            label="Download points as .txt",
            data=pts_to_export_str(pts_arr),
            file_name=f"points_f{f_val}_g{g_val}_p{p_val}_qf{qf[0]}_{qf[1]}_{qf[2]}_k{k}.txt",
            mime="text/plain",
            use_container_width=True,
        )
