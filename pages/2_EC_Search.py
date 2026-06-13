import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "pycode"))

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

from ecqf_tools import ec_look_up
from graphic_tools import ecfp_classic_plot, ecqf_mw_lattice_plot, ecqf_mwgen_arw_plot


def pts_to_export_str(pts: np.ndarray) -> str:
    lines = []
    for i in range((len(pts) // 3) + 1):
        chunk = pts[3 * i : 3 * (i + 1)]
        if len(chunk) == 0:
            break
        lines.append(",".join(f"[{str(v[0])[:10]},{str(v[1])[:10]}]" for v in chunk))
    return "\n".join(lines)


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
    p_input = st.number_input("p", value=prefill["p"] if prefill else 5, step=1, min_value=5, max_value=251)
    search = st.button("Look up curve", use_container_width=True)

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
    tab_classic    = st.container()
    tab_lattice_tab = None

with tab_classic:
    st.subheader("Classical picture  (affine points in $\\mathbb{F}_p^2$)")
    p = ecdata["ap"][1]
    if p > 100:
        st.warning(
            f"p = {p} is large — the classical plot has {p}² = {p*p} ambient points "
            "and may be slow to render. Showing it anyway."
        )
    fig, ax = ecfp_classic_plot(ecdata)
    # Hide axis ticks for large p to keep the plot readable
    if p >= 15:
        ax.set_xticks([])
        ax.set_yticks([])
    st.pyplot(fig)
    plt.close(fig)

if tab_lattice_tab is not None:
    with tab_lattice_tab:
        st.subheader("Lattice picture")

        k = st.number_input(
            "Frobenius power k  (computes $\\mathbb{F}_{p^k}$–rational points)",
            min_value=1, value=1, step=1,
            key="ec_search_k",
        )

        pic_mode = st.radio(
            "Picture type",
            ["Points", "Generators (arrows)"],
            horizontal=True,
            key="ec_pic_mode",
            help=(
                "Use **Generators** when the point group is large — "
                "it shows the generator directions instead of plotting every point."
            ),
        )

        if pic_mode == "Points":
            fig, ax, pts_arr = ecqf_mw_lattice_plot(ecdata, int(k))
            st.pyplot(fig)
            plt.close(fig)

            f_val, g_val = ecdata["coefs"]
            p_val        = ecdata["ap"][1]
            qf           = ecdata["qf"]
            st.markdown(f"**{len(pts_arr)} point(s)** in the fundamental parallelogram.")
            st.download_button(
                label="Download points as .txt",
                data=pts_to_export_str(pts_arr),
                file_name=(
                    f"points_f{f_val}_g{g_val}_p{p_val}"
                    f"_qf{qf[0]}_{qf[1]}_{qf[2]}_k{k}.txt"
                ),
                mime="text/plain",
                use_container_width=True,
            )
        else:
            fig, ax = ecqf_mwgen_arw_plot(ecdata, int(k))
            st.pyplot(fig)
            plt.close(fig)
