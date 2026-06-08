import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / "pycode"))

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

from ecqf_tools import ECQFIsogenyClass, ap_in_pc_data, abc_to_tau

st.set_page_config(page_title="ecfplat", layout="wide")

# ── Session state ────────────────────────────────────────────────────────────
if "isoclass" not in st.session_state:
    st.session_state.isoclass = None
if "selected_row" not in st.session_state:
    st.session_state.selected_row = None

# ── Sidebar: (a, p) input ────────────────────────────────────────────────────
with st.sidebar:
    st.title("ecfplat")
    st.markdown("Enter a pair *(a, p)* with *p* prime and *a² < 4p*.")
    a_input = st.number_input("a", value=22, step=1)
    p_input = st.number_input("p", value=1021, step=1, min_value=2)
    load = st.button("Load isogeny class", use_container_width=True)

    if load:
        a, p = int(a_input), int(p_input)
        if not ap_in_pc_data((a, p)):
            st.error(
                f"(a, p) = ({a}, {p}) is not in the precomputed data.\n\n"
                "Check that p is prime, a² < 4p, and the pair is within "
                "the precomputed range (4 ≤ p ≤ 1024)."
            )
            st.session_state.isoclass = None
            st.session_state.selected_row = None
        else:
            with st.spinner("Loading…"):
                st.session_state.isoclass = ECQFIsogenyClass(a, p)
                st.session_state.selected_row = None
            st.success(f"Loaded (a, p) = ({a}, {p})")

# ── Main area ────────────────────────────────────────────────────────────────
isoclass = st.session_state.isoclass

if isoclass is None:
    st.info("Enter a pair (a, p) in the sidebar and click **Load isogeny class**.")
    st.stop()

a, p = isoclass.ap
st.header(f"Isogeny class  (a, p) = ({a}, {p})")
st.caption(f"Discriminant: {isoclass.disc}   |   # lattice classes: {len(isoclass.qfs_all)}")

# ── DataFrame tab / Lattice tab ───────────────────────────────────────────────
tab_df, tab_lattice = st.tabs(["Isogeny class table", "Lattice picture"])

# ── Tab 1: DataFrame ──────────────────────────────────────────────────────────
with tab_df:
    st.markdown(
        "Click a row to select it, then switch to the **Lattice picture** tab."
    )
    df = isoclass.ecqf_df()

    # st.dataframe with row selection (Streamlit ≥ 1.35)
    event = st.dataframe(
        df,
        use_container_width=True,
        hide_index=True,
        on_select="rerun",
        selection_mode="single-row",
    )

    selected_rows = event.selection.rows
    if selected_rows:
        st.session_state.selected_row = selected_rows[0]

    if st.session_state.selected_row is not None:
        row = df.iloc[st.session_state.selected_row]
        st.success(
            f"Selected row {st.session_state.selected_row}:  "
            f"j = {row['js']},  fg = {row['fg']},  abc = {row['abc']}"
        )

# ── Tab 2: Lattice picture ────────────────────────────────────────────────────
with tab_lattice:
    if st.session_state.selected_row is None:
        st.info("Select a row in the **Isogeny class table** tab first.")
        st.stop()

    row = df.iloc[st.session_state.selected_row]
    qf = tuple(row["abc"])
    a_qf, b_qf, c_qf = qf

    st.subheader(f"Lattice for abc = {qf}")

    k = st.number_input(
        "Frobenius power k (computes Fp^k–rational points)",
        min_value=1, value=1, step=1
    )

    # ── Compute points ────────────────────────────────────────────────────────
    pts = isoclass.qf_to_mwgr_arr_single(int(k), qf)
    pts_arr = np.array(pts)

    # ── Draw parallelogram + points ───────────────────────────────────────────
    tau = abc_to_tau(qf)
    one = np.array([1.0, 0.0])
    verts = [
        np.array([0.0, 0.0]),
        one,
        one + tau,
        tau,
    ]
    xs = [v[0] for v in verts]
    ys = [v[1] for v in verts]

    fig, ax = plt.subplots(figsize=(5, 5 * tau[1]))
    poly = plt.Polygon(verts, facecolor=[0.85, 0.85, 0.95, 0.4], edgecolor="steelblue", linewidth=1.5)
    ax.add_patch(poly)
    if len(pts_arr) > 0:
        ax.scatter(pts_arr[:, 0], pts_arr[:, 1], s=18, color="steelblue", zorder=3)
    ax.set_xlim(min(xs) - 0.1, max(xs) + 0.1)
    ax.set_ylim(-0.1, max(ys) + 0.1)
    ax.set_aspect("equal")
    ax.set_title(f"abc = {qf},  k = {k},  {len(pts_arr)} point(s)")
    st.pyplot(fig)
    plt.close(fig)

    # ── Export points ─────────────────────────────────────────────────────────
    st.markdown(f"**{len(pts_arr)} point(s)** in the fundamental parallelogram.")

    lines = [f"[{pt[0]:.10f}, {pt[1]:.10f}]" for pt in pts_arr]
    export_text = "\n".join(lines)

    st.download_button(
        label="Download points as .txt",
        data=export_text,
        file_name=f"points_a{a}_p{p}_abc{a_qf}_{b_qf}_{c_qf}_k{k}.txt",
        mime="text/plain",
        use_container_width=True,
    )
