import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "pycode"))

import streamlit as st

from ecqf_tools import ECQFIsogenyClass, ap_in_pc_data

# ── Session state ─────────────────────────────────────────────────────────────
if "isoclass" not in st.session_state:
    st.session_state.isoclass = None
if "selected_row" not in st.session_state:
    st.session_state.selected_row = None

# ── Pre-fill from EC Search navigation ───────────────────────────────────────
prefill = st.session_state.pop("ic_prefill", None)

# ── Sidebar: (a, p) input ─────────────────────────────────────────────────────
with st.sidebar:
    st.title("Isogeny Class")
    st.markdown("Enter a pair *(a, p)* with *p* prime and *a² < 4p*.")
    a_input = st.number_input("a", value=prefill["a"] if prefill else -4, step=1)
    p_input = st.number_input("p", value=prefill["p"] if prefill else 5, step=1, min_value=2)
    load = st.button("Load isogeny class", use_container_width=True)

    # Auto-load when arriving via navigation link
    if prefill:
        load = True

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

# ── Main area ─────────────────────────────────────────────────────────────────
isoclass = st.session_state.isoclass

if isoclass is None:
    st.info("Enter a pair (a, p) in the sidebar and click **Load isogeny class**.")
    st.stop()

a, p = isoclass.ap
st.header(f"Isogeny class  (a, p) = ({a}, {p})")
st.caption(f"Discriminant: {isoclass.disc}   |   # lattice classes: {len(isoclass.qfs_all)}")

# ── Tabs ──────────────────────────────────────────────────────────────────────
tab_df, tab_isograph = st.tabs(["Isogeny class table", "Isogeny graph"])

# ── Tab 1: DataFrame ──────────────────────────────────────────────────────────
with tab_df:
    st.markdown("Click a row to select it and view the corresponding elliptic curve.")
    df = isoclass.ecqf_df()

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
        f_sel, g_sel = tuple(row["EC_coefs"])
        st.success(
            f"Selected row {st.session_state.selected_row}:  "
            f"ec_invs = {row['ec_invs']},  EC_coefs = {row['EC_coefs']},  qf_coefs = {row['qf_coefs']}"
        )
        if st.button("View this curve in EC Search →", use_container_width=True):
            st.session_state.ec_prefill = {"f": int(f_sel), "g": int(g_sel), "p": p}
            st.switch_page("pages/2_EC_Search.py")

# ── Tab 2: Isogeny graph (placeholder) ───────────────────────────────────────
with tab_isograph:
    st.info(
        "Isogeny graph data is not yet implemented. "
        "This tab will display ℓ-isogeny graphs between the curves in this class."
    )
