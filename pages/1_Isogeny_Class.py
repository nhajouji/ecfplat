import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "pycode"))

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

from ecqf_tools import ECQFIsogenyClass, ap_in_pc_data
from graphic_tools import isogeny_graph_figure
from nt import primeQ



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
    p_input = st.number_input("p", value=prefill["p"] if prefill else 5, step=1, min_value=2, max_value=251)
    load = st.button("Load isogeny class", use_container_width=True)

    if prefill:
        load = True

    if load:
        a, p = int(a_input), int(p_input)
        if not ap_in_pc_data((a, p)):
            st.error(
                f"(a, p) = ({a}, {p}) is not in the precomputed data.\n\n"
                "Check that p is prime, a² < 4p, and the pair is within "
                "the precomputed range (4 ≤ p ≤ 251)."
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

# ── Tab 2: Isogeny graph ──────────────────────────────────────────────────────
with tab_isograph:
    st.markdown(
        "Enter a prime ℓ to compute the degree-ℓ isogeny graph. "
        "Rows and columns of the adjacency matrix are in the same order as the table above."
    )

    l_input = st.number_input("ℓ (prime)", min_value=2, value=2, step=1)
    compute = st.button("Compute", use_container_width=True)

    if compute:
        l = int(l_input)
        if not primeQ(l):
            st.error(f"{l} is not prime. Please enter a prime ℓ.")
        else:
            with st.spinner(f"Computing degree-{l} isogeny graph…"):
                mat    = isoclass.adjacency_matrix(l)
                labels = [str(qf) for qf in isoclass.qfs_ordered]
                adj_df = pd.DataFrame(mat, index=labels, columns=labels)
            st.session_state["iso_l"]      = l
            st.session_state["iso_adj_df"] = adj_df

    if "iso_adj_df" in st.session_state and st.session_state.get("iso_l") is not None:
        l      = st.session_state["iso_l"]
        adj_df = st.session_state["iso_adj_df"]

        st.subheader(f"Degree-{l} adjacency matrix")
        st.dataframe(adj_df, use_container_width=True)

        st.subheader("Graph picture")

        # ── Label toggle ──────────────────────────────────────────────────────
        label_choice = st.radio(
            "Node labels",
            ["Index (0, 1, …)", "Quadratic form (a, b, c)", "ec_invs"],
            horizontal=True,
            key="iso_label_mode",
        )

        qfs_ordered = list(isoclass.qfs_ordered)
        if label_choice == "Index (0, 1, …)":
            qf_to_label = {qf: str(i) for i, qf in enumerate(qfs_ordered)}
        elif label_choice == "Quadratic form (a, b, c)":
            qf_to_label = {qf: str(qf) for qf in qfs_ordered}
        else:  # ec_invs
            df_local = isoclass.ecqf_df()
            qf_to_label = {
                tuple(row["qf_coefs"]): str(row["ec_invs"])
                for _, row in df_local.iterrows()
            }
            # Fallback for any qf not in df (shouldn't happen, but just in case)
            for qf in qfs_ordered:
                if qf not in qf_to_label:
                    qf_to_label[qf] = str(qfs_ordered.index(qf))

        fig, _ = isogeny_graph_figure(isoclass, l, qf_to_label)
        st.pyplot(fig)
        plt.close(fig)
