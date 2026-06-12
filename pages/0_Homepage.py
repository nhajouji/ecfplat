import streamlit as st

st.title("Elliptic Curves over F_p via CM Lattices: Tools, Data and Basic Pictures")
st.markdown(
    "Tools for visualizing various aspects of categories of elliptic curves over **F**_p "
    "using precomputed equivalences with a category of CM lattices."
)

st.markdown(
    "The lattice point data produced here is used as input for shader-rendered artwork "
    "by Nadir Hajouji and Steve Trettel. You can view the results at "
    "[elliptic-curves.art](https://elliptic-curves.art/)."
)

st.divider()

col1, col2 = st.columns(2)

with col1:
    st.subheader("Isogeny Class")
    st.markdown(
        "Start from a pair *(a, p)* to load an entire isogeny class. "
        "Browse the bijection table and inspect lattice pictures for each curve."
    )
    st.page_link("pages/1_Isogeny_Class.py", label="Go to Isogeny Class →")

with col2:
    st.subheader("EC Search")
    st.markdown(
        "Start from a specific elliptic curve *y² = x³ + fx + g* over **F**_p. "
        "Look up its trace of Frobenius, associated lattice, and view classical and lattice pictures."
    )
    st.page_link("pages/2_EC_Search.py", label="Go to EC Search →")
