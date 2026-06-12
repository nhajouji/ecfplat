import streamlit as st

st.set_page_config(page_title="ecfplat", layout="wide")


pg = st.navigation([
    st.Page("pages/0_Homepage.py",      title="Homepage",       icon="🏠"),
    st.Page("pages/3_Background.py",    title="Background",     icon="📖"),
    st.Page("pages/2_EC_Search.py",     title="EC Search",      icon="🔍"),
    st.Page("pages/1_Isogeny_Class.py", title="Isogeny Class",  icon="📋"),
])
pg.run()
