import streamlit as st
from pathlib import Path

st.set_page_config(page_title="ecfplat", layout="wide")

LOGO = Path(__file__).parent / "pages" / "ec_logo.jpg"
st.sidebar.image(str(LOGO), use_container_width=True)

pg = st.navigation([
    st.Page("pages/0_Homepage.py",      title="Homepage",       icon="🏠"),
    st.Page("pages/1_Isogeny_Class.py", title="Isogeny Class",  icon="📋"),
    st.Page("pages/2_EC_Search.py",     title="EC Search",      icon="🔍"),
])
pg.run()
