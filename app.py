import streamlit as st
from pathlib import Path
from PIL import Image

_LOGO = Path(__file__).parent / "pages" / "ec_logo.jpg"
_icon = Image.open(_LOGO)

st.set_page_config(page_title="ecfplat", page_icon=_icon, layout="wide")


pg = st.navigation([
    st.Page("pages/0_Homepage.py",      title="Homepage",       icon="🏠"),
    st.Page("pages/3_Background.py",    title="Background",     icon="📖"),
    st.Page("pages/5_Modular_Curves.py", title="Modular Curves", icon="🌀"),
    st.Page("pages/2_EC_Search.py",     title="EC Search",      icon="🔍"),
    st.Page("pages/1_Isogeny_Class.py", title="Isogeny Class",  icon="📋"),
    st.Page("pages/4_SS_Graph.py",      title="Supersingular Graph", icon="🕸️"),
])
pg.run()
