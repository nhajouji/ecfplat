import streamlit as st
from pathlib import Path
from PIL import Image

_LOGO = Path(__file__).parent / "pages" / "ec_logo.jpg"
_icon = Image.open(_LOGO)

st.set_page_config(page_title="ecfplat", page_icon=_icon, layout="wide")


pg = st.navigation([
    st.Page("pages/0_Homepage.py",   title="Homepage",   icon="🏠"),
    st.Page("pages/3_Background.py", title="Background", icon="📖"),
    st.Page("pages/1_Explorer.py",   title="Explorer",   icon="🌋"),
    st.Page("pages/4_SS_Graph.py",   title="Supersingular Graph", icon="🕸️"),
    st.Page("pages/5_Gallery.py",    title="Gallery",    icon="🖼️"),
])
pg.run()
