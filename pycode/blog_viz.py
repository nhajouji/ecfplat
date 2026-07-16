"""Interactive widgets for the Blog page.

Each builder returns a self-contained HTML document (canvas + JS) to be dropped
into a Streamlit page via ``st.components.v1.html``.  Unlike the other ``*_viz``
modules these mostly load prebuilt, hand-authored applets from
``pages/blog_assets/`` so a post can embed the exact widget it was written
around.  The applets carry their own theming and small validated data; they do
not touch the app's precomputed stores.
"""
from pathlib import Path

_ASSETS = Path(__file__).parent.parent / "pages" / "blog_assets"


def _load(name: str) -> str:
    return (_ASSETS / name).read_text(encoding="utf-8")


def sos_descent_html() -> str:
    """"The Shortest Vector" -- pick a prime p = 1 (mod 4) and watch the
    Euclidean algorithm in Z[i] descend to the shortest lattice vector, the
    sum of two squares p = a^2 + b^2.  (Post 1: sum of two squares.)"""
    return _load("sos-descent-explorer.html")


def three_languages_square_html() -> str:
    """"One square, three languages" -- the l=2 isogeny graph of discriminant
    -39 (class group Z/4) as a 4-cycle, labelled at once as a class group, a
    fusion graph of representations, and a system of particles.  Toggling the
    reflection across the real axis leaves the graph unchanged: the invisible
    mirror that motivates rigidity.  (Post 2: three languages.)"""
    return _load("three-languages-square.html")


def three_languages_prism_html() -> str:
    """The hexagonal prism of Z/2 x Z/6 -- discriminant -231, the non-cyclic
    example where the spanning set genuinely needs two generators (l=2 of order
    6 for the hexagons, l=3 of order 2 for the rungs).  Same three-track
    labelling; the vertical mirror swaps matter and antimatter.  (Post 2.)"""
    return _load("three-languages-prism.html")
