"""The Gallery — a picture wall.

Grid of thumbnails grouped by venue (ICM, then Bridges); every picture
links to its own detail view (?pic=slug) with the full-size render, a
short note, and links into the Explorer / Background. No teaching here —
that's what the Background is for.
"""

import sys
import base64
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "pycode"))

import streamlit as st

IMG = Path(__file__).parent / "gallery_img"

# ── the pictures ──────────────────────────────────────────────────────────────
# Each entry: slug (thumb/<slug>.jpg + full/<slug>.jpg), title, blurb
# (1-2 sentences, shown on the detail page and as hover text), links
# [(label, href)], optional credit override.

CREDIT_DEFAULT = "Nadir Hajouji & Steve Trettel"

ICM = [
    {"slug": "icm-a1", "title": "Trace 1 over 𝔽₅",
     "blurb": "The isogeny class (a, p) = (1, 5): a single curve with 5 "
              "rational points and CM by discriminant −19.",
     "links": [("open in the Explorer", "/Explorer?a=1&p=5"),
               ("disc −19", "/Explorer?d=-19")]},
    {"slug": "icm-a2", "title": "Trace 2 over 𝔽₅ — the column",
     "blurb": "Two curves joined by a 2-isogeny: j = 1728 above (disc −4), "
              "its descendant below (disc −16). Hung vertically so ascending "
              "and descending mean what they say.",
     "links": [("open in the Explorer", "/Explorer?a=2&p=5"),
               ("disc −16", "/Explorer?d=-16")]},
    {"slug": "icm-a3", "title": "Trace 3 over 𝔽₅",
     "blurb": "One curve, three rational points, CM by discriminant −11.",
     "links": [("open in the Explorer", "/Explorer?a=3&p=5"),
               ("disc −11", "/Explorer?d=-11")]},
    {"slug": "icm-a4", "title": "Trace 4 over 𝔽₅",
     "blurb": "j = 1728 again — the ℤ[i] lattice — now with trace 4 and just "
              "two rational points.",
     "links": [("open in the Explorer", "/Explorer?a=4&p=5"),
               ("disc −4", "/Explorer?d=-4")]},
    {"slug": "icm-a0", "title": "Trace 0 over 𝔽₅ — supersingular",
     "blurb": "The supersingular class (a, p) = (0, 5): y² = x³ + 1, with j = 0 "
              "and no ordinary CM lift — the one exceptional class over 𝔽₅. Made "
              "as a lenticular print that flips between the untwisted torus and "
              "its twisted form as you move past it.",
     "links": [("open in the Explorer", "/Explorer?a=0&p=5"),
               ("disc −20", "/Explorer?d=-20")]},
    {"slug": "icm-pcomplex", "title": "X₀(11) over ℂ",
     "blurb": "The complex points of X₀(11) form a torus; the real locus is "
              "traced in green, the marked point in red.",
     "links": [("modular curves in the Background", "/Background")]},
    {"slug": "icm-p23", "title": "X₀(11) mod 23",
     "blurb": "Reduced mod 23 the curve has trace −1 and discriminant −91 — "
              "class number 2, so Frobenius has two lifts.",
     "links": [("open in the Explorer", "/Explorer?a=-1&p=23"),
               ("disc −91", "/Explorer?d=-91")]},
    {"slug": "icm-p101", "title": "X₀(11) mod 101",
     "blurb": "Trace 2, discriminant −400, one hundred rational points; the "
              "eight lattice classes of this discriminant spread over three "
              "levels of the volcano.",
     "links": [("open in the Explorer", "/Explorer?a=2&p=101"),
               ("disc −400", "/Explorer?d=-400")]},
    {"slug": "icm-p107", "title": "X₀(11) mod 107 — the six lifts",
     "blurb": "Mod 107 the Frobenius order is ℤ[√−26], with class number 6: "
              "six lattices, six lifts, one curve. The class group is cyclic "
              "and the arrangement encodes it — the tori sit like the 6th "
              "roots of unity.",
     "links": [("open in the Explorer", "/Explorer?a=18&p=107"),
               ("disc −104", "/Explorer?d=-104")]},
    {"slug": "icm-x05", "title": "X₀(5) as a real surface",
     "blurb": "A slice of the universal curve over X₁(5) — every thread a "
              "real elliptic curve. Four fibers are highlighted: the two "
              "lattices with CM by √−5, and the pair X₁(11) ⇄ X₀(11), whose "
              "degree-5 map is itself a point of X₀(5).",
     "links": [("the supersingular class (0, 5)", "/Explorer?a=0&p=5"),
               ("disc −20", "/Explorer?d=-20"),
               ("modular curves in the Background", "/Background")]},
    {"slug": "icm-starscape",
     "title": "Elliptic Starscape in the complex projective plane",
     "blurb": "A guest piece by Elliot Kienzle (chessapig): the modular curve "
              "X₀(11) drawn as a night sky in ℂℙ².",
     "desc": (
        r"Illustration of the modular curve $X_0(11)$ — the solutions to "
        r"$y^2 - y = x^3 - x^2 - 10x - 20$ — as a subset of $\mathbb{CP}^2$. "
        r"The white points are projections of complex points, plotting the real "
        r"part of the point in a (nonstandard) affine chart. The blue curve shows "
        r"the projected real locus. The points were computed by intersecting "
        r"$X_0(11)$ with ~100,000 complex lines in $\mathbb{CP}^2$ of rational "
        r"slope, then sized according to the slope's denominator. The resulting "
        r"density of points is induced from the rotation-invariant metric on "
        r"$\mathbb{CP}^2$. Under this projection, a complex line appears as an "
        r"elliptical cloud. Notice how the cubic curve is a smoothing of three "
        r"complex lines: a strong diagonal through the central point, the wider "
        r"horizontal line encircling the central point, and the diffuse "
        r"background parallel to the screen."),
     "credit": "Elliot Kienzle (chessapig)",
     "links": []},
]

BRIDGES = [
    {"slug": "FrontCover", "title": "Two curves, two lattices",
     "blurb": "640 points of y² = x³ + 3x (mod 5) over 𝔽₆₂₅, and 2379 points "
              "of y² = x³ + 3 (mod 7) over 𝔽₂₄₀₁; each marked point "
              "highlighted.",
     "links": [("the square curve", "/Explorer?a=-4&p=5"),
               ("the hexagonal curve", "/Explorer?a=-5&p=7")]},
    {"slug": "4_Torus", "title": "y² = x³ + 3x (mod 5), embedded",
     "blurb": "The 𝔽₅ and 𝔽₂₅ points on the square torus; the cutaway shows "
              "the Cayley graph wrapping the surface.",
     "links": [("open in the Explorer", "/Explorer?a=-4&p=5")]},
    {"slug": "4_Gallery", "title": "…and up the tower",
     "blurb": "The same curve over 𝔽₁₂₅, 𝔽₆₂₅ and 𝔽₃₁₂₅ — fixed points of "
              "successive powers of the lifted Frobenius.",
     "links": [("open in the Explorer", "/Explorer?a=-4&p=5")]},
    {"slug": "3_Intro", "title": "The hexagonal torus",
     "blurb": "Two conformal embeddings of ℂ/(ℤ ⊕ ωℤ), and the 39 points of "
              "y² = x³ + 3 (mod 7) over 𝔽₄₉.",
     "links": [("open in the Explorer", "/Explorer?a=-5&p=7")]},
    {"slug": "3_Gallery", "title": "y² = x³ + 3 (mod 7)",
     "blurb": "Over 𝔽₃₄₃, 𝔽₂₄₀₁ and 𝔽₁₆₈₀₇ — the sixfold symmetry of ℤ[ω], "
              "invisible in coordinates, unmissable on the torus.",
     "links": [("open in the Explorer", "/Explorer?a=-5&p=7"),
               ("disc −3", "/Explorer?d=-3")]},
    {"slug": "7_Gallery", "title": "y² = x³ + 5x + 7 (mod 11)",
     "blurb": "Over 𝔽₁₂₁, 𝔽₁₃₃₁ and 𝔽₁₄₆₄₁. The Frobenius order has "
              "discriminant −28; its endomorphism ring sits one level up, at "
              "field discriminant −7.",
     "links": [("open in the Explorer", "/Explorer?a=-4&p=11"),
               ("disc −28", "/Explorer?d=-28")]},
    {"slug": "8_Gallery", "title": "y² = x³ + x + 3 (mod 11)",
     "blurb": "Over 𝔽₁₂₁, 𝔽₁₃₃₁ and 𝔽₁₄₆₄₁ — the ℤ[√−2] lattice.",
     "links": [("open in the Explorer", "/Explorer?a=-6&p=11"),
               ("disc −8", "/Explorer?d=-8")]},
    {"slug": "11_Gallery", "title": "y² = x³ + x + 1 (mod 5)",
     "blurb": "Over 𝔽₁₂₅, 𝔽₆₂₅ and 𝔽₃₁₂₅. Class number 1: a single lattice "
              "class carries the whole story.",
     "links": [("open in the Explorer", "/Explorer?a=-3&p=5"),
               ("disc −11", "/Explorer?d=-11")]},
    {"slug": "RealLifts", "title": "Real curves, lifted",
     "blurb": "y² = x³ + 3x, y² = x³ + 1 and y² = x³ − x in the plane and on "
              "the uniformized curve — the real locus is a circle (or two) on "
              "the torus, and infinity is just a point.",
     "links": []},
    {"slug": "hopf", "title": "The Hopf fibration",
     "blurb": "Point preimages are pairwise-linked circles; the preimage of a "
              "curve on S² is a flat torus in S³.",
     "links": [("how the embedding works — Background §3", "/Background")]},
    {"slug": "HopfTori", "title": "Hopf tori",
     "blurb": "Curves on S² and the flat tori they bound — the machine behind "
              "every embedding in this gallery.",
     "links": [("how the embedding works — Background §3", "/Background")]},
]

GROUPS = [
    ("From the ICM gallery",
     "One wall for a single prime — every isogeny class over 𝔽₅ — and one "
     "wall for a single equation, X₀(11), across many characteristics.",
     ICM),
    ("From the Bridges paper",
     "The renders of [*Elliptic Curves and the Hopf Fibration*]"
     "(https://arxiv.org/abs/2505.09627) (Bridges 2025).",
     BRIDGES),
]

ALL = {e["slug"]: e for e in ICM + BRIDGES}


# ── grid rendering ────────────────────────────────────────────────────────────

_GRID_CSS = """
<style>
.gal-grid { display:grid; grid-template-columns:repeat(auto-fill, minmax(165px, 1fr));
            gap:12px; margin:6px 0 18px 0; }
.gal-card, .gal-card:hover, .gal-card:visited {
            display:block; border-radius:10px; overflow:hidden;
            text-decoration:none !important; color:inherit !important;
            border:1px solid rgba(128,138,148,0.25); background:rgba(128,138,148,0.06); }
.gal-card img { width:100%; aspect-ratio:1/1; object-fit:cover; display:block;
                transition:transform .15s ease; }
.gal-card:hover img { transform:scale(1.045); }
.gal-cap { padding:6px 9px 7px 9px; font-size:.83rem; line-height:1.3; opacity:.9; }
</style>
"""


@st.cache_data(show_spinner=False)
def _grid_html(slugs: tuple) -> str:
    cards = []
    for slug in slugs:
        e = ALL[slug]
        b64 = base64.b64encode((IMG / "thumb" / f"{slug}.jpg").read_bytes()).decode()
        cards.append(
            f'<a class="gal-card" href="?pic={slug}" title="{e["blurb"]}">'
            f'<img src="data:image/jpeg;base64,{b64}" alt="{e["title"]}">'
            f'<div class="gal-cap">{e["title"]}</div></a>'
        )
    return '<div class="gal-grid">' + "".join(cards) + "</div>"


def wall_view():
    st.header("Gallery")
    st.markdown(
        "Renders by Nadir Hajouji & [Steve Trettel](https://stevejtrettel.site/) "
        "(one guest piece by Elliot Kienzle) — elliptic curves over finite "
        "fields, lifted to characteristic zero and embedded as tori in ℝ³. "
        "Click any picture; more at "
        "[elliptic-curves.art](https://elliptic-curves.art/)."
    )
    st.markdown(_GRID_CSS, unsafe_allow_html=True)
    for title, sub, entries in GROUPS:
        st.markdown(f"### {title}")
        st.markdown(sub)
        st.markdown(_grid_html(tuple(e["slug"] for e in entries)),
                    unsafe_allow_html=True)


def detail_view(slug: str):
    e = ALL[slug]
    st.markdown("[← back to the Gallery](?)")
    st.subheader(e["title"])
    gif = IMG / "full" / f"{slug}.gif"
    jpg = IMG / "full" / f"{slug}.jpg"
    if gif.exists():   # animated (e.g. the lenticular flip) — embed so it animates
        b64 = base64.b64encode(gif.read_bytes()).decode()
        st.markdown(
            f'<img src="data:image/gif;base64,{b64}" '
            f'style="width:100%;border-radius:8px;" alt="{e["title"]}">',
            unsafe_allow_html=True)
    elif jpg.exists():
        st.image(str(jpg), width="stretch")
    else:
        st.info("The render for this piece is coming soon.")
    st.markdown(e["blurb"])
    if e.get("desc"):
        st.markdown(e["desc"])
    if e["links"]:
        st.markdown(" · ".join(f"[{lbl}]({href})" for lbl, href in e["links"]))
    st.caption(f"Render: {e.get('credit', CREDIT_DEFAULT)}")


# ── router ────────────────────────────────────────────────────────────────────
# ?pic=<slug> is the canonical detail link; ?g=<id> is the short link printed on
# the ICM gallery QR cards (pycode/shortlinks.py), resolved here to a slug.
_pic = st.query_params.get("pic")
if not _pic:
    _g = st.query_params.get("g")
    if _g:
        from shortlinks import resolve
        _pic = resolve(_g)
if _pic in ALL:
    detail_view(_pic)
else:
    wall_view()
