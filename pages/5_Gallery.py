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
     "eqn": "y² = x³ + 3x + 2 &nbsp;(mod 5)",
     "blurb": "The isogeny class (a, p) = (1, 5): a single curve with 5 "
              "rational points and CM by discriminant −19.",
     "desc": (
        r"Exactly $p$ points: $\#E(\mathbb{F}_5) = 5$. The anomalous class — CM "
        r"by $\mathbb{Q}(\sqrt{-19})$, $j = (-96)^3$."),
     "idcard": {"cap": "isogeny class (1, 5)", "fd": "a1",
                "fd_cap": "the CM point — on the wall Re = ½",
                "rows": [("Frobenius π", "(1 + √−19)/2"),
                         ("⇒ implies", "𝔽₅ · 5 points · trace 1"),
                         ("endomorphisms", "ℤ[(1+√−19)/2] · h = 1"),
                         ("j (char 0)", "−884736 = (−96)³ → 4 (mod 5)")]},
     "seed": ["CM by ℚ(√−19); class number 1 — a single lattice, no hidden "
              "siblings.",
              "Frobenius root is a free choice: either of (1 ± √−19)/2.",
              "Lift j = −884736 = (−96)³; Weber's γ₂ = ∛j is the integer −96."],
     "links": [("open in the Explorer", "/Explorer?a=1&p=5"),
               ("disc −19", "/Explorer?d=-19")]},
    {"slug": "icm-a2", "title": "Trace 2 over 𝔽₅ — the column",
     "eqn": ("<div>y² = x³ + x <span class='ab'>(above)</span></div>"
             "<div>y² = x³ + x + 2 <span class='ab'>(below)</span></div>"
             "<div class='modline'>(mod 5)</div>"),
     "blurb": "Two curves joined by a 2-isogeny: j = 1728 above (disc −4), "
              "its descendant below (disc −16). Hung vertically so ascending "
              "and descending mean what they say.",
     "desc": (
        r"Two curves, four points each over $\mathbb{F}_5$ — the same count, told "
        r"apart only by their $2$-torsion. Stacked to trace the volcano."),
     "idcard": {"cap": "isogeny class (2, 5)", "fd": "a2",
                "fd_cap": "i (above) & 2i (below)",
                "rows": [("Frobenius π", "1 + 2i"),
                         ("⇒ implies", "𝔽₅ · 4 points · trace 2"),
                         ("above", "ℤ[i] (−4) · j = 1728 = 12³ → 3", "two"),
                         ("below", "ℤ[2i] (−16) · j = 287496 = 66³ → 1", "two b")]},
     "seed": ["Two curves, not two lifts: the Frobenius order (disc −16) is "
              "non-maximal — a 2-isogeny volcano with the maximal-order curve "
              "(disc −4, j = 1728) on top and its descendant (disc −16) below.",
              "Frobenius 1 + 2i is an associate of the trace-4 curve's 2 − i:  "
              "1 + 2i = i·(2 − i)."],
     "links": [("open in the Explorer", "/Explorer?a=2&p=5"),
               ("disc −16", "/Explorer?d=-16")]},
    {"slug": "icm-a3", "title": "Trace 3 over 𝔽₅",
     "eqn": "y² = x³ + 4x + 2 &nbsp;(mod 5)",
     "blurb": "One curve, three rational points, CM by discriminant −11.",
     "desc": (
        r"Three points over $\mathbb{F}_5$. Discriminant $-11$, class number "
        r"one: a single lattice, nothing to choose — the quiet opposite of the "
        r"six lifts down the wall."),
     "idcard": {"cap": "isogeny class (3, 5)", "fd": "a3",
                "fd_cap": "the CM point — on the wall Re = ½",
                "rows": [("Frobenius π", "(3 + √−11)/2"),
                         ("⇒ implies", "𝔽₅ · 3 points · trace 3"),
                         ("endomorphisms", "ℤ[(1+√−11)/2] · h = 1"),
                         ("j (char 0)", "−32768 = (−32)³ → 2 (mod 5)")]},
     "seed": ["CM by ℚ(√−11); class number 1.",
              "Three rational points (trace 3).",
              "Lift j = −32768 = (−32)³; γ₂ = ∛j = −32."],
     "links": [("open in the Explorer", "/Explorer?a=3&p=5"),
               ("disc −11", "/Explorer?d=-11")]},
    {"slug": "icm-a4", "title": "Trace 4 over 𝔽₅",
     "eqn": "y² = x³ + 2x &nbsp;(mod 5)",
     "blurb": "j = 1728 again — the ℤ[i] lattice — now with trace 4 and just "
              "two rational points.",
     "desc": (
        r"$j = 1728$ — the Gaussian lattice $\mathbb{Z}[i]$. Here "
        r"$y^2 = x^3 + 2x$; on the Bridges wall, $y^2 = x^3 + 3x$. Different "
        r"equations, one picture: quadratic twists, indistinguishable over "
        r"$\mathbb{F}_{5^6}$."),
     "idcard": {"cap": "isogeny class (4, 5)", "fd": "a4",
                "fd_cap": "the point i  (j = 1728)",
                "rows": [("Frobenius π", "2 − i"),
                         ("⇒ implies", "𝔽₅ · 2 points · trace 4"),
                         ("endomorphisms", "ℤ[i] (−4) · h = 1"),
                         ("j (char 0)", "1728 = 12³ → 3 (mod 5)")]},
     "seed": ["j = 1728: the ℤ[i] curve, with extra automorphisms (units ±1, "
              "±i).",
              "Frobenius 2 − i is an associate of the trace-2 curve's 1 + 2i:  "
              "2 − i = −i·(1 + 2i).",
              "Just two rational points — the leanest class on the wall."],
     "links": [("open in the Explorer", "/Explorer?a=4&p=5"),
               ("disc −4", "/Explorer?d=-4")]},
    {"slug": "icm-a0", "title": "Trace 0 over 𝔽₅ — supersingular",
     "eqn": "y² = x³ + 1 &nbsp;(mod 5)",
     "blurb": "The supersingular class (a, p) = (0, 5): y² = x³ + 1, j = 0 — the "
              "one exceptional class over 𝔽₅, made as a lenticular print that "
              "flips between the two lifts as you move past it.",
     "desc": (
        r"One supersingular $j$-invariant over $\mathbb{F}_5$ — but two curves "
        r"realize it, and two lattices lift it. A lenticular print, so both lifts "
        r"of Frobenius share one frame."),
     "idcard": {"cap": "supersingular class (0, 5)", "fd": "a0",
                "fd_cap": "two CM points = the two lifts",
                "rows": [("Frobenius π", "±√−5 &nbsp;<span class='note'>(two "
                          "orientations)</span>"),
                         ("⇒ implies", "𝔽₅ · 6 points · trace 0"),
                         ("endomorphisms", "order disc −20 · h = 2"),
                         ("j (char 0)", "632000 ± 282880√5 → 0 (mod 5)")]},
     "seed": ["The only supersingular class over 𝔽₅ (trace 0).",
              "CM order disc −20 = ℚ(√−5) has class number 2 — two lifts, "
              "j = 632000 ± 282880√5, conjugate over ℚ.",
              "The two roots ±√−5 are the two orientations = the two lifts = "
              "the flip of the lenticular.",
              "Both reduce mod 5 to j ≡ 0: two faces in char 0, one "
              "supersingular curve in char 5 (Deuring lifting)."],
     "links": [("open in the Explorer", "/Explorer?a=0&p=5"),
               ("disc −20", "/Explorer?d=-20")]},
    {"slug": "icm-pcomplex", "title": "X₀(11) over ℂ",
     "blurb": "The complex points of X₀(11) form a torus; the real locus is "
              "traced in green, the marked point in red.",
     "links": [("modular curves in the Background", "/Background")]},
    {"slug": "icm-p23", "title": "X₀(11) mod 23",
     "blurb": "Reduced mod 23 the curve has trace −1 and discriminant −91 — "
              "class number 2, so Frobenius has two lifts.",
     "desc": (
        r"There are two possible lifts for $X_0(11) \pmod{23}$, and they look "
        r"very different."),
     "links": [("open in the Explorer", "/Explorer?a=-1&p=23"),
               ("disc −91", "/Explorer?d=-91")]},
    {"slug": "icm-p101", "title": "X₀(11) mod 101",
     "blurb": "Reduced mod 101, X₀(11) is ordinary; the piece shows the curve "
              "whose endomorphism ring is ℤ[2i], with its single CM lift "
              "ℂ/(ℤ ⊕ 2iℤ).",
     "desc": (
        r"The reduction $X_0(11) \pmod{101}$ is an ordinary curve with "
        r"endomorphism ring isomorphic to $\mathbb{Z}[2i]$. This means there is "
        r"only one possibility for the CM lift — the elliptic curve "
        r"$\mathbb{C}/(\mathbb{Z}\oplus 2i\,\mathbb{Z})$. Note that this same CM "
        r"lift is used in the picture for $(a,p)=(2,5)$ in the other section, "
        r"drawn with a different conformal model in $\mathbb{R}^3$."),
     "links": [("open in the Explorer", "/Explorer?a=2&p=101"),
               ("disc −400", "/Explorer?d=-400")]},
    {"slug": "icm-p107", "title": "X₀(11) mod 107 — the six lifts",
     "blurb": "Mod 107 the Frobenius order is ℤ[√−26], with class number 6: "
              "six lattices, six lifts, one curve. The class group is cyclic "
              "and the arrangement encodes it — the tori sit like the 6th "
              "roots of unity.",
     "desc": (
        r"There are $6$ possible CM lifts for $X_0(11) \pmod{107}$, "
        r"corresponding to the $6$ ideal classes in the class group of "
        r"$\mathrm{End}(X_0(11) \bmod 107) \cong \mathbb{Z}[\sqrt{-26}]$. The "
        r"arrangement of the CM lifts is meant to reflect the structure of the "
        r"underlying class group."),
     "links": [("open in the Explorer", "/Explorer?a=18&p=107"),
               ("disc −104", "/Explorer?d=-104")]},
    {"slug": "icm-x05", "title": "X₀(5) as a real surface",
     "blurb": "A slice of the universal curve over X₁(5) — every thread a "
              "real elliptic curve. Four fibers are highlighted: the two "
              "lattices with CM by √−5, and the pair X₁(11) ⇄ X₀(11), whose "
              "degree-5 map is itself a point of X₀(5).",
     "desc": (
        r"A picture depicting some of the real fibers on the universal curve "
        r"over a portion of $X_0(5)$. Two of the fibers represent the isogeny "
        r"$X_1(11)\to X_0(11)$ (the subject of our *One equation, many fields* "
        r"section), and the two solid-colored fibers represent the lifts of "
        r"Frobenius of the supersingular curve over $\mathbb{F}_5$ — the curves "
        r"in the lenticular print in the *One field, many equations* section. "
        r"This piece unifies the two sections of the gallery."),
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

GALLERIES = [
    {"header": "From the ICM gallery",
     "intro": r"Two walls — $\mathbb{F}_5$ held fixed while the equation varies, "
              r"and $X_0(11)$ held fixed while the field varies — bridged by a "
              r"single surface:",
     "lead": ["icm-x05"],
     "sections": [
        ("One field, many equations",
         r"The pieces here are the isogeny classes over $\mathbb{F}_5$. We use "
         r"the Deuring–Pinkall construction to draw each curve inside a CM lift "
         r"of the elliptic curve over $\mathbb{F}_5$; the points shown are those "
         r"defined over $\mathbb{F}_{5^6}$.",
         ["icm-a0", "icm-a1", "icm-a2", "icm-a3", "icm-a4"]),
        ("One equation, many fields",
         r"The well-known (elliptic) modular curve $X_0(11)$ from different "
         r"perspectives — some pieces depict the reductions $X_0(11) \bmod p$, "
         r"others the characteristic-zero object itself.",
         ["icm-pcomplex", "icm-p23", "icm-p101", "icm-p107", "icm-starscape"]),
     ]},
    {"header": "From the Bridges paper",
     "intro": "The renders of [*Elliptic Curves and the Hopf Fibration*]"
              "(https://arxiv.org/abs/2505.09627) (Bridges 2025).",
     "lead": [],
     "sections": [(None, None, [e["slug"] for e in BRIDGES])]},
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
    if "draft" in st.query_params:   # NH's private editing overlay (?draft=1)
        todo = [e for e in ICM if not e.get("desc")]
        links = " · ".join(f"[{e['title']}](?pic={e['slug']}&draft=1)" for e in todo)
        st.info(f"**DRAFT MODE** — {len(todo)} ICM page(s) still need prose: {links}")
    st.markdown(_GRID_CSS, unsafe_allow_html=True)
    for g in GALLERIES:
        st.markdown(f"## {g['header']}")
        if g.get("intro"):
            st.markdown(g["intro"])
        if g.get("lead"):
            st.markdown(_grid_html(tuple(g["lead"])), unsafe_allow_html=True)
        for title, sub, slugs in g["sections"]:
            if title:
                st.markdown(f"### {title}")
            if sub:
                st.markdown(sub)
            st.markdown(_grid_html(tuple(slugs)), unsafe_allow_html=True)


_DETAIL_CSS = """
<style>
.pg-eqn{font-family:Palatino,'Palatino Linotype',Georgia,serif;font-weight:600;
        font-size:1.7rem;line-height:1.3;margin:.2rem 0 .1rem}
.pg-eqn .ab{color:#8b93a0;font-size:.92rem;font-weight:400}
.pg-eqn .modline{color:#8b93a0;font-size:1.05rem;margin-top:.1rem}
.pg-sub{color:#8b93a0;font-size:.95rem;margin-bottom:.3rem}
.gc-draft{display:inline-block;background:rgba(233,170,70,.16);color:#e9b45a;
        border:1px solid rgba(233,170,70,.4);border-radius:20px;padding:2px 11px;
        font-size:.72rem;letter-spacing:.04em;font-weight:600;margin:.3rem 0}
.gc-seed{color:#c7ccd4;font-size:.92rem;line-height:1.5;margin:.2rem 0 .4rem}
.gc-seed ul{margin:.3rem 0 0;padding-left:1.1rem}.gc-seed li{margin:.34rem 0}
.gc-card{margin:14px 0;background:rgba(128,138,148,.06);max-width:560px;
        border:1px solid rgba(128,138,148,.28);border-radius:13px;padding:15px 18px}
.gc-cap{font-size:.7rem;letter-spacing:.08em;text-transform:uppercase;
        color:#7f8894;margin-bottom:10px}
.gc-row{display:flex;justify-content:space-between;gap:14px;padding:6px 0;
        border-top:1px solid rgba(128,138,148,.14);font-size:.92rem}
.gc-row:first-of-type{border-top:none}
.gc-k{color:#8b93a0}.gc-v{text-align:right;font-family:Palatino,Georgia,serif}
.gc-v .note{color:#7f8894;font-size:.8rem;font-family:-apple-system,sans-serif}
.gc-row.two .gc-k{color:#8fd0ff}.gc-row.two.b .gc-k{color:#c9a0ff}
.gc-row.gc-pi .gc-v{color:#f0d9a8;font-size:1.02rem}
.gc-row.gc-j .gc-v{color:#f0d9a8}
.gc-fd{display:flex;gap:14px;align-items:center;margin-top:12px}
.gc-fd img{height:132px;background:rgba(255,255,255,.015);border-radius:8px}
.gc-fd span{font-size:.82rem;color:#8b93a0}
.gc-exp{margin-top:9px;font-size:.9rem}
</style>
"""


def _idcard_html(e: dict) -> str:
    c = e["idcard"]
    rows = ""
    for row in c["rows"]:
        lab, val = row[0], row[1]
        cls = row[2] if len(row) > 2 else ""
        if lab.startswith("Frobenius"):
            cls += " gc-pi"
        elif lab.startswith("j ("):
            cls += " gc-j"
        rows += (f'<div class="gc-row {cls}"><span class="gc-k">{lab}</span>'
                 f'<span class="gc-v">{val}</span></div>')
    fd = base64.b64encode((IMG / "fd" / f'{c["fd"]}.png').read_bytes()).decode()
    links = " · ".join(f'<a href="{href}">{lbl}</a>' for lbl, href in e["links"])
    return (f'<div class="gc-card"><div class="gc-cap">ID · {c["cap"]}</div>{rows}'
            f'<div class="gc-fd"><img src="data:image/png;base64,{fd}">'
            f'<span>{c["fd_cap"]}</span></div>'
            f'<div class="gc-exp">→ {links}</div></div>')


def detail_view(slug: str):
    e = ALL[slug]
    draft = "draft" in st.query_params
    st.markdown("[← back to the Gallery](?)")
    st.markdown(_DETAIL_CSS, unsafe_allow_html=True)
    if e.get("eqn"):
        st.markdown(f'<div class="pg-eqn">{e["eqn"]}</div>'
                    f'<div class="pg-sub">{e["title"]}</div>', unsafe_allow_html=True)
    else:
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
    # draft scaffolding — visible only with ?draft, never to the public / at the ICM
    if draft and not e.get("desc"):
        st.markdown('<span class="gc-draft">DRAFT · prose to be written by NH</span>',
                    unsafe_allow_html=True)
        if e.get("seed"):
            lis = "".join(f"<li>{s}</li>" for s in e["seed"])
            st.markdown(f'<div class="gc-seed">Seed material:<ul>{lis}</ul></div>',
                        unsafe_allow_html=True)
    if e.get("idcard"):
        st.markdown(_idcard_html(e), unsafe_allow_html=True)
    elif e["links"]:
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
