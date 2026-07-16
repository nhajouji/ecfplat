"""Short-link ids for the printed QR codes on the ICM gallery cards.

Each printed QR encodes   https://elliptic-curves.info/Gallery?g=<id>
which the gallery router (pages/5_Gallery.py) resolves to  ?pic=<slug>.

The <id>s are the contract with the *printed, physical* cards: once a card
is printed, its id can never change — but the slug it maps to can, and the
whole page can move, without reprinting anything. Keep ids short (denser,
more scannable QR) and STABLE.
"""

# short id  ->  gallery detail slug (key in pages/5_Gallery.py ALL)
LINKS = {
    # Wall 1 — every isogeny class over F5
    "a1":  "icm-a1",         # y^2 = x^3 + 3x + 2   (mod 5)
    "a2":  "icm-a2",         # y^2 = x^3 + x  /  + x + 2   (the column)
    "a3":  "icm-a3",         # y^2 = x^3 + 4x + 2   (mod 5)
    "a4":  "icm-a4",         # y^2 = x^3 + 2x   (mod 5)
    "a0":  "icm-a0",         # y^2 = x^3 + 1   (mod 5)  supersingular / lenticular
    # Wall 2 — X_0(11) across characteristics
    "xc":   "icm-pcomplex",  # X_0(11) over C
    "p23":  "icm-p23",       # X_0(11) mod 23
    "p101": "icm-p101",      # X_0(11) mod 101
    "p107": "icm-p107",      # X_0(11) mod 107
    "x05":  "icm-x05",       # real elliptic curves over X_0(5)
    "star": "icm-starscape", # guest piece (Elliot Kienzle)
}


def resolve(sid: str) -> str | None:
    """Return the gallery slug for a short id, or None if unknown."""
    return LINKS.get((sid or "").strip().lower())
