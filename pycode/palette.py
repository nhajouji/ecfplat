"""
palette.py
==========
Shared categorical colours keyed by lattice-class row index.

Both the Plotly fundamental-domain picture (`plotly_tools`) and the matplotlib
isogeny graph (`graphic_tools`) colour classes with `row_colors(n)`, so a given
class has the *same* colour everywhere it appears. The colour for a class is its
position in the dataframe / `qfs_ordered` (row index), not its discriminant.

No heavy dependencies (stdlib `colorsys` only) so it is safe to import from both
the matplotlib and Plotly sides.
"""

import colorsys

# Hue band to skip: yellows read poorly on a white background.
_YELLOW_GAP = (0.13, 0.20)


def row_colors(n: int, saturation: float = 0.62, lightness: float = 0.46):
    """`n` visually distinct hex colours, evenly spread around the hue wheel.

    The yellow band is skipped and lightness is held at a mid-dark value, so
    every colour stays legible on white (and takes white text well).
    """
    if n <= 0:
        return []
    gap_lo, gap_hi = _YELLOW_GAP
    gap = gap_hi - gap_lo
    span = 1.0 - gap
    out = []
    for i in range(n):
        h = (i / n) * span
        if h >= gap_lo:            # hop over the yellow band
            h += gap
        r, g, b = colorsys.hls_to_rgb(h, lightness, saturation)
        out.append("#{:02x}{:02x}{:02x}".format(
            int(r * 255 + 0.5), int(g * 255 + 0.5), int(b * 255 + 0.5)))
    return out


def text_color_for(hex_color: str) -> str:
    """Return black or white, whichever reads better on `hex_color`."""
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255.0
    return "#000000" if luminance > 0.6 else "#ffffff"
