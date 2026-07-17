#!/usr/bin/env python3
"""Offline high-quality render of the Frobenius-flow animation as a looping GIF.

Ports the maths of ``basics_viz.frobenius_flow_html`` (the Background-page applet)
and renders it with matplotlib at high resolution -- antialiased dots, comet
trails, and the angle colour-wheel -- for the pinned curve

    y^2 = x^3 + 3x  over  F_5,   Lambda = Z[i],   Frobenius = mult by alpha = -2 + i.

Each colour rides its own Frobenius orbit and returns home after (its degree)
beats, so the whole picture repeats after PERIOD = lcm(degrees) beats: 1 for F_5,
2 for F_25, 6 for F_125. We render exactly one PERIOD, so the GIF loops perfectly.

    python3 scripts/frobenius_flow_gif.py --view galaxy --level 3 --out frob_galaxy.mp4

Output format is chosen by the --out extension: .mp4 / .webm (full colour, crisp,
small -- needs imageio + imageio-ffmpeg) or .gif (256 colours, embeds anywhere --
needs only Pillow). Frames are drawn with matplotlib and supersampled (--ss) for
antialiasing. Self-contained: no Streamlit stack.
"""
import argparse
import colorsys
import math

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Circle
from PIL import Image

# ---- the curve / lattice -----------------------------------------------------
ALPHA = (-2, 1)                      # Frobenius multiplier, alpha = -2 + i
THETA = math.atan2(ALPHA[1], ALPHA[0])
RHO = math.hypot(*ALPHA)             # |alpha| = sqrt 5
TAUP = (0.4, 0.2)                    # rebased galaxy basis tau' = 0.4 + 0.2 i
DEN = (-2, 1)                        # c*tau0 + d for the galaxy transport
QABS = math.exp(-2 * math.pi * TAUP[1])
PANEL_BG = "#17191c"


def amul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def apow(phi):
    r = RHO ** phi
    return (r * math.cos(THETA * phi), r * math.sin(THETA * phi))


def r01(x):
    return x - math.floor(x)


def r05(x):
    return x - round(x)


def fixed_pts(b):
    br, bi = b
    N = br * br + bi * bi
    seen, out = set(), []
    for u in range(N):
        for v in range(N):
            x = (u * br + v * bi) % N
            y = (v * br - u * bi) % N
            if (x, y) not in seen:
                seen.add((x, y))
                out.append((x / N, y / N))
    return out


def key_of(q):
    return (round(q[0] * 720000), round(q[1] * 720000))


def build(level):
    """Return pts, deg, hue, perm, PERIOD, base_tab -- matching the applet."""
    pts, deg, idx = [], [], {}
    an = (1, 0)
    for n in range(1, level + 1):
        an = amul(an, ALPHA)
        for q in fixed_pts((an[0] - 1, an[1])):
            k = key_of(q)
            if k not in idx:
                idx[k] = len(pts)
                pts.append(q)
                deg.append(n)
    hue = []
    for (u, v) in pts:
        uc, vc = r05(u), r05(v)
        hue.append(-1.0 if (uc == 0 and vc == 0)
                   else (math.degrees(math.atan2(vc, uc)) + 360) % 360)
    perm = []
    for q in pts:
        w = amul(ALPHA, q)
        perm.append(idx[key_of((r01(w[0]), r01(w[1])))])

    def lcm(a, b):
        return a // math.gcd(a, b) * b

    period = 1
    for d in set(deg):
        period = lcm(period, d)
    base_tab = [list(range(len(pts)))]
    for _ in range(1, period):
        base_tab.append([perm[j] for j in base_tab[-1]])
    return pts, deg, hue, perm, period, base_tab


# ---- geometry (460-unit coordinate system, same as the canvas) ---------------
W = 460.0
CX = CY = W / 2
MG = 40.0
S = W - 2 * MG          # torus half-cell scale
R = W / 2 - MG          # galaxy radius


def to_torus(z):
    return (CX + r05(z[0]) * S, CY - r05(z[1]) * S)


def to_galaxy(z):
    zr, zi = z
    dr, di = DEN
    n2 = dr * dr + di * di
    wr = (zr * dr + zi * di) / n2
    wi = (zi * dr - zr * di) / n2
    vs = wi / TAUP[1]                      # unreduced (1, tau') coords of z'
    vp, up = r01(vs), r01(wr - vs * TAUP[0])   # reduce each independently
    rr = R * math.exp(-2 * math.pi * vp * TAUP[1])
    th = 2 * math.pi * (up + vp * TAUP[0])
    return (CX + rr * math.cos(th), CY - rr * math.sin(th))


def strip_of(z, view):
    if view == "torus":
        return (z[0], z[1])
    zr, zi = z
    dr, di = DEN
    n2 = dr * dr + di * di
    wr = (zr * dr + zi * di) / n2
    wi = (zi * dr - zr * di) / n2
    vs = wi / TAUP[1]
    return (wr - vs * TAUP[0], vs)


def head_red(u, v, view):          # reduce the unreduced strip head so it lands on the dot
    if view == "torus":
        return (r05(u), r05(v))
    return (r01(u), r01(v))


def strip_to_screen(u, v, view):
    if view == "torus":
        return (CX + u * S, CY - v * S)
    rr = R * math.exp(-2 * math.pi * v * TAUP[1])
    th = 2 * math.pi * (u + v * TAUP[0])
    return (CX + rr * math.cos(th), CY - rr * math.sin(th))


def project(z, view):
    return to_torus(z) if view == "torus" else to_galaxy(z)


def hsl(h, s, l, a=1.0):
    if h < 0:                                   # identity marker
        return (0.92, 0.92, 0.92, a)
    r, g, b = colorsys.hls_to_rgb((h % 360) / 360.0, l, s)
    return (r, g, b, a)


# ---- one frame ---------------------------------------------------------------
TRAILMAX = 52.0
TRAILSTEP = 0.010
TRAILSTEPS = 22


def draw_frame(ax, tau, data, view, sc):
    pts, deg, hue, perm, period, base_tab = data
    ax.clear()
    ax.set_xlim(0, W)
    ax.set_ylim(0, W)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_facecolor(PANEL_BG)

    beat = math.floor(tau)
    bt = base_tab[beat % period]
    phi = tau - beat
    af = apow(phi)

    # frame guides
    if view == "torus":
        for dx, dy in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            ax.add_patch(plt.Rectangle((CX + (dx - .5) * S, CY - (dy + .5) * S), S, S,
                                       fill=False, ec=(1, 1, 1, .10), lw=.8 * sc))
        ax.add_patch(plt.Rectangle((CX - .5 * S, CY - .5 * S), S, S,
                                   fill=False, ec=(1, 1, 1, .26), lw=1.1 * sc))
    else:
        for rad, al in ((R, .18), (R * QABS, .18)):
            ax.add_patch(Circle((CX, CY), rad, fill=False, ec=(1, 1, 1, al), lw=1.0 * sc))
        vv = np.linspace(0, 1, 200)
        rr = R * np.exp(-2 * np.pi * vv * TAUP[1])
        th = 2 * np.pi * vv * TAUP[0]
        ax.plot(CX + rr * np.cos(th), CY - rr * np.sin(th), color=(1, 1, 1, .12), lw=1.0 * sc)

    # shadows: every resting spot, coloured by its home point's hue
    ax.add_collection(PatchCollection(
        [Circle(project(pts[i], view), 3.4) for i in range(len(pts))],
        facecolors=[hsl(hue[i], .52, .52, .30) for i in range(len(pts))],
        edgecolors="none", match_original=False))

    # trails: short comet tails, capped by screen length, fading to nothing
    segs, cols, lws = [], [], []
    for i in range(len(pts)):
        base = pts[bt[i]]
        sh = strip_of(amul(af, base), view)
        hr = head_red(sh[0], sh[1], view)
        head = strip_to_screen(hr[0], hr[1], view)
        path, acc, prev = [head], 0.0, head
        for k in range(1, TRAILSTEPS + 1):
            if acc >= TRAILMAX:
                break
            ss = strip_of(amul(apow(phi - k * TRAILSTEP), base), view)
            cur = strip_to_screen(hr[0] + (ss[0] - sh[0]), hr[1] + (ss[1] - sh[1]), view)
            acc += math.hypot(cur[0] - prev[0], cur[1] - prev[1])
            path.append(cur)
            prev = cur
        n = len(path)
        for k in range(1, n):
            f = 1 - (k - 1) / n
            segs.append([path[k - 1], path[k]])
            cols.append(hsl(hue[i], .68, .60, .5 * f * f))
            lws.append((0.6 + 1.7 * f) * sc)
    if segs:
        ax.add_collection(LineCollection(segs, colors=cols, linewidths=lws,
                                         capstyle="round"))

    # dots: each walks its Frobenius orbit
    circ, cols = [], []
    sel = None
    for i in range(len(pts)):
        x, y = project(amul(af, pts[bt[i]]), view)
        circ.append(Circle((x, y), 6.0 if hue[i] < 0 else 4.6))
        cols.append(hsl(hue[i], .68, .58, .98))
    ax.add_collection(PatchCollection(circ, facecolors=cols, edgecolors="none"))


def render(view, level, out, px, fps, fpb, ss=2, crf=19):
    data = build(level)
    period = data[4]
    sc = px * ss / W                             # size scale vs the 460 canvas (supersampled)
    nframes = fpb * period
    dpi = 100
    big = px * ss
    fig = plt.figure(figsize=(big / dpi, big / dpi), dpi=dpi)
    fig.patch.set_facecolor(PANEL_BG)
    ax = fig.add_axes([0, 0, 1, 1])

    def frame_img(f):                            # render frame f -> downsampled RGB PIL image
        draw_frame(ax, f / nframes * period, data, view, sc)
        fig.canvas.draw()
        im = Image.fromarray(np.asarray(fig.canvas.buffer_rgba()), "RGBA").convert("RGB")
        return im.resize((px, px), Image.LANCZOS) if ss != 1 else im

    ext = out.rsplit(".", 1)[-1].lower()
    if ext == "gif":                             # PIL needs all frames in hand
        frames = [frame_img(f) for f in range(nframes)]
        pal = frames[0].convert("P", palette=Image.ADAPTIVE, colors=256)
        q = [im.quantize(palette=pal, dither=Image.NONE) for im in frames]
        q[0].save(out, save_all=True, append_images=q[1:], loop=0,
                  duration=round(1000 / fps), optimize=True, disposal=2)
    else:                                        # stream frames to ffmpeg (constant memory)
        import imageio.v2 as imageio
        if ext == "webm":
            codec, params = "libvpx-vp9", ["-crf", str(crf + 10), "-b:v", "0"]
        else:                                    # H.264: higher crf = smaller file
            codec, params = "libx264", ["-crf", str(crf), "-preset", "slow",
                                        "-pix_fmt", "yuv420p"]
        writer = imageio.get_writer(out, fps=fps, codec=codec, macro_block_size=None,
                                    ffmpeg_params=params)
        for f in range(nframes):
            writer.append_data(np.asarray(frame_img(f)))
        writer.close()
    plt.close(fig)
    print(f"wrote {out}: {nframes} frames, {px}px (ss{ss}), {nframes/fps:.0f}s, "
          f"{fpb/fps:.1f}s/beat, period {period} beats")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--view", choices=["torus", "galaxy"], default="galaxy")
    ap.add_argument("--level", type=int, default=3, help="1=F5, 2=F25, 3=F125")
    ap.add_argument("--out", default=None, help="extension picks format: .mp4 / .webm / .gif")
    ap.add_argument("--px", type=int, default=800)
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--fpb", type=int, default=24, help="frames per beat (ignored if --beat-sec set)")
    ap.add_argument("--beat-sec", type=float, default=None,
                    help="seconds of playback per Frobenius beat (relaxing ~= 30); overrides --fpb")
    ap.add_argument("--ss", type=int, default=2, help="supersample factor for antialiasing")
    ap.add_argument("--crf", type=int, default=19,
                    help="H.264 quality: lower=better/bigger, higher=smaller (~18 great, ~26 small)")
    a = ap.parse_args()
    out = a.out or f"frob_{a.view}_L{a.level}.mp4"
    fpb = round(a.fps * a.beat_sec) if a.beat_sec else a.fpb
    render(a.view, a.level, out, a.px, a.fps, fpb, a.ss, a.crf)


if __name__ == "__main__":
    main()
