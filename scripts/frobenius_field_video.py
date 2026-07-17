#!/usr/bin/env python3
"""Offline export of the big Frobenius-flow field E(F_5^m) as a looping MP4/GIF.

Matches the WebGL point-sprite prototype (shaders/frobenius_galaxy_F5_6.html):
all of E(F_5^m) = ker(alpha^m - 1) for the curve y^2 = x^3 + 3x over F_5
(alpha = -2 + i), drawn in the galaxy model, coloured by angle, with point size
set by field degree (simpler subfield = bigger).  Colour is position-based, so
the picture returns to itself after ONE Frobenius beat -> a short seamless loop.

    python3 scripts/frobenius_field_video.py --m 6 --out frob_F5_6.mp4

m = 6 gives E(F_5^6) = 15,860 points (degrees 1,2,3,6: the subfields F5, F25,
F125 inside their compositum).  numpy + matplotlib + Pillow (+ imageio for mp4).
"""
import argparse, colorsys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
from PIL import Image

ALPHA = (-2.0, 1.0)
ARGA  = math.atan2(1.0, -2.0)
ABSA  = math.sqrt(5.0)
TAUP  = (0.4, 0.2)
DEN   = (-2.0, 1.0)                 # transport denominator = alpha
QABS  = math.exp(-2 * math.pi * TAUP[1])
W = 460.0; CX = CY = W / 2; R = W / 2 - 40; BG = (0.055, 0.06, 0.067)

def amul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])

# ---- build E(F_5^m) = ker(alpha^m - 1): points (u,v) on the torus + degree ----
def build_field(m):
    an = (1.0, 0.0)
    for _ in range(m): an = amul(an, ALPHA)
    beta = (an[0]-1.0, an[1]); N = round(beta[0]**2 + beta[1]**2)
    cr, ci = beta[0], -beta[1]                       # conj(beta)
    fr = lambda x: x - math.floor(x)
    seen, us, vs = set(), [], []
    B = int(math.isqrt(N)) + 2
    for x in range(-B, B+1):
        for y in range(-B, B+1):
            u = fr((x*cr - y*ci)/N); v = fr((x*ci + y*cr)/N)
            key = (round(u*N), round(v*N))
            if key not in seen:
                seen.add(key); us.append(u); vs.append(v)
    divs = [d for d in range(1, m+1) if m % d == 0]
    betad = {}
    for d in divs:
        ad = (1.0, 0.0)
        for _ in range(d): ad = amul(ad, ALPHA)
        betad[d] = (ad[0]-1.0, ad[1])
    def degof(u, v):
        for d in divs:
            b = betad[d]
            if abs(b[0]*u-b[1]*v - round(b[0]*u-b[1]*v)) < 1e-6 and \
               abs(b[0]*v+b[1]*u - round(b[0]*v+b[1]*u)) < 1e-6:
                return d
        return m
    deg = np.array([degof(u, v) for u, v in zip(us, vs)])
    return np.array(us), np.array(vs), deg, N

# ---- galaxy transform (vectorised): torus point q -> screen (X,Y) and angle ----
def galaxy(qx, qy):
    dr, di = DEN; n2 = dr*dr + di*di
    wr = (qx*dr + qy*di)/n2; wi = (qy*dr - qx*di)/n2
    vs = wi/TAUP[1]; vp = vs - np.floor(vs)
    up = wr - vs*TAUP[0]; up = up - np.floor(up)
    rr = R*np.exp(-2*np.pi*vp*TAUP[1]); th = 2*np.pi*(up + vp*TAUP[0])
    return CX + rr*np.cos(th), CY - rr*np.sin(th), th

def render(m, out, px, fps, frames, ss, crf):
    gx, gy, deg, N = build_field(m)
    print(f"E(F_5^{m}) = {N} points; degree counts "
          f"{ {int(d): int((deg==d).sum()) for d in sorted(set(deg))} }")
    size_map = {1: 42.0, 2: 26.0, 3: 14.0, 6: 2.6, 4: 8.0, 5: 4.0}   # area; smaller subfield = bigger
    base = np.array([size_map.get(int(d), 3.0) for d in deg]) * ss * ss
    alpha = np.where(deg <= 3, 1.0, 0.6)
    order = np.argsort(-deg)                          # draw dust first, subfield dots on top

    dpi = 100; big = px*ss
    fig = plt.figure(figsize=(big/dpi, big/dpi), dpi=dpi); fig.patch.set_facecolor(BG)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, W)
    ax.set_aspect("equal"); ax.axis("off")

    def frame_img(f):
        phi = f/frames
        r = ABSA**phi; t = ARGA*phi; afx, afy = r*math.cos(t), r*math.sin(t)
        qx = afx*gx - afy*gy; qy = afx*gy + afy*gx     # q = alpha^phi * g
        X, Y, th = galaxy(qx, qy)
        rgb = hsv_to_rgb(np.stack([(th/(2*np.pi)) % 1.0,
                                   np.full_like(th, 0.72), np.full_like(th, 0.98)], -1))
        ax.clear(); ax.set_xlim(0, W); ax.set_ylim(0, W); ax.axis("off")
        ax.scatter(X[order], Y[order], s=base[order], c=rgb[order],
                   alpha=alpha[order], linewidths=0, edgecolors="none")
        fig.canvas.draw()
        im = Image.fromarray(np.asarray(fig.canvas.buffer_rgba()), "RGBA").convert("RGB")
        return im.resize((px, px), Image.LANCZOS) if ss != 1 else im

    ext = out.rsplit(".", 1)[-1].lower()
    if ext == "gif":
        fr = [frame_img(f) for f in range(frames)]
        pal = fr[0].convert("P", palette=Image.ADAPTIVE, colors=256)
        q = [im.quantize(palette=pal, dither=Image.NONE) for im in fr]
        q[0].save(out, save_all=True, append_images=q[1:], loop=0,
                  duration=round(1000/fps), optimize=True, disposal=2)
    else:
        import imageio.v2 as imageio
        wtr = imageio.get_writer(out, fps=fps, codec="libx264", macro_block_size=None,
                                 ffmpeg_params=["-crf", str(crf), "-preset", "slow",
                                                "-pix_fmt", "yuv420p"])
        for f in range(frames):
            wtr.append_data(np.asarray(frame_img(f)))
        wtr.close()
    plt.close(fig)
    print(f"wrote {out}: {frames} frames, {px}px (ss{ss}), {frames/fps:.1f}s loop")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, default=6, help="field is F_5^m")
    ap.add_argument("--out", default="frob_F5_6.mp4", help=".mp4 / .gif")
    ap.add_argument("--px", type=int, default=900)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--frames", type=int, default=360, help="frames for one beat (loop length = frames/fps s)")
    ap.add_argument("--ss", type=int, default=2)
    ap.add_argument("--crf", type=int, default=20)
    a = ap.parse_args()
    render(a.m, a.out, a.px, a.fps, a.frames, a.ss, a.crf)

if __name__ == "__main__":
    main()
