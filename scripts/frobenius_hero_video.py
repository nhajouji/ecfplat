#!/usr/bin/env python3
"""Standalone export of the glowing homepage hero (pycode/hero_viz.py) to MP4/GIF.

Faithful reproduction of the canvas hero: the E(F_5^level) points of
y^2 = x^3 + 3x over F_5 (alpha = -2 + i) on the galaxy model, with ADDITIVE glow
(radial-gradient halos + comet tails), colour by angle, one seamless Frobenius
beat.  Glow is composited into a float buffer and clipped, matching the browser's
'lighter' blend -- something matplotlib's scatter can't do.

    python3 scripts/frobenius_hero_video.py --out frob_hero_F25.mp4

level 2 = F_25 (the 20-point hero).
"""
import argparse, math
import numpy as np
from PIL import Image

ALPHA = (-2.0, 1.0)
ARGA = math.atan2(1.0, -2.0)
ABSA = math.sqrt(5.0)
TAUP = (0.4, 0.2)
DEN = (-2.0, 1.0)
BG = np.array([0.055, 0.06, 0.067])

def amul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])

def gen_field(m):                                  # E(F_5^m) = ker(alpha^m - 1)
    an = (1.0, 0.0)
    for _ in range(m): an = amul(an, ALPHA)
    beta = (an[0]-1.0, an[1]); N = round(beta[0]**2 + beta[1]**2)
    cr, ci = beta[0], -beta[1]
    fr = lambda x: x - math.floor(x)
    seen, pts = set(), []
    B = int(math.isqrt(N)) + 2
    for x in range(-B, B+1):
        for y in range(-B, B+1):
            u = fr((x*cr - y*ci)/N); v = fr((x*ci + y*cr)/N)
            key = (round(u*N), round(v*N))
            if key not in seen:
                seen.add(key); pts.append((u, v))
    return pts

def apow(phi):
    r = ABSA**phi; t = ARGA*phi
    return (r*math.cos(t), r*math.sin(t))

def galaxy(q, CX, CY, R):                           # torus point -> screen (x,y) + angle
    zr, zi = q; dr, di = DEN; n2 = dr*dr + di*di
    wr = (zr*dr + zi*di)/n2; wi = (zi*dr - zr*di)/n2
    vs = wi/TAUP[1]; vp = vs - math.floor(vs); up = (wr - vs*TAUP[0]) % 1.0
    rr = R*math.exp(-2*math.pi*vp*TAUP[1]); th = 2*math.pi*(up + vp*TAUP[0])
    return CX + rr*math.cos(th), CY - rr*math.sin(th), th

def hue(a):                                         # angle -> rgb in [0,1] (matches hero_viz)
    h = (a/(2*math.pi)) % 1.0 * 6.0
    x = 1 - abs(h % 2 - 1)
    c = [[1,x,0],[x,1,0],[0,1,x],[0,x,1],[x,0,1],[1,0,x]][int(h) % 6]
    return np.array([0.28 + 0.72*v for v in c])

def splat(buf, x, y, rad, col, a):                  # additive radial-gradient sprite
    H, W = buf.shape[:2]
    x0, x1 = max(0, int(x-rad)), min(W, int(x+rad)+1)
    y0, y1 = max(0, int(y-rad)), min(H, int(y+rad)+1)
    if x1 <= x0 or y1 <= y0: return
    dx = np.arange(x0, x1) - x; dy = np.arange(y0, y1) - y
    d = np.hypot(dy[:, None], dx[None, :])
    g = np.clip(1 - d/rad, 0, 1) * a
    buf[y0:y1, x0:x1, :] += g[:, :, None] * col

def frame(G, phi, W, H, sc):
    CX, CY, R = W/2, H/2, H/2 - 40*sc
    buf = np.tile(BG, (H, W, 1)).astype(np.float64)
    af = apow(phi)
    for g in G:
        for j in (4, 3, 2, 1):                      # glowing comet tail
            x, y, th = galaxy(amul(apow(phi - j*0.011), g), CX, CY, R)
            f = 1 - j/5.0
            splat(buf, x, y, (5+7*f)*sc, hue(th), 0.16*f)
        x, y, th = galaxy(amul(af, g), CX, CY, R)   # head: two haloes + bright core
        col = hue(th)
        splat(buf, x, y, 34*sc, col, 0.20)
        splat(buf, x, y, 12*sc, col, 0.55)
        splat(buf, x, y, 4.2*sc, col, 0.95)
    return Image.fromarray((np.clip(buf, 0, 1)*255).astype(np.uint8), "RGB")

def render(level, out, W, H, fps, frames, crf):
    G = gen_field(level)
    print(f"E(F_5^{level}) = {len(G)} points")
    sc = H/400.0                                    # hero sizes were tuned at height 400
    ext = out.rsplit(".", 1)[-1].lower()
    if ext == "gif":
        fr = [frame(G, f/frames, W, H, sc) for f in range(frames)]
        pal = fr[0].convert("P", palette=Image.ADAPTIVE, colors=256)
        q = [im.quantize(palette=pal, dither=Image.NONE) for im in fr]
        q[0].save(out, save_all=True, append_images=q[1:], loop=0,
                  duration=round(1000/fps), optimize=True, disposal=2)
    else:
        import imageio.v2 as imageio
        w = imageio.get_writer(out, fps=fps, codec="libx264", macro_block_size=None,
                               ffmpeg_params=["-crf", str(crf), "-preset", "slow", "-pix_fmt", "yuv420p"])
        for f in range(frames):
            w.append_data(np.asarray(frame(G, f/frames, W, H, sc)))
        w.close()
    print(f"wrote {out}: {frames} frames, {W}x{H}, {frames/fps:.1f}s loop")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, default=2, help="2 = F_25 hero")
    ap.add_argument("--out", default="frob_hero_F25.mp4")
    ap.add_argument("--w", type=int, default=1000)
    ap.add_argument("--h", type=int, default=800)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--frames", type=int, default=600, help="one beat (loop length = frames/fps s)")
    ap.add_argument("--crf", type=int, default=18)
    a = ap.parse_args()
    render(a.level, a.out, a.w, a.h, a.fps, a.frames, a.crf)

if __name__ == "__main__":
    main()
