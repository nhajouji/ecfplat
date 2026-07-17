# Shaders

GLSL prototypes of the Frobenius-flow animation (see the Background page applet
`basics_viz.frobenius_flow_html` and `scripts/frobenius_flow_gif.py`).

- **`frobenius_galaxy.glsl`** — paste into a new Shadertoy shader's *Image* tab and
  hit Alt+Enter. Shadertoy provides `iResolution`/`iTime` and calls `mainImage`.
  Curve `y^2 = x^3 + 3x` over `F_5`, `alpha = -2+i`, points `= ker(alpha^3-1)`
  (130), flow `z -> alpha^phi z` in the multiplicative "galaxy" model. Set
  `#define VIEW_GALAXY 0` for the flat torus. Tunables: `N`, the `iTime*0.08`
  speed, the `0.016`/`0.011` dot & trail radii.
- **`frobenius_galaxy_preview.html`** — the same shader in a standalone WebGL2
  page (open it in a browser, no Shadertoy account needed). Useful for local
  tweaking; the `<script id=frag>` block is the shader.

- **`frobenius_galaxy_F5_6.html`** — the big one: `E(F_5^6) = ker(alpha^6-1)` = **15,860
  points**, drawn as GL point-sprites (one vertex per point) so it stays real-time at
  this scale — this is the technique that does NOT fit Shadertoy (no vertex buffers).
  Point size encodes field degree (simpler subfield = bigger): deg1 F5 (10) > deg2 F25
  (10) > deg3 F125 (120) > deg6 (15,720). Shows F25 and F125 as incomparable subfields
  inside the compositum F_5^6. Points + degrees are enumerated on the CPU (box+dedup;
  beta6=-118-44i is not primitive, so the cyclic k/beta trick does NOT apply here).
  Open in a browser; no data files. Tunables: sizes in the vertex shader, SPEED.
