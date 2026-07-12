"""Interactive canvas widgets ported from the MAA MathFest talk slides.

Each builder returns a self-contained HTML document (canvas + JS, KaTeX from a
CDN) meant to be dropped into a Streamlit page via ``st.components.v1.html``.
The widgets are illustrations only -- they carry their own small, validated data
and do not touch the app's precomputed stores, so they render for any visitor
regardless of the ``p <= 1024`` interactive data cap.

Ported graphics:
  * ``cm_torus_html`` -- the interactive CM torus (slide 7): p / trace / class
    sliders drive the N = p+1-a fixed points of x |-> alpha on a fundamental
    domain of C/Lambda, with KaTeX stats.
  * ``volcano_html`` -- the disc -368 2-isogeny volcano (slide 11): pills switch
    the labels between lattice classes and j-invariants mod p across many (a,p)
    sharing the same discriminant -- one graph, many labelings.
  * ``gallery_html`` -- the gallery and all its labelings (slide 14): the disc
    -368 graph anchored at (a,p) = (6,101) with 2- and 3-isogeny edges;
    rotate/reflect/reset act by the class group x {+-1}, restamping every label.
"""

# Shared dark panel + KaTeX head. The panel gives the widget its own theming so
# it looks like an intentional figure on either a light or dark Streamlit theme.
_HEAD = r"""
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css">
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.js"></script>
<style>
  :root { --accent: #4da3d8; --muted: #9aa4ad; --panel: #17191c; --ink: #d7d9dc; }
  * { box-sizing: border-box; }
  body { margin: 0; font-family: -apple-system, "Avenir Next", Helvetica, Arial, sans-serif; }
  .panel {
    background: var(--panel); color: var(--ink);
    border: 1px solid #2b2f34; border-radius: 12px;
    padding: 16px 18px; max-width: 760px; margin: 0 auto;
  }
  .controls { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 10px 22px; margin-bottom: 12px; }
  .ctl label { display: block; font-size: .8rem; color: var(--muted); margin-bottom: 2px; }
  .ctl .val { color: var(--accent); font-weight: 600; }
  .ctl input[type=range] { width: 100%; accent-color: var(--accent); }
  .stage { display: grid; grid-template-columns: auto minmax(150px, 1fr); gap: 16px; align-items: center; }
  canvas { max-width: 100%; height: auto; background: transparent; }
  .stats { font-size: .98rem; line-height: 1.85; }
  .stats > div { color: var(--ink); }
  @media (max-width: 560px) { .controls, .stage { grid-template-columns: 1fr; } }
</style>
"""


def cm_torus_html() -> str:
    """The interactive CM-torus widget (slide 7 of the talk)."""
    return _HEAD + r"""
<div class="panel">
  <div class="controls">
    <div class="ctl">
      <label>prime <span class="val" id="pOut">101</span></label>
      <input type="range" id="pSlider" min="0" max="24" value="17">
    </div>
    <div class="ctl">
      <label>trace <span class="val">a = <span id="aOut">6</span></span></label>
      <input type="range" id="aSlider" min="1" max="19" value="6">
    </div>
    <div class="ctl">
      <label>lattice class <span class="val" id="clsOut">1 / 1</span></label>
      <input type="range" id="clsSlider" min="0" max="0" value="0">
    </div>
  </div>
  <div class="stage">
    <canvas id="torus" width="340" height="310"></canvas>
    <div class="stats" id="stats"></div>
  </div>
</div>
<script>
"use strict";
const ACCENT = "#4da3d8";
const PRIMES = [11,13,17,19,23,29,31,37,41,43,47,53,59,61,71,79,89,101,113,127,149,173,199,229,257];
const pS = document.getElementById("pSlider"), aS = document.getElementById("aSlider"),
      cS = document.getElementById("clsSlider");
pS.max = PRIMES.length - 1; pS.value = PRIMES.indexOf(101);

const gcd = (x, y) => { x = Math.abs(x); y = Math.abs(y); while (y) [x, y] = [y, x % y]; return x; };
const mod = (x, n) => ((x % n) + n) % n;

// primitive reduced forms (A,B,C) of discriminant D < 0
function reducedForms(D) {
  const forms = [];
  const B0 = mod(D, 2);
  for (let B = B0; B * B <= -D / 3 + 1e-9; B += 2) {
    const AC = (B * B - D) / 4;
    for (let A = Math.max(B, 1); A * A <= AC; A++) {
      if (AC % A) continue;
      const C = AC / A;
      if (gcd(gcd(A, B), C) !== 1) continue;
      forms.push([A, B, C]);
      if (B > 0 && B < A && A < C) forms.push([A, -B, C]);
    }
  }
  return forms.sort((f, g) => f[0] - g[0] || Math.abs(f[1]) - Math.abs(g[1]) || g[1] - f[1]);
}

// default to a visually interesting (non-rectangular) class
function pickDefaultClass(forms) {
  let best = -1, bestA = 0;
  forms.forEach(([A, B], i) => { if (B !== 0 && Math.abs(B) < A && A > bestA) { bestA = A; best = i; } });
  if (best >= 0) return best;
  forms.forEach(([A, B], i) => { if (best < 0 && B !== 0) best = i; });
  return best >= 0 ? best : 0;
}

let curDisc = null;
function drawTorus() {
  const p = PRIMES[+pS.value];
  const aMax = Math.floor(2 * Math.sqrt(p) - 1e-9);
  aS.max = aMax; if (+aS.value > aMax) aS.value = aMax;
  const a = +aS.value;
  document.getElementById("pOut").textContent = p;
  document.getElementById("aOut").textContent = a;

  const D = a * a - 4 * p, N = p + 1 - a;
  const forms = reducedForms(D);
  if (D !== curDisc) { curDisc = D; cS.max = forms.length - 1; cS.value = pickDefaultClass(forms); }
  const cls = Math.min(+cS.value, forms.length - 1);
  const [A, B, C] = forms[cls];
  document.getElementById("clsOut").textContent = `(${A}, ${B}, ${C}) · ${cls + 1} / ${forms.length}`;

  // Lambda = <1, tau>, tau = (-B + sqrt(D)) / 2A
  const tx = -B / (2 * A), ty = Math.sqrt(-D) / (2 * A);
  // mult by (alpha - 1) on basis (1, tau) is the integer matrix [[s,u],[t,v]]
  const s = (a - 2 + B) / 2, t = A, u = -C, v = (a - 2 - B) / 2;
  // fixed points = columns of the adjugate generate a subgroup of (Z/N)^2
  const g1x = mod(v, N), g1y = mod(-t, N), g2x = mod(-u, N), g2y = mod(s, N);
  const seen = new Set(), pts = [];
  outer:
  for (let i = 0; i < N; i++) {
    for (let j = 0; j < N; j++) {
      const x = mod(i * g1x + j * g2x, N), y = mod(i * g1y + j * g2y, N);
      const key = x * N + y;
      if (!seen.has(key)) { seen.add(key); pts.push([x, y]); if (pts.length === N) break outer; }
    }
  }

  const cv = document.getElementById("torus"), ctx = cv.getContext("2d");
  ctx.clearRect(0, 0, cv.width, cv.height);
  const pad = 30;
  const scale = Math.min((cv.width - 2 * pad) / (1 + Math.abs(tx)), (cv.height - 2 * pad) / ty);
  const ox = pad + (tx < 0 ? -tx * scale : 0), oy = cv.height - pad;
  const X = (uu, vv) => ox + (uu + vv * tx) * scale;
  const Y = (uu, vv) => oy - vv * ty * scale;

  ctx.strokeStyle = "#555"; ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.moveTo(X(0,0), Y(0,0)); ctx.lineTo(X(1,0), Y(1,0));
  ctx.lineTo(X(1,1), Y(1,1)); ctx.lineTo(X(0,1), Y(0,1)); ctx.closePath(); ctx.stroke();

  pts.forEach(([x, y], k) => {
    ctx.fillStyle = `hsl(${205 + 130 * k / N}, 65%, ${58 - 18 * k / N}%)`;
    ctx.beginPath(); ctx.arc(X(x / N, y / N), Y(x / N, y / N), Math.max(3.5, 9 - N / 40), 0, 7); ctx.fill();
  });

  const stats = document.getElementById("stats");
  stats.innerHTML = "";
  const line = tex => {
    const d = document.createElement("div");
    stats.appendChild(d);
    if (window.katex) katex.render(tex, d, { throwOnError: false });
    else d.textContent = tex;
  };
  line(`\\chi(x) = x^2 - ${a}x + ${p}`);
  line(`\\#\\,\\mathrm{Fix}([\\alpha]) = \\chi(1) = ${N}`);
  line(`\\mathrm{disc} = a^2 - 4p = ${D},\\;\\; h = ${forms.length}`);
  line(`\\tau = ${tx === 0 ? "" : (tx > 0 ? tx.toFixed(3) + " + " : "-" + Math.abs(tx).toFixed(3) + " + ")}${ty.toFixed(3)}\\,i`);
}
pS.addEventListener("input", drawTorus);
aS.addEventListener("input", drawTorus);
cS.addEventListener("input", drawTorus);

function boot() { drawTorus(); }
if (window.katex) boot();
else {
  const t = setInterval(() => { if (window.katex) { clearInterval(t); boot(); } }, 40);
  boot();  // draw immediately; stats re-render once KaTeX lands
}
</script>
"""


# Dark panel + pills/buttons styling for the two isogeny-graph widgets.
_GRAPH_HEAD = r"""
<style>
  :root { --accent: #4da3d8; --muted: #9aa4ad; --panel: #17191c; --ink: #d7d9dc; }
  * { box-sizing: border-box; }
  body { margin: 0; font-family: -apple-system, "Avenir Next", Helvetica, Arial, sans-serif; }
  .panel {
    background: var(--panel); color: var(--ink);
    border: 1px solid #2b2f34; border-radius: 12px;
    padding: 14px 16px; max-width: 780px; margin: 0 auto;
  }
  .bar { display: flex; flex-wrap: wrap; gap: 8px; margin-bottom: 10px; }
  .pill, .btn {
    font: 500 .9rem/1 -apple-system, "Avenir Next", Helvetica, sans-serif;
    color: var(--ink); background: #22262b; border: 1px solid #33383e;
    border-radius: 999px; padding: 7px 13px; cursor: pointer; transition: all .12s;
  }
  .pill:hover, .btn:hover { border-color: var(--accent); }
  .pill.on { background: var(--accent); border-color: var(--accent); color: #0c1013; }
  .btn { border-radius: 8px; }
  canvas#volcano, canvas#galvol { display: block; width: 100%; height: auto; background: transparent; }
  .cards { display: grid; grid-template-columns: repeat(6, 1fr); gap: 6px; margin-top: 10px; }
  .card { border: 1px solid #2e2e2e; border-radius: 7px; padding: 5px 4px;
          display: flex; flex-direction: column; align-items: center; gap: 3px; }
  .card canvas { width: 100%; height: auto; }
  .cap { font-size: .78rem; color: var(--muted); text-align: center; line-height: 1.25; }
  @media (max-width: 520px) { .cards { grid-template-columns: repeat(3, 1fr); } }
</style>
"""

# disc -368 2-isogeny volcano; form<->j bijections from ecfplat
# ecqf_ord_pcbij_ext.json, curve-side adjacency verified via Phi_2 mod p.
_VOLCANO_JSON = r"""{"disc": -368, "verts": [{"form": [1, 1, 6], "level": 0}, {"form": [2, -1, 3], "level": 0}, {"form": [2, 1, 3], "level": 0}, {"form": [1, 0, 23], "level": 1}, {"form": [3, -2, 8], "level": 1}, {"form": [3, 2, 8], "level": 1}, {"form": [1, 0, 92], "level": 2}, {"form": [3, -2, 31], "level": 2}, {"form": [3, 2, 31], "level": 2}, {"form": [4, 0, 23], "level": 2}, {"form": [9, -8, 12], "level": 2}, {"form": [9, 8, 12], "level": 2}], "edges": [[0, 1], [0, 2], [0, 3], [1, 2], [1, 5], [2, 4], [3, 6], [3, 9], [4, 8], [4, 10], [5, 7], [5, 11]], "cats": [{"label": "lattices", "js": null}, {"label": "(6, 101)", "js": ["30", "28", "65", "98", "60", "27", "22", "45", "39", "74", "62", "68"]}, {"label": "(18, 173)", "js": ["23", "14", "45", "53", "148", "128", "1", "112", "37", "113", "101", "77"]}, {"label": "(30, 317)", "js": ["300", "107", "232", "216", "171", "70", "46", "205", "170", "71", "108", "296"]}, {"label": "(54, 821)", "js": ["629", "819", "157", "15", "24", "815", "101", "417", "139", "426", "134", "602"]}, {"label": "(66, 1181)", "js": ["895", "1112", "822", "165", "1090", "810", "21", "922", "337", "658", "1013", "516"]}, {"label": "(78, 1613)", "js": ["371", "1327", "310", "810", "1586", "901", "162", "1306", "627", "1020", "343", "1470"]}]}"""


def volcano_html() -> str:
    """The disc -368 volcano with label-switching pills (slide 11)."""
    return _GRAPH_HEAD + r"""
<div class="panel">
  <div id="volcats" class="bar"></div>
  <canvas id="volcano" width="700" height="600"></canvas>
</div>
<script>
"use strict";
const VOLCANO = """ + _VOLCANO_JSON + r""";
let volPos = null;
function volLayout(cv, sc = 1) {
  const cx = cv.width / 2, cy = cv.height / 2 + 24;
  const r1 = 95 * sc, r2 = 200 * sc, r3 = 300 * sc;
  const pos = new Array(VOLCANO.verts.length);
  const crater = VOLCANO.verts.map((v, i) => v.level === 0 ? i : -1).filter(i => i >= 0);
  crater.forEach((c, k) => {
    const th = -Math.PI / 2 + k * 2 * Math.PI / crater.length;
    pos[c] = { th, r: r1 };
  });
  VOLCANO.edges.forEach(([i, j]) => {
    const [c, m] = [i, j].sort((x, y) => VOLCANO.verts[x].level - VOLCANO.verts[y].level);
    if (VOLCANO.verts[c].level === 0 && VOLCANO.verts[m].level === 1) pos[m] = { th: pos[c].th, r: r2 };
  });
  const kids = {};
  VOLCANO.edges.forEach(([i, j]) => {
    const [m, f] = [i, j].sort((x, y) => VOLCANO.verts[x].level - VOLCANO.verts[y].level);
    if (VOLCANO.verts[m].level === 1 && VOLCANO.verts[f].level === 2) {
      kids[m] = (kids[m] || 0) + 1;
      pos[f] = { th: pos[m].th + (kids[m] === 1 ? -0.30 : 0.30), r: r3 };
    }
  });
  return pos.map(({ th, r }) => ({ x: cx + r * Math.cos(th), y: cy + r * Math.sin(th), th, r }));
}
function drawVolcano(catIdx) {
  const cv = document.getElementById("volcano"), ctx = cv.getContext("2d");
  if (!volPos) volPos = volLayout(cv);
  const cat = VOLCANO.cats[catIdx];
  ctx.clearRect(0, 0, cv.width, cv.height);
  ctx.strokeStyle = "#4a4a4a"; ctx.lineWidth = 2.5;
  VOLCANO.edges.forEach(([i, j]) => {
    ctx.beginPath(); ctx.moveTo(volPos[i].x, volPos[i].y); ctx.lineTo(volPos[j].x, volPos[j].y); ctx.stroke();
  });
  const lvlColor = ["#e0b64f", "#4da3d8", "#8fb4c7"];
  VOLCANO.verts.forEach((v, i) => {
    ctx.fillStyle = lvlColor[v.level];
    ctx.beginPath(); ctx.arc(volPos[i].x, volPos[i].y, 12, 0, 7); ctx.fill();
  });
  ctx.font = "500 22px 'Avenir Next', Helvetica, sans-serif"; ctx.fillStyle = "#cfcdc8";
  ctx.textAlign = "center"; ctx.textBaseline = "middle";
  VOLCANO.verts.forEach((v, i) => {
    const label = cat.js ? cat.js[i] : "(" + v.form.join(",") + ")";
    const { x, y, th } = volPos[i];
    let lx, ly;
    if (v.level === 0) { lx = x - 44 * Math.cos(th); ly = y - 44 * Math.sin(th); }
    else if (v.level === 1) { lx = x - 36 * Math.sin(th); ly = y + 36 * Math.cos(th); }
    else { lx = x + 46 * Math.cos(th); ly = y + 46 * Math.sin(th); }
    ctx.fillText(label, lx, ly);
  });
  ctx.textAlign = "left";
  ctx.fillStyle = "#777"; ctx.font = "18px 'Avenir Next', Helvetica, sans-serif";
  ctx.fillText(cat.js ? "j-invariants mod " + cat.label.split(", ")[1].slice(0, -1)
                      : "lattice classes ⟨1, τ⟩ — reduced forms", 20, 26);
  document.querySelectorAll("#volcats .pill").forEach((b, k) => b.classList.toggle("on", k === catIdx));
}
function buildVolcanoUI() {
  const bar = document.getElementById("volcats");
  VOLCANO.cats.forEach((cat, k) => {
    const b = document.createElement("button");
    b.className = "pill";
    b.textContent = cat.label === "lattices" ? "𝒜(α) · lattices" : "(a,p) = " + cat.label;
    b.addEventListener("click", () => drawVolcano(k));
    bar.appendChild(b);
  });
  drawVolcano(0);
}
buildVolcanoUI();
</script>
"""


# disc -368 anchored at (a,p)=(6,101); rot/refl generate all 12 valid labelings.
_GALLERY_JSON = r"""{"a": 6, "p": 101, "verts": [{"form": [1, 1, 6], "level": 0}, {"form": [2, -1, 3], "level": 0}, {"form": [2, 1, 3], "level": 0}, {"form": [1, 0, 23], "level": 1}, {"form": [3, -2, 8], "level": 1}, {"form": [3, 2, 8], "level": 1}, {"form": [1, 0, 92], "level": 2}, {"form": [3, -2, 31], "level": 2}, {"form": [3, 2, 31], "level": 2}, {"form": [4, 0, 23], "level": 2}, {"form": [9, -8, 12], "level": 2}, {"form": [9, 8, 12], "level": 2}], "edges2": [[0, 1], [0, 2], [0, 3], [1, 2], [1, 5], [2, 4], [3, 6], [3, 9], [4, 8], [4, 10], [5, 7], [5, 11]], "edges3": [[0, 1], [0, 2], [1, 2], [3, 4], [3, 5], [4, 5], [6, 7], [6, 8], [7, 10], [8, 11], [9, 10], [9, 11]], "rot": [1, 2, 0, 5, 3, 4, 7, 10, 6, 11, 9, 8], "refl": [0, 2, 1, 3, 5, 4, 6, 8, 7, 9, 11, 10], "js": [30, 28, 65, 98, 60, 27, 22, 45, 39, 74, 62, 68], "eqs": [[28, 36], [45, 10], [98, 87], [1, 15], [68, 68], [68, 3], [82, 72], [56, 10], [26, 28], [10, 32], [8, 31], [88, 90]]}"""


def gallery_html() -> str:
    """The gallery and its 12 labelings, with rotate/reflect/reset (slide 14)."""
    return _GRAPH_HEAD + r"""
<div class="panel">
  <div class="bar">
    <button class="btn" id="galL">⟲ rotate</button>
    <button class="btn" id="galR">rotate ⟳</button>
    <button class="btn" id="galC">reflect</button>
    <button class="btn" id="gal0">reset</button>
  </div>
  <canvas id="galvol" width="620" height="400"></canvas>
  <div class="cards" id="galcards"></div>
</div>
<script>
"use strict";
const GALLERY = """ + _GALLERY_JSON + r""";
const mod = (x, n) => ((x % n) + n) % n;
let galLab = null, galPos = null, galRotInv = null;

function galLayout(cv, sc) {
  const cx = cv.width / 2, cy = cv.height / 2 + 12;
  const r1 = 78 * sc, r2 = 150 * sc, r3 = 210 * sc;
  const pos = new Array(GALLERY.verts.length);
  const crater = GALLERY.verts.map((v, i) => v.level === 0 ? i : -1).filter(i => i >= 0);
  crater.forEach((c, k) => { pos[c] = { th: -Math.PI / 2 + k * 2 * Math.PI / crater.length, r: r1 }; });
  GALLERY.edges2.forEach(([i, j]) => {
    const [c, m] = [i, j].sort((x, y) => GALLERY.verts[x].level - GALLERY.verts[y].level);
    if (GALLERY.verts[c].level === 0 && GALLERY.verts[m].level === 1) pos[m] = { th: pos[c].th, r: r2 };
  });
  const kids = {};
  GALLERY.edges2.forEach(([i, j]) => {
    const [m, f] = [i, j].sort((x, y) => GALLERY.verts[x].level - GALLERY.verts[y].level);
    if (GALLERY.verts[m].level === 1 && GALLERY.verts[f].level === 2) {
      kids[m] = (kids[m] || 0) + 1;
      pos[f] = { th: pos[m].th + (kids[m] === 1 ? -0.30 : 0.30), r: r3 };
    }
  });
  return pos.map(({ th, r }) => ({ x: cx + r * Math.cos(th), y: cy + r * Math.sin(th), th, r }));
}

const galFloorCycle = (() => {
  const floor = GALLERY.verts.map((v, i) => v.level === 2 ? i : -1).filter(i => i >= 0);
  const nbr = {}; floor.forEach(v => nbr[v] = []);
  GALLERY.edges3.forEach(([i, j]) => { if (i in nbr && j in nbr) { nbr[i].push(j); nbr[j].push(i); } });
  const cyc = [floor[0]]; let prev = null;
  while (cyc.length < floor.length) {
    const nxt = nbr[cyc[cyc.length - 1]].find(w => w !== prev);
    prev = cyc[cyc.length - 1]; cyc.push(nxt);
  }
  return cyc;
})();

function drawMiniTorus(cv, form) {
  const [A, B, C] = form, a = GALLERY.a, p = GALLERY.p, D = B * B - 4 * A * C, N = p + 1 - a;
  const tx = -B / (2 * A), ty = Math.sqrt(-D) / (2 * A);
  const s = (a - 2 + B) / 2, t = A, u = -C, v = (a - 2 - B) / 2;
  const g1x = mod(v, N), g1y = mod(-t, N), g2x = mod(-u, N), g2y = mod(s, N);
  const seen = new Set(), pts = [];
  outer:
  for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
    const x = mod(i * g1x + j * g2x, N), y = mod(i * g1y + j * g2y, N);
    const key = x * N + y;
    if (!seen.has(key)) { seen.add(key); pts.push([x, y]); if (pts.length === N) break outer; }
  }
  const ctx = cv.getContext("2d"), pad = 6;
  ctx.clearRect(0, 0, cv.width, cv.height);
  const scale = Math.min((cv.width - 2 * pad) / (1 + Math.abs(tx)), (cv.height - 2 * pad) / ty);
  const ox = pad + (tx < 0 ? -tx * scale : 0), oy = cv.height - (cv.height - ty * scale) / 2;
  const X = (uu, vv) => ox + (uu + vv * tx) * scale, Y = (uu, vv) => oy - vv * ty * scale;
  ctx.strokeStyle = "#444"; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(X(0,0), Y(0,0)); ctx.lineTo(X(1,0), Y(1,0));
  ctx.lineTo(X(1,1), Y(1,1)); ctx.lineTo(X(0,1), Y(0,1)); ctx.closePath(); ctx.stroke();
  pts.forEach(([x, y], k) => {
    ctx.fillStyle = `hsl(${205 + 130 * k / N}, 60%, ${55 - 15 * k / N}%)`;
    ctx.beginPath(); ctx.arc(X(x / N, y / N), Y(x / N, y / N), 1.8, 0, 7); ctx.fill();
  });
}

function drawGalVol() {
  const cv = document.getElementById("galvol"), ctx = cv.getContext("2d");
  if (!galPos) galPos = galLayout(cv, 1.0);
  ctx.clearRect(0, 0, cv.width, cv.height);
  const cx = cv.width / 2, cy = cv.height / 2 + 12;
  ctx.strokeStyle = "#6b5a2e"; ctx.lineWidth = 2; ctx.setLineDash([6, 6]);
  GALLERY.edges3.forEach(([i, j]) => {
    const a = galPos[i], b = galPos[j];
    const mx = (a.x + b.x) / 2, my = (a.y + b.y) / 2;
    const dx = mx - cx, dy = my - cy, n = Math.hypot(dx, dy) || 1;
    ctx.beginPath(); ctx.moveTo(a.x, a.y);
    ctx.quadraticCurveTo(mx + 30 * dx / n, my + 30 * dy / n, b.x, b.y); ctx.stroke();
  });
  ctx.setLineDash([]);
  ctx.strokeStyle = "#4a4a4a"; ctx.lineWidth = 2.2;
  GALLERY.edges2.forEach(([i, j]) => {
    ctx.beginPath(); ctx.moveTo(galPos[i].x, galPos[i].y); ctx.lineTo(galPos[j].x, galPos[j].y); ctx.stroke();
  });
  const lvlColor = ["#e0b64f", "#4da3d8", "#8fb4c7"];
  GALLERY.verts.forEach((v, i) => {
    ctx.fillStyle = lvlColor[v.level];
    ctx.beginPath(); ctx.arc(galPos[i].x, galPos[i].y, 8, 0, 7); ctx.fill();
  });
  ctx.font = "500 15px 'Avenir Next', Helvetica, sans-serif"; ctx.fillStyle = "#cfcdc8";
  ctx.textAlign = "center"; ctx.textBaseline = "middle";
  GALLERY.verts.forEach((v, i) => {
    const { x, y, th } = galPos[i];
    let lx, ly;
    if (v.level === 0) { lx = x - 28 * Math.cos(th); ly = y - 28 * Math.sin(th); }
    else if (v.level === 1) { lx = x - 24 * Math.sin(th); ly = y + 24 * Math.cos(th); }
    else { lx = x + 30 * Math.cos(th); ly = y + 30 * Math.sin(th); }
    ctx.fillText("j=" + GALLERY.js[galLab[i]], lx, ly);
  });
  ctx.textAlign = "left"; ctx.fillStyle = "#777"; ctx.font = "16px 'Avenir Next', Helvetica, sans-serif";
  ctx.fillText("one labeling of " + GALLERY.verts.length + " classes · (a, p) = (" + GALLERY.a + ", " + GALLERY.p + ")", 14, 22);
}

function galRefreshCards() {
  galFloorCycle.forEach((v, k) => {
    const cap = document.getElementById("galcap" + k);
    const [f, g] = GALLERY.eqs[galLab[v]];
    cap.textContent = "j = " + GALLERY.js[galLab[v]];
    cap.parentElement.title = "y² = x³ + " + f + "x + " + g;  // full equation on hover
  });
}
function galApply(fn) { galLab = fn(galLab); drawGalVol(); galRefreshCards(); }
function buildGalleryUI() {
  galLab = GALLERY.verts.map((_, i) => i);
  galRotInv = new Array(12); GALLERY.rot.forEach((w, v) => galRotInv[w] = v);
  const grid = document.getElementById("galcards");
  galFloorCycle.forEach((v, k) => {
    const card = document.createElement("div"); card.className = "card";
    const mini = document.createElement("canvas");
    mini.width = 150; mini.height = 92;
    const cap = document.createElement("div"); cap.id = "galcap" + k; cap.className = "cap";
    card.appendChild(mini); card.appendChild(cap); grid.appendChild(card);
    drawMiniTorus(mini, GALLERY.verts[v].form);
  });
  document.getElementById("galL").addEventListener("click", () => galApply(l => l.map((_, w) => l[GALLERY.rot[w]])));
  document.getElementById("galR").addEventListener("click", () => galApply(l => l.map((_, w) => l[galRotInv[w]])));
  document.getElementById("galC").addEventListener("click", () => galApply(l => l.map((_, w) => l[GALLERY.refl[w]])));
  document.getElementById("gal0").addEventListener("click", () => galApply(() => GALLERY.verts.map((_, i) => i)));
  drawGalVol(); galRefreshCards();
}
buildGalleryUI();
</script>
"""
