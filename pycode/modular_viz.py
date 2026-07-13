"""Interactive widgets for the Modular Curves pages.

Original (not talk-ported) canvas/JS illustrations, embedded via
``st.components.v1.html``. Self-contained, no external data.

Conventions (see the page prose): SL(2,Z) acts on the LEFT of the upper
half-plane, γτ = (aτ+b)/(cτ+d), so the moduli space is written Γ\\H (group on
the left). These widgets are KaTeX-free -- all labels are plain Unicode.

Widgets:
  * ``moduli_applet_html`` -- the moduli space of lattices. One applet, two
    modes toggled: drive the shape τ in Γ\\H and watch the lattice Z+Zτ, or drag
    a basis (ω1, ω2) in C and watch τ land (and reduce) in the fundamental
    domain.
"""

_HEAD = r"""
<style>
  :root { --accent:#4da3d8; --gold:#e0b64f; --muted:#9aa4ad; --panel:#17191c; --ink:#d7d9dc; }
  * { box-sizing: border-box; }
  body { margin:0; font-family:-apple-system,"Avenir Next",Helvetica,Arial,sans-serif; }
  .panel { background:var(--panel); color:var(--ink); border:1px solid #2b2f34;
           border-radius:12px; padding:14px 16px; max-width:760px; margin:0 auto; }
  .modebar { display:flex; gap:8px; margin-bottom:12px; }
  .seg { font:600 .9rem/1 -apple-system,"Avenir Next",Helvetica,sans-serif; color:var(--ink);
         background:#22262b; border:1px solid #33383e; border-radius:8px; padding:8px 14px; cursor:pointer; transition:all .12s; }
  .seg:hover { border-color:var(--accent); }
  .seg.on { background:var(--accent); border-color:var(--accent); color:#0c1013; }
  .stage { display:grid; grid-template-columns:1fr 1fr; gap:14px; }
  .cell { display:flex; flex-direction:column; gap:5px; }
  .cell .cap { font-size:.82rem; color:var(--muted); text-align:center; }
  canvas { width:100%; height:auto; background:transparent; border:1px solid #23272c; border-radius:8px; touch-action:none; }
  .info { margin-top:11px; font-size:.92rem; color:var(--ink); text-align:center; line-height:1.5; }
  .hint { margin-top:4px; font-size:.8rem; color:var(--muted); text-align:center; }
  @media (max-width:560px){ .stage { grid-template-columns:1fr; } }
</style>
"""


def moduli_applet_html() -> str:
    """The moduli-of-lattices applet (Γ\\H <-> a lattice in C), one toggle."""
    return _HEAD + r"""
<div class="panel">
  <div class="modebar">
    <button class="seg on" id="segTau">drive the shape τ</button>
    <button class="seg" id="segBasis">drag a basis</button>
  </div>
  <div class="stage">
    <div class="cell">
      <div class="cap">the shape τ ∈ Γ∖ℍ</div>
      <canvas id="mvH" width="340" height="340"></canvas>
    </div>
    <div class="cell">
      <div class="cap">a lattice in ℂ</div>
      <canvas id="mvL" width="340" height="340"></canvas>
    </div>
  </div>
  <div class="info" id="mvInfo"></div>
  <div class="hint" id="mvHint"></div>
</div>
<script>
"use strict";
const ACCENT="#4da3d8", GOLD="#e0b64f", INK="#d7d9dc", MUT="#9aa4ad", GRID="#3a3f45";

// ---------- complex helpers ----------
const C=(re,im)=>({re,im});
const cabs2=z=>z.re*z.re+z.im*z.im;
const cdiv=(a,b)=>{const d=cabs2(b); return C((a.re*b.re+a.im*b.im)/d,(a.im*b.re-a.re*b.im)/d);};

// reduce τ ∈ ℍ to the standard fundamental domain (|Re|<=1/2, |τ|>=1)
function reduceFD(t){
  let z=C(t.re,t.im), guard=0;
  while(guard++<200){
    const n=Math.round(z.re); z=C(z.re-n, z.im);
    if(cabs2(z) < 1-1e-9){ const d=cabs2(z); z=C(-z.re/d, z.im/d); }  // S: τ ↦ -1/τ
    else break;
  }
  return z;
}

// ---------- state ----------
let mode="tau";
let tau=C(0.28,0.95);
let b1=C(1,0), b2=C(0.28,0.95);
let dragging=null, handleA=[0,0], handleB=[0,0];

const hCv=document.getElementById("mvH"), hCtx=hCv.getContext("2d");
const lCv=document.getElementById("mvL"), lCtx=lCv.getContext("2d");
const info=document.getElementById("mvInfo"), hint=document.getElementById("mvHint");

// lattice canvas geometry
const S=66, cx=lCv.width/2, cy=lCv.height*0.60;
const toPx=z=>[cx+S*z.re, cy-S*z.im];

// ℍ canvas mapping
const HW={xmin:-1.7,xmax:1.7,ymin:0.0,ymax:2.7};
const hX=re=>(re-HW.xmin)/(HW.xmax-HW.xmin)*hCv.width;
const hY=im=>hCv.height-(im-HW.ymin)/(HW.ymax-HW.ymin)*hCv.height;
const hInv=(px,py)=>C(HW.xmin+px/hCv.width*(HW.xmax-HW.xmin), HW.ymin+(hCv.height-py)/hCv.height*(HW.ymax-HW.ymin));

function dot(ctx,x,y,r){ctx.beginPath();ctx.arc(x,y,r,0,7);ctx.fill();}
function arrow(ctx,p,q,color){
  ctx.strokeStyle=color; ctx.fillStyle=color; ctx.lineWidth=2;
  ctx.beginPath(); ctx.moveTo(p[0],p[1]); ctx.lineTo(q[0],q[1]); ctx.stroke();
  const ang=Math.atan2(q[1]-p[1],q[0]-p[0]), h=8;
  ctx.beginPath(); ctx.moveTo(q[0],q[1]);
  ctx.lineTo(q[0]-h*Math.cos(ang-0.4), q[1]-h*Math.sin(ang-0.4));
  ctx.lineTo(q[0]-h*Math.cos(ang+0.4), q[1]-h*Math.sin(ang+0.4));
  ctx.closePath(); ctx.fill();
}

function drawH(active, rep){
  hCtx.clearRect(0,0,hCv.width,hCv.height);
  // fundamental domain
  hCtx.beginPath();
  hCtx.moveTo(hX(-0.5),0);
  hCtx.lineTo(hX(-0.5),hY(Math.sqrt(3)/2));
  for(let a=120;a>=60;a-=2){const r=a*Math.PI/180; hCtx.lineTo(hX(Math.cos(r)),hY(Math.sin(r)));}
  hCtx.lineTo(hX(0.5),0);
  hCtx.closePath();
  hCtx.fillStyle="rgba(77,163,216,0.10)"; hCtx.fill();
  hCtx.strokeStyle="rgba(77,163,216,0.45)"; hCtx.lineWidth=1.5; hCtx.stroke();
  // real axis + unit circle
  hCtx.strokeStyle=GRID; hCtx.lineWidth=1;
  hCtx.beginPath(); hCtx.moveTo(0,hY(0)); hCtx.lineTo(hCv.width,hY(0)); hCtx.stroke();
  hCtx.strokeStyle="rgba(255,255,255,0.10)";
  hCtx.beginPath();
  for(let a=0;a<=180;a+=2){const r=a*Math.PI/180; const px=hX(Math.cos(r)),py=hY(Math.sin(r)); a===0?hCtx.moveTo(px,py):hCtx.lineTo(px,py);}
  hCtx.stroke();
  // reduced representative (gold ring) + connector
  if(rep){
    const same=Math.abs(rep.re-active.re)<1e-3 && Math.abs(rep.im-active.im)<1e-3;
    if(!same){
      hCtx.strokeStyle="rgba(224,182,79,0.5)"; hCtx.setLineDash([4,4]);
      hCtx.beginPath(); hCtx.moveTo(hX(active.re),hY(active.im)); hCtx.lineTo(hX(rep.re),hY(rep.im)); hCtx.stroke();
      hCtx.setLineDash([]);
      hCtx.strokeStyle=GOLD; hCtx.lineWidth=2;
      hCtx.beginPath(); hCtx.arc(hX(rep.re),hY(rep.im),6,0,7); hCtx.stroke();
    }
  }
  // the active τ
  hCtx.fillStyle=ACCENT; dot(hCtx,hX(active.re),hY(active.im),6.5);
  hCtx.fillStyle=INK; hCtx.font="13px system-ui";
  hCtx.fillText("τ", hX(active.re)+9, hY(active.im)-8);
}

function drawLattice(bb1,bb2,interactive){
  lCtx.clearRect(0,0,lCv.width,lCv.height);
  // fundamental parallelogram 0, b1, b1+b2, b2
  const p0=toPx(C(0,0)), pa=toPx(bb1), pb=toPx(bb2),
        pab=toPx(C(bb1.re+bb2.re,bb1.im+bb2.im));
  lCtx.fillStyle="rgba(77,163,216,0.10)";
  lCtx.beginPath(); lCtx.moveTo(p0[0],p0[1]); lCtx.lineTo(pa[0],pa[1]);
  lCtx.lineTo(pab[0],pab[1]); lCtx.lineTo(pb[0],pb[1]); lCtx.closePath(); lCtx.fill();
  // lattice points
  lCtx.fillStyle="rgba(215,217,220,0.55)";
  for(let m=-6;m<=6;m++)for(let n=-6;n<=6;n++){
    if(m===0&&n===0)continue;
    const z=C(m*bb1.re+n*bb2.re, m*bb1.im+n*bb2.im);
    const px=cx+S*z.re, py=cy-S*z.im;
    if(px<-2||px>lCv.width+2||py<-2||py>lCv.height+2)continue;
    lCtx.beginPath(); lCtx.arc(px,py,2.2,0,7); lCtx.fill();
  }
  // basis arrows
  arrow(lCtx,p0,pa,ACCENT); arrow(lCtx,p0,pb,GOLD);
  // origin
  lCtx.fillStyle=INK; dot(lCtx,p0[0],p0[1],3.5);
  // labels
  lCtx.fillStyle=INK; lCtx.font="13px system-ui";
  lCtx.fillText(mode==="tau"?"1":"ω₁", pa[0]+7, pa[1]+4);
  lCtx.fillText(mode==="tau"?"τ":"ω₂", pb[0]+7, pb[1]+4);
  // drag handles (basis mode)
  handleA=pa; handleB=pb;
  if(interactive){
    lCtx.fillStyle=ACCENT; dot(lCtx,pa[0],pa[1],7);
    lCtx.fillStyle=GOLD;   dot(lCtx,pb[0],pb[1],7);
  }
}

const fmt=z=>{const a=z.re.toFixed(3), b=Math.abs(z.im).toFixed(3), s=z.im<0?"−":"+"; return `${a} ${s} ${b}i`;};

function render(){
  if(mode==="tau"){
    const rep=reduceFD(tau);
    drawH(tau, rep, true);
    drawLattice(C(1,0), tau, false);
    const same=Math.abs(rep.re-tau.re)<1e-3 && Math.abs(rep.im-tau.im)<1e-3;
    info.innerHTML = `τ = ${fmt(tau)} &nbsp;→&nbsp; lattice ℤ + ℤτ. ` +
      (same ? `This τ is the canonical shape (in the fundamental domain).`
            : `Its canonical shape is <b style="color:${GOLD}">${fmt(rep)}</b> — the same lattice, so the same point of Γ∖ℍ.`);
    hint.textContent = "drag anywhere in the left panel to move τ";
  } else {
    let t=cdiv(b2,b1), ob1=b1, ob2=b2;
    if(t.im<0){ t=cdiv(b1,b2); ob1=b2; ob2=b1; }   // reorient so Im τ > 0
    const rep=reduceFD(t);
    drawH(t, rep, false);
    drawLattice(b1,b2, true);
    const same=Math.abs(rep.re-t.re)<1e-3 && Math.abs(rep.im-t.im)<1e-3;
    info.innerHTML = `τ = ω₂ / ω₁ = ${fmt(t)}` +
      (same ? ` — already the canonical shape.`
            : ` &nbsp;≡&nbsp; <b style="color:${GOLD}">${fmt(rep)}</b> in the fundamental domain.`);
    hint.textContent = "drag the blue (ω₁) and gold (ω₂) arrowheads";
  }
}

// ---------- interaction ----------
function evtToCanvas(cv,e){const r=cv.getBoundingClientRect(); return {x:(e.clientX-r.left)*cv.width/r.width, y:(e.clientY-r.top)*cv.height/r.height};}
const near=(p,h)=>Math.hypot(p.x-h[0],p.y-h[1])<15;
function setTauFromPx(p){const z=hInv(p.x,p.y); tau=C(z.re, Math.max(0.05,z.im)); render();}

hCv.addEventListener("pointerdown",e=>{
  if(mode!=="tau")return;
  dragging="tau"; hCv.setPointerCapture(e.pointerId); setTauFromPx(evtToCanvas(hCv,e)); e.preventDefault();
});
lCv.addEventListener("pointerdown",e=>{
  if(mode!=="basis")return;
  const p=evtToCanvas(lCv,e);
  if(near(p,handleA))dragging="a"; else if(near(p,handleB))dragging="b";
  if(dragging){ lCv.setPointerCapture(e.pointerId); e.preventDefault(); }
});
function moveHandler(e){
  if(!dragging)return;
  if(dragging==="tau"){ setTauFromPx(evtToCanvas(hCv,e)); }
  else{ const p=evtToCanvas(lCv,e); const z=C((p.x-cx)/S, -(p.y-cy)/S); if(dragging==="a")b1=z; else b2=z; render(); }
}
hCv.addEventListener("pointermove",moveHandler);
lCv.addEventListener("pointermove",moveHandler);
window.addEventListener("pointerup",()=>{dragging=null;});

function setMode(m){
  if(m===mode){render();return;}
  if(m==="basis"){ b1=C(1,0); b2=C(tau.re,tau.im); }
  else { let t=cdiv(b2,b1); if(t.im<0)t=cdiv(b1,b2); tau=t; }
  mode=m;
  document.getElementById("segTau").classList.toggle("on",m==="tau");
  document.getElementById("segBasis").classList.toggle("on",m==="basis");
  render();
}
document.getElementById("segTau").addEventListener("click",()=>setMode("tau"));
document.getElementById("segBasis").addEventListener("click",()=>setMode("basis"));

render();
</script>
"""


import json as _json


def hilbert_applet_html(entries) -> str:
    """Singular moduli / Hilbert class polynomials: CM points in the fundamental
    domain, coloured by discriminant. Click a point to light up its
    same-endomorphism-ring siblings and load H_D. ``entries`` is a list of
    ``{d, h, pts:[[re,im]...], poly:<html>, j0}`` (see the page builder)."""
    data = _json.dumps(entries)
    return _HEAD + r"""
<div class="panel">
  <div class="stage" style="grid-template-columns: 300px 1fr;">
    <div class="cell">
      <div class="cap">CM points in Γ∖ℍ</div>
      <canvas id="hcm" width="300" height="470"></canvas>
    </div>
    <div class="cell" style="justify-content:flex-start;">
      <div class="cap" style="text-align:left;">the singular modulus / its class polynomial</div>
      <div id="hpanel" style="border:1px solid #23272c;border-radius:8px;padding:14px;min-height:150px;"></div>
      <div class="hint" style="text-align:left;margin-top:8px;">each colour is one discriminant D; click a point to select its class</div>
    </div>
  </div>
</div>
<script>
"use strict";
const ENTRIES = """ + data + r""";
const ACC="#4da3d8", INK="#d7d9dc", MUT="#9aa4ad";
const cv=document.getElementById("hcm"), ctx=cv.getContext("2d");
const panel=document.getElementById("hpanel");

// view window in ℍ
let Imax=1.5; ENTRIES.forEach(e=>e.pts.forEach(p=>{ if(p[1]>Imax)Imax=p[1]; }));
Imax+=0.35;
const Rmin=-0.62, Rmax=0.62, Imin=0.42;
const X=re=>(re-Rmin)/(Rmax-Rmin)*cv.width;
const Y=im=>cv.height-(im-Imin)/(Imax-Imin)*cv.height;
const colOf=i=>`hsl(${(i*137.5)%360}, 62%, 60%)`;

let sel = ENTRIES.reduce((b,e,i)=> e.h>ENTRIES[b].h ? i : b, 0);  // default: a high-h class

function draw(){
  ctx.clearRect(0,0,cv.width,cv.height);
  // fundamental domain boundary
  ctx.strokeStyle="rgba(77,163,216,0.35)"; ctx.lineWidth=1.4;
  ctx.beginPath();
  ctx.moveTo(X(-0.5),0); ctx.lineTo(X(-0.5),Y(Math.sqrt(3)/2));
  for(let a=120;a>=60;a-=2){const r=a*Math.PI/180; ctx.lineTo(X(Math.cos(r)),Y(Math.sin(r)));}
  ctx.lineTo(X(0.5),0);
  ctx.stroke();
  // points
  ENTRIES.forEach((e,i)=>{
    const on = i===sel;
    e.pts.forEach(p=>{
      const x=X(p[0]), y=Y(p[1]);
      ctx.beginPath(); ctx.arc(x,y, on?5.5:3.6, 0,7);
      ctx.fillStyle = on ? colOf(i) : `hsla(${(i*137.5)%360},45%,55%,0.45)`;
      ctx.fill();
      if(on){ ctx.strokeStyle="#fff"; ctx.lineWidth=1.6; ctx.stroke(); }
    });
  });
}

function showPanel(){
  const e=ENTRIES[sel];
  let html = `<div style="font-size:1.02rem;color:${INK};margin-bottom:6px;">`
    + `discriminant <b>D = ${e.d}</b> &nbsp;·&nbsp; class number <b>h = ${e.h}</b></div>`;
  html += `<div style="color:${MUT};font-size:.86rem;margin-bottom:10px;">`
    + (e.h===1 ? `one CM point — a single singular modulus`
              : `${e.h} CM points, all sharing this endomorphism ring`) + `</div>`;
  html += `<div style="color:${MUT};font-size:.86rem;">Hilbert class polynomial</div>`;
  html += `<div style="color:${INK};font-size:1.05rem;line-height:1.7;word-break:break-word;margin-top:2px;">`
    + `H<sub>${e.d}</sub>(x) = ${e.poly}</div>`;
  if(e.j0!==null && e.j0!==undefined)
    html += `<div style="color:${ACC};font-size:.95rem;margin-top:8px;">singular modulus &nbsp; j = ${e.j0}</div>`;
  panel.innerHTML=html;
}

cv.addEventListener("pointerdown",e=>{
  const r=cv.getBoundingClientRect();
  const px=(e.clientX-r.left)*cv.width/r.width, py=(e.clientY-r.top)*cv.height/r.height;
  let best=-1, bd=14*14;
  ENTRIES.forEach((en,i)=>en.pts.forEach(p=>{
    const dx=X(p[0])-px, dy=Y(p[1])-py, d2=dx*dx+dy*dy;
    if(d2<bd){bd=d2; best=i;}
  }));
  if(best>=0){ sel=best; draw(); showPanel(); }
});

draw(); showPanel();
</script>
"""


def x0_subgroup_html() -> str:
    """§9.1 applet: a point of X_0(l) = a curve + an order-l subgroup.

    Left panel: the Gamma_0(l) fundamental domain -- l+1 tiles, one per subgroup;
    drag the marker to choose a point of X_0(l), which sets tau AND the subgroup
    at once. Right panel: E = C/<1,tau> with its l^2 torsion points, coloured by
    which of the l+1 cyclic order-l subgroups they lie in; click a point to light
    its subgroup.

    Geometry: tile delta_k = S T^k = (0 -1; 1 k) carries the base subgroup
    <(1,0)> to <(k,1)>, and the point (shape s, subgroup C_k) sits at
    delta_k . s = -1/(s+k); the base tile F is the subgroup C_inf = <(1,0)>.
    """
    return _HEAD + r"""
<div class="panel">
  <div class="modebar">
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin-right:4px;">a point of X₀(ℓ), &nbsp; ℓ =</span>
    <button class="seg ellbtn" data-l="2">2</button>
    <button class="seg ellbtn" data-l="3">3</button>
    <button class="seg ellbtn on" data-l="5">5</button>
    <button class="seg ellbtn" data-l="7">7</button>
  </div>
  <div class="stage">
    <div class="cell">
      <div class="cap">a point of X₀(ℓ) — the Γ₀(ℓ) tiling</div>
      <canvas id="mvX" width="340" height="340"></canvas>
    </div>
    <div class="cell">
      <div class="cap">E[ℓ] and its ℓ+1 order-ℓ subgroups</div>
      <canvas id="mvE" width="340" height="340"></canvas>
    </div>
  </div>
  <div class="info" id="mvInfo"></div>
  <div class="hint" id="mvHint"></div>
</div>
<script>
"use strict";
const ACC="#4da3d8", INK="#d7d9dc", MUT="#9aa4ad", GRID="#3a3f45";

// ---- complex ----
const C=(re,im)=>({re,im});
const cabs2=z=>z.re*z.re+z.im*z.im;
const cmul=(a,b)=>C(a.re*b.re-a.im*b.im, a.re*b.im+a.im*b.re);
const cinv=z=>{const d=cabs2(z);return C(z.re/d,-z.im/d);};
const mob=(a,b,c,d,z)=>cmul(C(a*z.re+b,a*z.im), cinv(C(c*z.re+d,c*z.im))); // (az+b)/(cz+d)
function modinv(x,m){x=((x%m)+m)%m;for(let i=1;i<m;i++)if((x*i)%m===1)return i;return 0;}
function reduceFD(z){let t=C(z.re,z.im);for(let g=0;g<80;g++){const n=Math.round(t.re);t=C(t.re-n,t.im);if(cabs2(t)<1-1e-9){const d=cabs2(t);t=C(-t.re/d,t.im/d);}else break;}return t;}
const inFD=z=>Math.abs(z.re)<=0.5+1e-6 && cabs2(z)>=1-1e-6 && z.im>0;

// ---- state ----
let ell=5;
let s=C(0.20,1.25);   // shape tau in F
let sub=ell;          // subgroup index: 0..ell-1 => C_k=<(k,1)>, ell => C_inf=<(1,0)>
let dragging=false;

const xCv=document.getElementById("mvX"), xCtx=xCv.getContext("2d");
const eCv=document.getElementById("mvE"), eCtx=eCv.getContext("2d");
const info=document.getElementById("mvInfo"), hint=document.getElementById("mvHint");

const nSub=()=>ell+1;
const subHue=i=>`${Math.round(360*i/nSub())}`;
const subColor=(i,on)=>`hsl(${subHue(i)}, ${on?70:38}%, ${on?62:50}%)`;
// centred display offset kk for subgroup index sIdx (0..ell-1): petals fan symmetrically for odd ell
const kkOf=sIdx => (sIdx <= (ell-1)/2 ? sIdx : sIdx-ell);

// ---- X(1) canvas mapping ----
const XW={xmin:-1.7,xmax:1.7,ymin:0.0,ymax:2.7};
const xX=re=>(re-XW.xmin)/(XW.xmax-XW.xmin)*xCv.width;
const xY=im=>xCv.height-(im-XW.ymin)/(XW.ymax-XW.ymin)*xCv.height;
const xInv=(px,py)=>C(XW.xmin+px/xCv.width*(XW.xmax-XW.xmin), XW.ymin+(xCv.height-py)/xCv.height*(XW.ymax-XW.ymin));

function fdPath(ctx, map){ // trace the level-1 fundamental domain boundary via map(z)->[px,py]
  ctx.beginPath();
  let first=true;
  const push=z=>{const p=map(z); if(first){ctx.moveTo(p[0],p[1]);first=false;} else ctx.lineTo(p[0],p[1]);};
  for(let im=XW.ymax; im>=Math.sqrt(3)/2; im-=0.12) push(C(-0.5,im));   // left wall down
  for(let a=120;a>=60;a-=3) push(C(Math.cos(a*Math.PI/180),Math.sin(a*Math.PI/180))); // arc
  for(let im=Math.sqrt(3)/2; im<=XW.ymax; im+=0.12) push(C(0.5,im));   // right wall up
}

function drawX(){
  xCtx.clearRect(0,0,xCv.width,xCv.height);
  xCtx.strokeStyle=GRID; xCtx.lineWidth=1;
  xCtx.beginPath(); xCtx.moveTo(0,xY(0)); xCtx.lineTo(xCv.width,xY(0)); xCtx.stroke();
  // Gamma_0(l) tiling: base F (C_inf) + petals delta_kk . F, coloured by subgroup
  for(let sIdx=0;sIdx<ell;sIdx++){
    const kk=kkOf(sIdx);
    fdPath(xCtx, z=>{const w=mob(0,-1,1,kk,z); return [xX(w.re),xY(w.im)];});
    xCtx.fillStyle=`hsla(${subHue(sIdx)},45%,52%,${sub===sIdx?0.34:0.14})`; xCtx.fill();
    xCtx.strokeStyle=subColor(sIdx, sub===sIdx); xCtx.lineWidth=sub===sIdx?1.8:1; xCtx.stroke();
  }
  // base tile (C_inf)
  fdPath(xCtx, z=>[xX(z.re),xY(z.im)]);
  xCtx.fillStyle=`hsla(${subHue(ell)},45%,52%,${sub===ell?0.30:0.12})`; xCtx.fill();
  xCtx.strokeStyle=subColor(ell, sub===ell); xCtx.lineWidth=sub===ell?1.8:1; xCtx.stroke();
  // marker at delta . s
  const m = (sub===ell) ? s : mob(0,-1,1,kkOf(sub),s);
  xCtx.fillStyle="#fff"; xCtx.beginPath(); xCtx.arc(xX(m.re),xY(m.im),5.5,0,7); xCtx.fill();
  xCtx.strokeStyle=subColor(sub,true); xCtx.lineWidth=2; xCtx.stroke();
}

function drawE(){
  eCtx.clearRect(0,0,eCv.width,eCv.height);
  // fit the parallelogram {x*1 + y*tau} into the canvas
  const pad=26;
  const xs=[0,1,s.re,1+s.re], ys=[0,s.im];
  const xmin=Math.min(...xs), xmax=Math.max(...xs), ymax=Math.max(...ys);
  const sc=Math.min((eCv.width-2*pad)/(xmax-xmin), (eCv.height-2*pad)/ymax);
  const ox=pad-xmin*sc, oy=eCv.height-pad;
  const P=z=>[ox+z.re*sc, oy-z.im*sc];
  const lat=(x,y)=>C(x+y*s.re, y*s.im);   // x*1 + y*tau
  // parallelogram
  const c0=P(lat(0,0)),c1=P(lat(1,0)),c11=P(lat(1,1)),c01=P(lat(0,1));
  eCtx.strokeStyle="#4a4f55"; eCtx.lineWidth=1.5;
  eCtx.beginPath(); eCtx.moveTo(...c0); eCtx.lineTo(...c1); eCtx.lineTo(...c11); eCtx.lineTo(...c01); eCtx.closePath(); eCtx.stroke();
  // selected subgroup drawn as a torus geodesic 0, g, 2g, ... : a straight line
  // of slope g wound onto the torus, wrapping across the edges (rather than
  // jumping back to 0). Uses a minimal-winding generator v/ell.
  {
    const v = (sub===ell) ? [1,0] : [kkOf(sub),1];   // generator step v/ell
    const wind = Math.abs(v[0]) + Math.abs(v[1]);
    const NS = Math.max(240, 70*ell*wind);
    eCtx.strokeStyle=subColor(sub,true); eCtx.lineWidth=1.7;
    eCtx.beginPath();
    let pen=false, px=0, py=0;
    for(let i=0;i<=NS;i++){
      const t=i/NS*ell;
      const x=((t*v[0]/ell)%1+1)%1, y=((t*v[1]/ell)%1+1)%1;
      const q=P(lat(x,y));
      if(pen && (Math.abs(x-px)>0.5 || Math.abs(y-py)>0.5)) pen=false;  // wrapped -> lift pen
      if(!pen){ eCtx.moveTo(q[0],q[1]); pen=true; } else eCtx.lineTo(q[0],q[1]);
      px=x; py=y;
    }
    eCtx.stroke();
  }
  // torsion points
  window._eHit=[];
  for(let a=0;a<ell;a++)for(let b=0;b<ell;b++){
    const p=P(lat(a/ell,b/ell));
    const zero=(a===0&&b===0);
    const sIdx = zero? -1 : (b===0? ell : (a*modinv(b,ell))%ell);
    const on = (sIdx===sub) || zero;
    eCtx.fillStyle = zero ? INK : subColor(sIdx,on);
    eCtx.beginPath(); eCtx.arc(p[0],p[1], on?5:3.4, 0,7); eCtx.fill();
    if(on && !zero){ eCtx.strokeStyle="rgba(255,255,255,0.6)"; eCtx.lineWidth=1; eCtx.stroke(); }
    if(!zero) window._eHit.push([p[0],p[1],sIdx]);
  }
  // origin label
  eCtx.fillStyle=INK; eCtx.font="12px system-ui"; eCtx.fillText("0", P(lat(0,0))[0]-12, P(lat(0,0))[1]+4);
}

function subName(){ return sub===ell ? "⟨1/ℓ⟩" : `⟨(${kkOf(sub)}·1+τ)/ℓ⟩`; }
function render(){
  drawX(); drawE();
  const idx = sub===ell ? (ell+1) : (sub+1);
  info.innerHTML = `E = ℂ/⟨1, τ⟩, &nbsp; τ = ${s.re.toFixed(2)} + ${s.im.toFixed(2)}i `
    + `&nbsp;·&nbsp; subgroup C = <b style="color:${subColor(sub,true)}">${subName()}</b> `
    + `(${idx} of ℓ+1 = ${ell+1}) &nbsp;→&nbsp; a point of X₀(${ell}).`;
  hint.textContent = "drag the marker on the left to choose a point of X₀(ℓ); or click a torsion point on the right";
}

// ---- interaction ----
function evt(cv,e){const r=cv.getBoundingClientRect();return {x:(e.clientX-r.left)*cv.width/r.width, y:(e.clientY-r.top)*cv.height/r.height};}

// right panel: pick a subgroup by clicking a torsion point
eCv.addEventListener("pointerdown",e=>{
  const p=evt(eCv,e); let best=-1,bd=16*16;
  for(const [hx,hy,sIdx] of (window._eHit||[])){const dx=hx-p.x,dy=hy-p.y,d2=dx*dx+dy*dy; if(d2<bd){bd=d2;best=sIdx;}}
  if(best>=0){ sub=best; render(); }
});

// left panel: drag a point of X_0(l); locate its tile -> (tau, subgroup)
function setX0FromPx(p){
  const z=xInv(p.x,p.y);
  // find the tile whose delta^{-1} z lands in F
  let bestSub=ell, bestCand=reduceFD(z), bestScore=inFD(z)?2+z.im:z.im;
  for(let sIdx=0;sIdx<ell;sIdx++){
    const kk=kkOf(sIdx);
    const cand=mob(kk,1,-1,0,z);           // delta_kk^{-1} = (kk 1; -1 0)
    const score=(inFD(cand)?2:0)+cand.im;
    if(score>bestScore){bestScore=score;bestSub=sIdx;bestCand=cand;}
  }
  sub=bestSub; s=reduceFD(bestCand); if(s.im<0.06)s=C(s.re,0.06); render();
}
xCv.addEventListener("pointerdown",e=>{dragging=true; xCv.setPointerCapture(e.pointerId); setX0FromPx(evt(xCv,e)); e.preventDefault();});
xCv.addEventListener("pointermove",e=>{if(dragging)setX0FromPx(evt(xCv,e));});
window.addEventListener("pointerup",()=>{dragging=false;});

document.querySelectorAll(".ellbtn").forEach(b=>b.addEventListener("click",()=>{
  ell=+b.dataset.l; sub=ell;
  document.querySelectorAll(".ellbtn").forEach(x=>x.classList.toggle("on",x===b));
  render();
}));

render();
</script>
"""


# Endomorphism points on X_0(l): CM points where the l-isogeny is an endo
# (j(tau)=j(l tau) -- the Fricke-fixed / same-j locus). Precomputed with
# qfs.qf_x0_endos over |D| <= 80, as [re, im, disc] on the Gamma_0(l) tiling.
_X0_ENDOS = {
  "2": [[0.5,0.5,-4],[0.375,0.3307,-7],[0.25,0.6614,-7],[0.0,0.7071,-8]],
  "3": [[0.5,0.2887,-3],[0.3333,0.4714,-8],[-0.3333,0.4714,-8],[0.1667,0.5528,-11],[-0.1667,0.5528,-11],[0.0,0.5774,-12]],
  "5": [[0.4,0.2,-4],[-0.4,0.2,-4],[0.3,0.3317,-11],[-0.3,0.3317,-11],[0.2,0.4,-16],[-0.2,0.4,-16],[0.1,0.4359,-19],[-0.1,0.4359,-19],[0.0,0.4472,-20],[0.3333,0.1491,-20]],
  "7": [[0.3571,0.1237,-3],[-0.3571,0.1237,-3],[0.25,0.0945,-7],[0.2857,0.2474,-12],[-0.2857,0.2474,-12],[0.2143,0.3113,-19],[-0.2143,0.3113,-19],[0.1429,0.3499,-24],[-0.1429,0.3499,-24],[0.2857,0.1166,-24],[-0.2857,0.1166,-24],[0.0714,0.3712,-27],[-0.0714,0.3712,-27],[0.0,0.378,-28]]
}


def x0_fricke_html() -> str:
    """§9.2 applet: the j-maps and the Fricke involution.

    Left: the Gamma_0(l) tiling (coloured by subgroup, as on the previous tab),
    with a draggable point tau, its Fricke image w(tau) = -1/(l tau) reduced back
    into the tiling, and the interior real locus of j (dashed). Right: the two
    lattices of the isogeny E -> E' -- domain Lambda = <1,tau> and the index-l
    superlattice Lambda' in a second colour. (_X0_ENDOS is kept for the later
    algebraic before/after comparison; not drawn here.)
    """
    return _HEAD + r"""
<div class="panel">
  <div class="modebar">
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin-right:4px;">an isogeny in X₀(ℓ), &nbsp; ℓ =</span>
    <button class="seg ellbtn" data-l="2">2</button>
    <button class="seg ellbtn" data-l="3">3</button>
    <button class="seg ellbtn on" data-l="5">5</button>
    <button class="seg ellbtn" data-l="7">7</button>
  </div>
  <div class="stage">
    <div class="cell">
      <div class="cap">X₀(ℓ): the isogeny τ and its Fricke dual 𝔉ℓ(τ)</div>
      <canvas id="mvX" width="340" height="340"></canvas>
    </div>
    <div class="cell">
      <div class="cap">the isogeny Λ ⊂ Λ′ &nbsp;(domain · superlattice)</div>
      <canvas id="mvE" width="340" height="340"></canvas>
    </div>
  </div>
  <div class="info" id="mvInfo"></div>
  <div class="hint" id="mvHint"></div>
</div>
<script>
"use strict";
const DOM="#4da3d8", SUP="#e0b64f", INK="#d7d9dc", MUT="#9aa4ad", GRID="#3a3f45";

const C=(re,im)=>({re,im});
const cabs2=z=>z.re*z.re+z.im*z.im;
const cmul=(a,b)=>C(a.re*b.re-a.im*b.im, a.re*b.im+a.im*b.re);
const cinv=z=>{const d=cabs2(z);return C(z.re/d,-z.im/d);};
const mob=(a,b,c,d,z)=>cmul(C(a*z.re+b,a*z.im), cinv(C(c*z.re+d,c*z.im)));

let ell=5, s=C(0.20,1.30), sub=0, dragging=false;
const kkOf=sIdx => (sIdx <= (ell-1)/2 ? sIdx : sIdx-ell);
const nSub=()=>ell+1;
const subHue=i=>`${Math.round(360*i/nSub())}`;
const subColor=(i,on)=>`hsl(${subHue(i)}, ${on?70:38}%, ${on?62:50}%)`;

const xCv=document.getElementById("mvX"), xCtx=xCv.getContext("2d");
const eCv=document.getElementById("mvE"), eCtx=eCv.getContext("2d");
const info=document.getElementById("mvInfo"), hint=document.getElementById("mvHint");

// left window: the full Gamma_0(l) fundamental domain (as on the previous tab)
const XW={xmin:-1.7,xmax:1.7,ymin:0.0,ymax:2.7};
const xX=re=>(re-XW.xmin)/(XW.xmax-XW.xmin)*xCv.width;
const xY=im=>xCv.height-(im-XW.ymin)/(XW.ymax-XW.ymin)*xCv.height;
const xInv=(px,py)=>C(XW.xmin+px/xCv.width*(XW.xmax-XW.xmin), XW.ymin+(xCv.height-py)/xCv.height*(XW.ymax-XW.ymin));

function reduceMat(z){
  let t=C(z.re,z.im), g=[[1,0],[0,1]];
  for(let it=0;it<120 && t.im>1e-9;it++){
    const n=Math.round(t.re);
    if(n!==0){ t=C(t.re-n,t.im); g=[[g[0][0]-n*g[1][0],g[0][1]-n*g[1][1]],[g[1][0],g[1][1]]]; }
    if(cabs2(t)<1-1e-9){ const d=cabs2(t); t=C(-t.re/d,t.im/d); g=[[-g[1][0],-g[1][1]],[g[0][0],g[0][1]]]; }
    else break;
  }
  return {s:t,g:g};
}
// FD representative of the X_0(l) point of any w in H: its tile + position
function x0pos(w){
  const r=reduceMat(w), g=r.g;
  const gi=[[g[1][1],-g[0][1]],[-g[1][0],g[0][0]]];   // g^{-1}
  const m=x=>(((x%ell)+ell)%ell);
  let sb=ell;
  if(m(gi[1][0])!==0){ for(let k=0;k<ell;k++) if(m(gi[1][0]*kkOf(k)-gi[1][1])===0){sb=k;break;} }
  const pos = sb===ell ? r.s : mob(0,-1,1,kkOf(sb),r.s);
  return {sub:sb, s:r.s, pos};
}

function fdPath(ctx, map){
  ctx.beginPath(); let first=true;
  const push=z=>{const p=map(z); if(first){ctx.moveTo(p[0],p[1]);first=false;} else ctx.lineTo(p[0],p[1]);};
  for(let im=XW.ymax+0.6; im>=Math.sqrt(3)/2; im-=0.1) push(C(-0.5,im));
  for(let a=120;a>=60;a-=3) push(C(Math.cos(a*Math.PI/180),Math.sin(a*Math.PI/180)));
  for(let im=Math.sqrt(3)/2; im<=XW.ymax+0.6; im+=0.1) push(C(0.5,im));
}

function drawX(){
  xCtx.clearRect(0,0,xCv.width,xCv.height);
  xCtx.strokeStyle=GRID; xCtx.lineWidth=1;
  xCtx.beginPath(); xCtx.moveTo(0,xY(0)); xCtx.lineTo(xCv.width,xY(0)); xCtx.stroke();
  // Gamma_0(l) tiling coloured by subgroup (the boundary = the real locus of j)
  for(let k=0;k<ell;k++){
    fdPath(xCtx, z=>{const w=mob(0,-1,1,kkOf(k),z);return [xX(w.re),xY(w.im)];});
    xCtx.fillStyle=`hsla(${subHue(k)},45%,52%,${sub===k?0.30:0.13})`; xCtx.fill();
    xCtx.strokeStyle=subColor(k, sub===k); xCtx.lineWidth=sub===k?1.7:1; xCtx.stroke();
  }
  fdPath(xCtx, z=>[xX(z.re),xY(z.im)]);
  xCtx.fillStyle=`hsla(${subHue(ell)},45%,52%,${sub===ell?0.28:0.11})`; xCtx.fill();
  xCtx.strokeStyle=subColor(ell, sub===ell); xCtx.lineWidth=sub===ell?1.7:1; xCtx.stroke();
  // interior real locus: delta_k . (imaginary axis) -- j real, but not on the
  // boundary, so drawn dashed to set it apart
  xCtx.strokeStyle="rgba(228,231,235,0.34)"; xCtx.lineWidth=1; xCtx.setLineDash([4,4]);
  const imAxis=mapfn=>{ xCtx.beginPath(); let f=true;
    for(let t=1;t<=9;t+=0.05){ const w=mapfn(C(0,t)); const px=xX(w.re),py=xY(w.im);
      if(f){xCtx.moveTo(px,py);f=false;}else xCtx.lineTo(px,py);} xCtx.stroke(); };
  imAxis(z=>z);
  for(let k=0;k<ell;k++){ const kk=kkOf(k); imAxis(z=>mob(0,-1,1,kk,z)); }
  xCtx.setLineDash([]);
  // tau and its Fricke dual
  const tau = sub===ell ? s : mob(0,-1,1,kkOf(sub),s);
  const w = cmul(C(-1,0), cinv(C(ell*tau.re, ell*tau.im)));   // -1/(l tau)
  const fr = x0pos(w).pos;
  xCtx.fillStyle="rgba(224,182,79,0.6)"; xCtx.beginPath(); xCtx.arc(xX(fr.re),xY(fr.im),6,0,7); xCtx.fill();
  xCtx.strokeStyle=SUP; xCtx.lineWidth=1.6; xCtx.stroke();
  xCtx.fillStyle=INK; xCtx.font="12px system-ui"; xCtx.fillText("𝔉(τ)", xX(fr.re)+9, xY(fr.im)-6);
  xCtx.fillStyle="#fff"; xCtx.beginPath(); xCtx.arc(xX(tau.re),xY(tau.im),6,0,7); xCtx.fill();
  xCtx.strokeStyle=DOM; xCtx.lineWidth=2; xCtx.stroke();
  xCtx.fillStyle=INK; xCtx.fillText("τ", xX(tau.re)+9, xY(tau.im)-6);
}

function drawE(){
  eCtx.clearRect(0,0,eCv.width,eCv.height);
  const gv = sub===ell ? [1/ell,0] : [kkOf(sub)/ell, 1/ell];  // subgroup generator in (a,b) coords
  const cpx=(a,b)=>C(a + b*s.re, b*s.im);
  const pad=14;
  const sc=(eCv.width-2*pad)/3.3;         // fixed scale: ~3 cells across -> denser
  const ox=eCv.width*0.34, oy=eCv.height-pad;
  const P=z=>[ox+z.re*sc, oy-z.im*sc];
  const vis=p=>p[0]>=-2&&p[0]<=eCv.width+2&&p[1]>=-2&&p[1]<=eCv.height+2;
  // superlattice extra points (Lambda' \ Lambda)
  eCtx.fillStyle="rgba(224,182,79,0.9)";
  for(let mm=-3;mm<=6;mm++)for(let nn=-3;nn<=9;nn++)for(let i=1;i<ell;i++){
    const p=P(cpx(mm+i*gv[0], nn+i*gv[1])); if(!vis(p))continue;
    eCtx.beginPath(); eCtx.arc(p[0],p[1],2.0,0,7); eCtx.fill();
  }
  // domain lattice points
  eCtx.fillStyle=DOM;
  for(let mm=-3;mm<=6;mm++)for(let nn=-3;nn<=9;nn++){
    const p=P(cpx(mm,nn)); if(!vis(p))continue;
    eCtx.beginPath(); eCtx.arc(p[0],p[1],3.2,0,7); eCtx.fill();
  }
  // one fundamental parallelogram of Lambda
  const c0=P(cpx(0,0)),c1=P(cpx(1,0)),c11=P(cpx(1,1)),c01=P(cpx(0,1));
  eCtx.strokeStyle="#565b61"; eCtx.lineWidth=1.2;
  eCtx.beginPath(); eCtx.moveTo(...c0);eCtx.lineTo(...c1);eCtx.lineTo(...c11);eCtx.lineTo(...c01);eCtx.closePath(); eCtx.stroke();
}

function subName(){
  if(sub===ell) return "⟨1/ℓ⟩";
  const kk=kkOf(sub);
  const g = kk===0 ? "τ" : (kk>0 ? `${kk}+τ` : `τ−${-kk}`);
  return `⟨(${g})/ℓ⟩`;
}
function render(){
  drawX(); drawE();
  info.innerHTML = `<b style="color:${DOM}">E = ℂ/Λ</b>, Λ = ⟨1, τ⟩ &nbsp;→&nbsp; `
    + `<b style="color:${SUP}">E′ = ℂ/Λ′</b> (degree ℓ, subgroup ${subName()}). `
    + `Fricke sends τ to <b style="color:${SUP}">𝔉ℓ(τ) = −1/(ℓτ)</b> — the dual isogeny E′ → E.`;
  hint.textContent = "drag anywhere on the left to move the isogeny τ; the dashed curves are the interior real locus of j (Re τ = 0 and its Γ₀(ℓ)-images)";
}

function evt(cv,e){const r=cv.getBoundingClientRect();return {x:(e.clientX-r.left)*cv.width/r.width, y:(e.clientY-r.top)*cv.height/r.height};}
function setFromPx(p){ const r=x0pos(xInv(p.x,p.y)); sub=r.sub; s=r.s; if(s.im<0.06)s=C(s.re,0.06); render(); }
xCv.addEventListener("pointerdown",e=>{dragging=true; xCv.setPointerCapture(e.pointerId); setFromPx(evt(xCv,e)); e.preventDefault();});
xCv.addEventListener("pointermove",e=>{if(dragging)setFromPx(evt(xCv,e));});
window.addEventListener("pointerup",()=>{dragging=false;});
document.querySelectorAll(".ellbtn").forEach(b=>b.addEventListener("click",()=>{
  ell=+b.dataset.l; sub=0; s=C(0.20,1.30);
  document.querySelectorAll(".ellbtn").forEach(x=>x.classList.toggle("on",x===b));
  render();
}));

render();
</script>
"""


# Genus-0 j-maps for X_0(l) from data/jcoefs.json: j_l(x) = a1(x)^1 * a3(x)^3 *
# a_{-1}(x)^{-1} (keys = exponents, lists = coeffs low->high). Fricke x -> m/x,
# m = l^{12/(l-1)}.  Only the odd genus-0 primes are tabulated.
_JCOEFS = {
  "3":  {"1": [27, -1], "3": [-3, 1],  "-1": [0, 1]},
  "5":  {"1": [1],      "3": [5, 10, 1], "-1": [0, 1]},
  "7":  {"1": [-49, 13, -1], "3": [1, -5, 1], "-1": [0, 1]},
  "13": {"1": [-13, 5, -1], "3": [1, -19, 20, -7, 1], "-1": [0, 1]},
}
_MVAL = {"3": 729, "5": 125, "7": 49, "13": 13}

# Endomorphism (diagonal) points of X_0(l): the 2l roots of j_l(x) = j_l(m/x).
# ALL lie exactly on the self-dual circle |x| = sqrt(m) (verified to 1e-10;
# l=5 matches the exact values -11±2i, -2±11i, ±5√5, -9±2i√11, 7±2i√19 from
# the x05pic.nb notebook). Precomputed via numpy polyroots; baked so the page
# stays import-light.
_ENDO_PTS = {
  "3": [[-5.0, -26.532998], [23.0, -14.142136], [27.0, 0.0], [23.0, 14.142136],
        [-5.0, 26.532998], [-27.0, 0.0]],
  "5": [[-11.0, -2.0], [-9.0, -6.63325], [-2.0, -11.0], [7.0, -8.717798],
        [11.18034, 0.0], [7.0, 8.717798], [-2.0, 11.0], [-9.0, 6.63325],
        [-11.0, 2.0], [-11.18034, 0.0]],
  "7": [[-5.5, -4.330127], [-1.742641, -6.779617], [2.5, -6.538348],
        [5.5, -4.330127], [6.5, -2.598076], [6.742641, -1.880637], [7.0, 0.0],
        [6.742641, 1.880637], [6.5, 2.598076], [5.5, 4.330127],
        [2.5, 6.538348], [-1.742641, 6.779617], [-5.5, 4.330127], [-7.0, 0.0]],
  "13": [[-3.342329, -1.352344], [-2.598076, -2.5], [-1.5, -3.278719],
         [-0.232051, -3.598076], [1.0, -3.464102], [2.0, -3.0],
         [2.5, -2.598076], [2.598076, -2.5], [2.842329, -2.21837],
         [3.0, -2.0], [3.232051, -1.598076], [3.5, -0.866025],
         [3.605551, 0.0], [3.5, 0.866025], [3.232051, 1.598076], [3.0, 2.0],
         [2.842329, 2.21837], [2.598076, 2.5], [2.5, 2.598076], [2.0, 3.0],
         [1.0, 3.464102], [-0.232051, 3.598076], [-1.5, 3.278719],
         [-2.598076, 2.5], [-3.342329, 1.352344], [-3.605551, 0.0]],
}



# Shared JS for the two genus-0 X_0(l) applets below: complex arithmetic, the
# j_l colouring, and makeDial() -- a closure that owns one canvas and draws the
# ambient-plane picture (dim backdrop, Belyi-weighted real locus, self-dual
# circle, endomorphism points, cusp). Both applets call dial.drawBase() and
# layer their own interaction on top.
_G0_LIB = r"""
const DOM="#4da3d8", SUP="#e0b64f", INK="#d7d9dc", MUT="#9aa4ad";
const C=(re,im)=>({re,im});
const cadd=(a,b)=>C(a.re+b.re,a.im+b.im);
const cmulr=(z,r)=>C(z.re*r,z.im*r);
const cmulc=(a,b)=>C(a.re*b.re-a.im*b.im, a.re*b.im+a.im*b.re);
const cabs2=z=>z.re*z.re+z.im*z.im;
const cinv=z=>{const d=cabs2(z);return C(z.re/d,-z.im/d);};
function peval(coeffs,x){let r=C(0,0),xp=C(1,0);for(const c of coeffs){r=cadd(r,cmulr(xp,c));xp=cmulc(xp,x);}return r;}
function cpow(z,n){if(n===0)return C(1,0); if(n<0)return cinv(cpow(z,-n)); let r=C(1,0);for(let i=0;i<n;i++)r=cmulc(r,z);return r;}
function hsv2rgb(h,s,v){const i=Math.floor(h*6),f=h*6-i,p=v*(1-s),q=v*(1-f*s),t=v*(1-(1-f)*s);
  let r,g,b; switch(((i%6)+6)%6){case 0:r=v;g=t;b=p;break;case 1:r=q;g=v;b=p;break;case 2:r=p;g=v;b=t;break;case 3:r=p;g=q;b=v;break;case 4:r=t;g=p;b=v;break;default:r=v;g=p;b=q;} return [r,g,b];}

function makeDial(cvId, JCO, MVAL, ENDO){
  let ell="5";
  const cv=document.getElementById(cvId), ctx=cv.getContext("2d");
  const rad=()=>2.2*Math.sqrt(MVAL[ell]);
  const X=re=>(re+rad())/(2*rad())*cv.width;
  const Y=im=>(rad()-im)/(2*rad())*cv.height;
  function jmap(x){const e=JCO[ell]; let r=C(1,0); for(const k in e){ r=cmulc(r,cpow(peval(e[k],x),+k)); } return r;}
  function jclass(j){ return (j.re>=0 && j.re<=1728) ? 1 : (j.re<0 ? 0 : 2); }

  // ambient backdrop: the j_l-colouring of the plane, dimmed (it is scenery)
  const bgCache={};
  function bg(){
    if(bgCache[ell])return bgCache[ell];
    const W=cv.width,H=cv.height,img=ctx.createImageData(W,H),R=rad();
    for(let py=0;py<H;py++)for(let px=0;px<W;px++){
      const j=jmap(C(-R+px/W*2*R, R-py/H*2*R));
      const mag=Math.hypot(j.re,j.im);
      let hue=(Math.atan2(j.im,j.re)/(2*Math.PI)); hue=((hue%1)+1)%1;
      const val=0.22+0.26*(0.5+0.5*Math.sin(Math.log(mag+1e-12)*1.7));
      const rgb=hsv2rgb(hue,0.34,Math.max(0,Math.min(1,val))), i=(py*W+px)*4;
      img.data[i]=rgb[0]*255; img.data[i+1]=rgb[1]*255; img.data[i+2]=rgb[2]*255; img.data[i+3]=255;
    }
    bgCache[ell]=img; return img;
  }

  // real locus of j_l via marching squares on Im j = 0, Belyi-weighted:
  // class 1 = preimage of [0,1728] (the dessin, prominent), 0 = j<0, 2 = j>1728
  const belyiCache={};
  function belyi(){
    if(belyiCache[ell])return belyiCache[ell];
    const N=150, R=rad(), h=2*R/N, W=N+1, im=new Float64Array(W*W);
    for(let iy=0;iy<W;iy++)for(let ix=0;ix<W;ix++)
      im[iy*W+ix]=jmap(C(-R+ix*h, R-iy*h)).im;
    const segs={0:[],1:[],2:[]};
    const px=ix=>-R+ix*h, py=iy=>R-iy*h;
    for(let iy=0;iy<N;iy++){
      const y0=py(iy), y1=py(iy+1);
      if(y0>=-1e-12 && y1<=1e-12) continue;          // the axis itself: drawn separately
      for(let ix=0;ix<N;ix++){
        const a=im[iy*W+ix], b=im[iy*W+ix+1], c=im[(iy+1)*W+ix+1], d=im[(iy+1)*W+ix];
        const pts=[];
        if(a*b<0) pts.push([px(ix)+h*a/(a-b), y0]);
        if(b*c<0) pts.push([px(ix+1), y0+(y1-y0)*b/(b-c)]);
        if(d*c<0) pts.push([px(ix)+h*d/(d-c), y1]);
        if(a*d<0) pts.push([px(ix), y0+(y1-y0)*a/(a-d)]);
        for(let k=0;k+1<pts.length;k+=2){
          const mx=(pts[k][0]+pts[k+1][0])/2, my=(pts[k][1]+pts[k+1][1])/2;
          segs[jclass(jmap(C(mx,my)))].push([pts[k][0],pts[k][1],pts[k+1][0],pts[k+1][1]]);
        }
      }
    }
    belyiCache[ell]=segs; return segs;
  }
  // same weighting along the real axis (the actual X_0(l))
  const axisCache={};
  function axisRuns(){
    if(axisCache[ell])return axisCache[ell];
    const R=rad(), M=420, runs=[]; let cur=null;
    for(let i=0;i<=M;i++){
      const v=-R+2*R*i/M;
      const cls=Math.abs(v)<1e-9?2:jclass(jmap(C(v,0)));
      if(cur && cur.c===cls) cur.b=v; else {cur={a:v,b:v,c:cls}; runs.push(cur);}
    }
    axisCache[ell]=runs; return runs;
  }

  function drawBase(){
    ctx.putImageData(bg(),0,0);
    const R=rad(), m=MVAL[ell];
    // imaginary axis, faint
    ctx.strokeStyle="rgba(255,255,255,0.12)"; ctx.lineWidth=1;
    ctx.beginPath(); ctx.moveTo(X(0),0); ctx.lineTo(X(0),cv.height); ctx.stroke();
    // Belyi-weighted real locus: faint classes first, dessin on top
    const segs=belyi();
    const styles={0:["rgba(150,185,255,0.30)",1],2:["rgba(255,255,255,0.22)",1],1:["#ffd47f",2.2]};
    for(const cls of [0,2,1]){
      ctx.strokeStyle=styles[cls][0]; ctx.lineWidth=styles[cls][1]; ctx.beginPath();
      for(const s of segs[cls]){ ctx.moveTo(X(s[0]),Y(s[1])); ctx.lineTo(X(s[2]),Y(s[3])); }
      ctx.stroke();
    }
    // the real line, same weighting but heavier
    for(const r of axisRuns()){
      const st={0:["rgba(150,185,255,0.55)",2],2:["rgba(255,255,255,0.45)",2],1:["#ffd47f",3.5]}[r.c];
      ctx.strokeStyle=st[0]; ctx.lineWidth=st[1];
      ctx.beginPath(); ctx.moveTo(X(r.a),Y(0)); ctx.lineTo(X(r.b),Y(0)); ctx.stroke();
    }
    // self-dual circle |x| = sqrt(m)
    ctx.strokeStyle="rgba(255,255,255,0.5)"; ctx.setLineDash([5,4]); ctx.lineWidth=1.2;
    ctx.beginPath(); ctx.arc(X(0),Y(0), Math.sqrt(m)/(2*R)*cv.width, 0, 7); ctx.stroke(); ctx.setLineDash([]);
    // the 2l endomorphism (CM) points -- all on the circle
    ctx.fillStyle="#ef6f6f";
    for(const p of ENDO[ell]){ ctx.beginPath(); ctx.arc(X(p[0]),Y(p[1]),3,0,7); ctx.fill(); }
    // cusp x=0
    ctx.fillStyle=INK; ctx.beginPath(); ctx.arc(X(0),Y(0),3,0,7); ctx.fill();
    ctx.font="12px system-ui"; ctx.fillText("0", X(0)+6, Y(0)+14);
  }

  return {cv, ctx, X, Y, rad, jmap, drawBase,
          setEll(l){ell=l;}, getEll(){return ell;}, m(){return MVAL[ell];}};
}
"""


def genus0_plane_html() -> str:
    """S9.3 applet A: the algebraic model X_0(l) = P^1 as a gallery of pictures.

    One canvas per selected l in {3,5,7,13}: the x-plane dimly domain-coloured
    by j_l (ambient space only), the real locus of j_l weighted Belyi-style
    (preimage of [0,1728] -- the dessin -- bright; preimages of j<=0 and
    j>=1728 faint), the 2l endomorphism points on the self-dual circle
    |x| = sqrt(m), and the real x-line = the locus we care about. No slider --
    the interactive l=5 walk lives in the companion applet below."""
    import json as _j
    return _HEAD + r"""
<div class="panel">
  <div class="modebar">
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin-right:4px;">algebraic X₀(ℓ) = ℙ¹, &nbsp; ℓ =</span>
    <button class="seg ellbtn" data-l="3">3</button>
    <button class="seg ellbtn on" data-l="5">5</button>
    <button class="seg ellbtn" data-l="7">7</button>
    <button class="seg ellbtn" data-l="13">13</button>
  </div>
  <div class="stage" style="grid-template-columns: 420px;">
    <div class="cell">
      <div class="cap">the ambient x-plane; X₀(ℓ) is the real line inside it</div>
      <canvas id="gp" width="420" height="420"></canvas>
      <div class="hint" style="margin-top:6px;">bright locus: j<sub>ℓ</sub> ∈ [0,1728] (the dessin) · faint blue: j ≤ 0 · faint white: j ≥ 1728<br>red dots: the 2ℓ endomorphism points — all on the self-dual circle |x| = √m (dashed)</div>
    </div>
  </div>
</div>
<script>
"use strict";
const JCO=""" + _j.dumps(_JCOEFS) + r""", MVAL=""" + _j.dumps(_MVAL) + r""", ENDO=""" + _j.dumps(_ENDO_PTS) + r""";
""" + _G0_LIB + r"""
const dial=makeDial("gp", JCO, MVAL, ENDO);
document.querySelectorAll(".ellbtn").forEach(b=>b.addEventListener("click",()=>{
  dial.setEll(b.dataset.l);
  document.querySelectorAll(".ellbtn").forEach(z=>z.classList.toggle("on",z===b));
  dial.drawBase();
}));
dial.drawBase();
</script>
"""


def x05_isogeny_html() -> str:
    """S9.3 applet B: walk the real points of X_0(5) and watch the 5-isogeny.

    l = 5 only. Drag x along the real line; the universal-curve section
    t(x) = (11+x-sqrt(x^2+22x+125))/2 (radicand positive on all of R) into the
    Tate form E_t: y^2+(1-t)xy-ty = x^3-tx^2 recovers the DOMAIN curve with its
    kernel <(0,0)> = {(0,0),(0,t),(t,t^2),(t,0)}; the same recipe at the Fricke
    partner 125/x recovers the CODOMAIN E/C. Marker snaps to the self-dual
    points +-sqrt(125), where E = E/C (CM). From x05pic.nb."""
    import json as _j
    return _HEAD + r"""
<div class="panel">
  <div class="stage" style="grid-template-columns: 380px 1fr;">
    <div class="cell">
      <div class="cap">X₀(5): drag x along the real line</div>
      <canvas id="gd" width="380" height="380"></canvas>
      <div class="hint" style="margin-top:6px;">gold = Fricke partner m/x, &nbsp;m = 125</div>
    </div>
    <div class="cell" style="justify-content:flex-start;">
      <div class="cap" style="text-align:left;">the 5-isogeny at x</div>
      <div id="gdpanel" style="border:1px solid #23272c;border-radius:8px;padding:12px 14px;font-size:.93rem;line-height:1.65;"></div>
      <div style="margin-top:10px;">
        <div style="display:flex;gap:10px;">
          <div><canvas id="gcd" width="180" height="168" style="border:1px solid #23272c;border-radius:8px;"></canvas>
               <div id="gcdcap" class="hint" style="text-align:center;margin-top:3px;"></div></div>
          <div><canvas id="gcc" width="180" height="168" style="border:1px solid #23272c;border-radius:8px;"></canvas>
               <div id="gcccap" class="hint" style="text-align:center;margin-top:3px;"></div></div>
        </div>
        <div class="hint" style="text-align:left;margin-top:6px;">E<sub>t</sub>: y² + (1−t)xy − ty = x³ − tx², &nbsp;t(x) = ½(11 + x − √(x²+22x+125)); dots = the kernel ⟨(0,0)⟩.</div>
      </div>
    </div>
  </div>
</div>
<script>
"use strict";
const JCO=""" + _j.dumps(_JCOEFS) + r""", MVAL=""" + _j.dumps(_MVAL) + r""", ENDO=""" + _j.dumps(_ENDO_PTS) + r""";
""" + _G0_LIB + r"""
const dial=makeDial("gd", JCO, MVAL, ENDO);   // ell stays "5"
const panel=document.getElementById("gdpanel");
const cvD=document.getElementById("gcd"), cvC=document.getElementById("gcc");
const capD=document.getElementById("gcdcap"), capC=document.getElementById("gcccap");
let xr=-6.0, dragging=false;

function fmtR(v){ if(!isFinite(v)) return "∞";
  return Math.abs(v)>=1e6?v.toExponential(2):(Math.abs(v)>=100?v.toFixed(0):v.toFixed(2)); }

function draw(){
  dial.drawBase();
  const ctx=dial.ctx, X=dial.X, Y=dial.Y, m=dial.m(), R=dial.rad();
  const fr=m/xr;
  if(Math.abs(fr)<=R){
    ctx.fillStyle=SUP; ctx.beginPath(); ctx.arc(X(fr),Y(0),6,0,7); ctx.fill();
    ctx.strokeStyle="#000"; ctx.lineWidth=1; ctx.stroke();
    ctx.fillStyle=INK; ctx.fillText("m/x", X(fr)+8, Y(0)-8);
  }
  ctx.fillStyle="#fff"; ctx.beginPath(); ctx.arc(X(xr),Y(0),6,0,7); ctx.fill();
  ctx.strokeStyle=DOM; ctx.lineWidth=2; ctx.stroke();
  ctx.fillStyle=INK; ctx.fillText("x", X(xr)+8, Y(0)-8);
}

// section X_0(5) -> X_1(5): t(x) solves t-11-1/t = x, radicand always > 0 on R
const tsec=v=>(11+v-Math.sqrt(v*v+22*v+125))/2;
function drawCurve(cc,t,stroke){
  const g=cc.getContext("2d"), Wp=cc.width, Hp=cc.height;
  g.clearRect(0,0,Wp,Hp);
  const a1=1-t, a2=-t, a3=-t;
  const q=xx=>a1*xx+a3, F=xx=>q(xx)*q(xx)+4*(xx*xx*xx+a2*xx*xx);
  // real roots of the branch cubic F (Cauchy bound + scan/bisect)
  const co=[a3*a3, 2*a1*a3, a1*a1+4*a2, 4];
  const B=1+Math.max(Math.abs(co[0]),Math.abs(co[1]),Math.abs(co[2]))/4;
  const roots=[]; let pv=F(-B);
  for(let i=1;i<=400;i++){ const xv=-B+2*B*i/400, v=F(xv);
    if(pv===0) roots.push(-B+2*B*(i-1)/400);
    if(pv*v<0){ let lo=-B+2*B*(i-1)/400, hi=xv;
      for(let k=0;k<40;k++){ const mid=(lo+hi)/2; (F(lo)*F(mid)<=0)?hi=mid:lo=mid; }
      roots.push((lo+hi)/2); }
    pv=v; }
  const ker=[[0,0],[0,t],[t,t*t],[t,0]];
  let lo=Math.min(...roots,0,t), hi=Math.max(...roots,0,t);
  const span=Math.max(hi-lo,1e-3); lo-=0.12*span; hi+=0.65*span;
  // sample both branches, collect y-range
  const NS=260, ys=[], pts=[];
  for(let i=0;i<=NS;i++){ const xv=lo+(hi-lo)*i/NS, f=F(xv);
    if(f>=0){ const s=Math.sqrt(f), y1=(-q(xv)+s)/2, y2=(-q(xv)-s)/2;
      pts.push([xv,y1,y2]); ys.push(y1,y2); } else pts.push(null); }
  for(const k of ker) ys.push(k[1]);
  let ylo=Math.min(...ys), yhi=Math.max(...ys);
  const ypad=0.1*Math.max(yhi-ylo,1e-3); ylo-=ypad; yhi+=ypad;
  const PX=xx=>(xx-lo)/(hi-lo)*(Wp-12)+6, PY=yy=>(yhi-yy)/(yhi-ylo)*(Hp-12)+6;
  g.strokeStyle="rgba(255,255,255,0.10)"; g.lineWidth=1;
  g.beginPath(); g.moveTo(PX(lo),PY(0)); g.lineTo(PX(hi),PY(0));
  g.moveTo(PX(0),PY(ylo)); g.lineTo(PX(0),PY(yhi)); g.stroke();
  g.strokeStyle=stroke; g.lineWidth=1.8;
  for(const br of [1,2]){ g.beginPath(); let pen=false;
    for(const p of pts){ if(!p){pen=false;continue;}
      const px=PX(p[0]), py=PY(p[br]);
      pen?g.lineTo(px,py):g.moveTo(px,py); pen=true; }
    g.stroke(); }
  g.fillStyle="#ef6f6f";
  for(const k of ker){ g.beginPath(); g.arc(PX(k[0]),PY(k[1]),3,0,7); g.fill(); }
}

function render(){
  const m=dial.m(), sm=Math.sqrt(m);
  if(Math.abs(Math.abs(xr)-sm)<0.015*dial.rad()) xr=(xr<0?-sm:sm);  // snap to +-sqrt(m)
  draw();
  const fr=m/xr, jx=dial.jmap(C(xr,0)).re, jfr=dial.jmap(C(fr,0)).re;
  const selfdual=Math.abs(Math.abs(xr)-sm)<1e-9;
  let html =
    `<div style="color:${MUT};font-size:.85rem;">x = ${fmtR(xr)} &nbsp;·&nbsp; m/x = ${fmtR(fr)}</div>`
    + `<div style="margin-top:8px;"><b style="color:${DOM}">domain</b> &nbsp;j₅(x) = ${fmtR(jx)}</div>`
    + `<div style="margin-top:4px;"><b style="color:${SUP}">codomain</b> &nbsp;j₅(m/x) = ${fmtR(jfr)}</div>`;
  if(selfdual) html += `<div style="color:#ef6f6f;font-size:.85rem;margin-top:6px;">x = ±√m: self-dual — E ≅ E/C, a real curve with a 5-endomorphism (CM).</div>`;
  panel.innerHTML=html;
  const t=tsec(xr), tf=tsec(fr);
  drawCurve(cvD,t,DOM); drawCurve(cvC,tf,SUP);
  capD.innerHTML=`domain: t = ${fmtR(t)}`;
  capC.innerHTML=`codomain: t′ = ${fmtR(tf)}`;
}

function evt(e){const cv=dial.cv, r=cv.getBoundingClientRect();
  let v=-dial.rad()+(e.clientX-r.left)/r.width*2*dial.rad();
  const eps=0.02*dial.rad(); if(Math.abs(v)<eps) v=(v<0?-eps:eps);   // keep off the cusp x=0
  return v;}
dial.cv.addEventListener("pointerdown",e=>{dragging=true; dial.cv.setPointerCapture(e.pointerId); xr=evt(e); render(); e.preventDefault();});
dial.cv.addEventListener("pointermove",e=>{if(dragging){xr=evt(e); render();}});
window.addEventListener("pointerup",()=>{dragging=false;});

render();
</script>
"""
