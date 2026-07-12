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
