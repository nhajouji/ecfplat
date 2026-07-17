"""Canvas applets for the Background page's basics chapters (§1-§3).

Replacements for the original matplotlib/plotly + slider applets, in the same
canvas-in-components.html style as modular_viz (§7-§9): drag/click interfaces,
no Streamlit reruns, no matplotlib mathtext (whose pyparsing parser is not
thread-safe under concurrent reruns).
"""

from modular_viz import _HEAD


def torus_group_html() -> str:
    """S1 applet: the group law on C/Lambda is literally vector addition.

    Replaces the tau sliders + plotly click-grid + 'compute sum' button. Drag
    the tau corner to reshape the lattice; click to place z1 then z2 (drag
    either afterwards; a third click clears). The sum z1+z2 updates LIVE as
    the parallelogram construction, and when it lands outside the fundamental
    domain a gold arrow wraps it back in -- the torus identification."""
    return _HEAD + r"""
<div class="panel">
  <div class="stage" style="grid-template-columns: 1fr;">
    <div class="cell">
      <div class="cap">the group law on ℂ/Λ — vector addition mod Λ</div>
      <canvas id="tg" width="640" height="430"></canvas>
      <div id="tgout" style="margin-top:8px;font-size:.93rem;"></div>
      <div class="hint" style="margin-top:4px;">drag the <b style="color:#e0b64f">τ corner</b> to reshape the lattice · click to place z₁, z₂ (drag to move; a third click clears) · the sum updates live, wrapping back into the parallelogram when it leaves</div>
    </div>
  </div>
</div>
<script>
"use strict";
const INK="#d7d9dc", MUT="#9aa4ad", DOM="#4da3d8", SUP="#e0b64f", RED="#ef6f6f";
const cv=document.getElementById("tg"), ctx=cv.getContext("2d");
const out=document.getElementById("tgout");
let tau={x:0.2,y:1.2};
let pts=[];                       // 0..2 points as fractional (s,t) in [0,1)^2
let drag=null;                    // 'tau' | 0 | 1 | null

// view: fit the fundamental parallelogram + margin, aspect-true
let VX0,VX1,VY0,VY1;
function fitView(){
  const xs=[0,1,1+tau.x,tau.x], ys=[0,0,tau.y,tau.y];
  const m=0.62;
  let x0=Math.min(...xs)-m, x1=Math.max(...xs)+m,
      y0=Math.min(...ys)-m, y1=Math.max(...ys)+m;
  const need=(x1-x0)/(y1-y0), have=cv.width/cv.height;
  if(need<have){const c=(x0+x1)/2,w=(y1-y0)*have; x0=c-w/2; x1=c+w/2;}
  else {const c=(y0+y1)/2,h=(x1-x0)/have; y0=c-h/2; y1=c+h/2;}
  VX0=x0; VX1=x1; VY0=y0; VY1=y1;
}
const PX=x=>(x-VX0)/(VX1-VX0)*cv.width;
const PY=y=>(VY1-y)/(VY1-VY0)*cv.height;
const toXY=(s,t)=>({x:s+t*tau.x, y:t*tau.y});
function toST(x,y){const t=y/tau.y; return {s:x-t*tau.x, t};}

function arrow(x1,y1,x2,y2,col,w){
  ctx.strokeStyle=col; ctx.fillStyle=col; ctx.lineWidth=w;
  ctx.beginPath(); ctx.moveTo(PX(x1),PY(y1)); ctx.lineTo(PX(x2),PY(y2)); ctx.stroke();
  const a=Math.atan2(PY(y2)-PY(y1),PX(x2)-PX(x1));
  ctx.beginPath(); ctx.moveTo(PX(x2),PY(y2));
  ctx.lineTo(PX(x2)-9*Math.cos(a-0.42),PY(y2)-9*Math.sin(a-0.42));
  ctx.lineTo(PX(x2)-9*Math.cos(a+0.42),PY(y2)-9*Math.sin(a+0.42));
  ctx.closePath(); ctx.fill();
}
function cell(dm,dn,fill,stroke,lw){
  const vs=[[0,0],[1,0],[1+tau.x,tau.y],[tau.x,tau.y]];
  ctx.beginPath();
  vs.forEach((v,i)=>{const x=v[0]+dm+dn*tau.x, y=v[1]+dn*tau.y;
    i?ctx.lineTo(PX(x),PY(y)):ctx.moveTo(PX(x),PY(y));});
  ctx.closePath();
  if(fill){ctx.fillStyle=fill; ctx.fill();}
  if(stroke){ctx.strokeStyle=stroke; ctx.lineWidth=lw; ctx.stroke();}
}
function dot(x,y,r,fill,ring){
  ctx.fillStyle=fill; ctx.beginPath(); ctx.arc(PX(x),PY(y),r,0,7); ctx.fill();
  if(ring){ctx.strokeStyle=ring; ctx.lineWidth=2; ctx.stroke();}
}
const f2=v=>v.toFixed(2);

function draw(){
  fitView();
  ctx.clearRect(0,0,cv.width,cv.height);
  // neighbouring lattice translates, then the fundamental cell
  for(let dm=-2;dm<=2;dm++)for(let dn=-2;dn<=2;dn++)
    if(dm||dn) cell(dm,dn,"rgba(255,255,255,0.035)","rgba(255,255,255,0.13)",0.8);
  cell(0,0,"rgba(77,163,216,0.14)",DOM,1.8);
  // lattice points
  ctx.fillStyle="rgba(255,255,255,0.5)";
  for(let dm=-3;dm<=3;dm++)for(let dn=-3;dn<=3;dn++){
    const x=dm+dn*tau.x, y=dn*tau.y;
    if(x>VX0&&x<VX1&&y>VY0&&y<VY1){ctx.beginPath();ctx.arc(PX(x),PY(y),2,0,7);ctx.fill();}
  }
  ctx.font="12px system-ui";
  // identity + tau corner (draggable)
  dot(0,0,4,"#9aa860"); ctx.fillStyle="#c6d38f"; ctx.fillText("0", PX(0)-14, PY(0)+4);
  dot(tau.x,tau.y,6,SUP,"#000"); ctx.fillStyle=INK; ctx.fillText("τ", PX(tau.x)-16, PY(tau.y)-6);
  ctx.fillStyle=MUT; ctx.fillText("1", PX(1)+6, PY(0)+14);
  // the two points and the live sum
  const P=pts.map(p=>toXY(p.s,p.t));
  if(pts.length===2){
    const rs=pts[0].s+pts[1].s, rt=pts[0].t+pts[1].t;
    const raw=toXY(rs,rt), wrapped=(rs>=1||rt>=1),
          red=toXY(((rs%1)+1)%1,((rt%1)+1)%1);
    // parallelogram construction, dotted
    ctx.setLineDash([4,4]); ctx.strokeStyle="rgba(255,255,255,0.5)"; ctx.lineWidth=1.2;
    ctx.beginPath(); ctx.moveTo(PX(P[0].x),PY(P[0].y)); ctx.lineTo(PX(raw.x),PY(raw.y));
    ctx.lineTo(PX(P[1].x),PY(P[1].y)); ctx.stroke(); ctx.setLineDash([]);
    if(wrapped){
      dot(raw.x,raw.y,5,"rgba(224,182,79,0.35)");
      arrow(raw.x,raw.y,red.x,red.y,"rgba(224,182,79,0.8)",1.4);
    }
    dot(red.x,red.y,6,SUP); ctx.fillStyle=SUP;
    ctx.fillText("z₁+z₂", PX(red.x)+8, PY(red.y)-6);
  }
  if(P[0]){arrow(0,0,P[0].x,P[0].y,DOM,1.6); dot(P[0].x,P[0].y,6,DOM);
           ctx.fillStyle=DOM; ctx.fillText("z₁", PX(P[0].x)+8, PY(P[0].y)-6);}
  if(P[1]){arrow(0,0,P[1].x,P[1].y,RED,1.6); dot(P[1].x,P[1].y,6,RED);
           ctx.fillStyle=RED; ctx.fillText("z₂", PX(P[1].x)+8, PY(P[1].y)-6);}
}
function render(){
  draw();
  let s=`τ = ${f2(tau.x)} ${tau.y<0?"−":"+"} ${f2(Math.abs(tau.y))}i`;
  pts.forEach((p,i)=>{s+=` &nbsp;·&nbsp; <b style="color:${i?RED:DOM}">z${i?"₂":"₁"}</b> (s,t) = (${f2(p.s)}, ${f2(p.t)})`;});
  if(pts.length===2){
    const rs=pts[0].s+pts[1].s, rt=pts[0].t+pts[1].t;
    s+=` &nbsp;·&nbsp; <b style="color:${SUP}">z₁+z₂</b> (${f2(((rs%1)+1)%1)}, ${f2(((rt%1)+1)%1)})`;
    if(rs>=1||rt>=1) s+=` <span style="color:${MUT}">(wrapped)</span>`;
  }
  out.innerHTML=s;
}
function evXY(e){const b=cv.getBoundingClientRect();
  return {x:VX0+(e.clientX-b.left)/b.width*(VX1-VX0),
          y:VY1-(e.clientY-b.top)/b.height*(VY1-VY0)};}
function near(p,x,y,px){const dx=PX(p.x)-PX(x),dy=PY(p.y)-PY(y);
  return dx*dx+dy*dy<px*px;}
cv.addEventListener("pointerdown",e=>{
  const q=evXY(e); cv.setPointerCapture(e.pointerId);
  if(near(tau,q.x,q.y,14)){drag="tau";}
  else{
    const P=pts.map(p=>toXY(p.s,p.t));
    const hit=P.findIndex(p=>near(p,q.x,q.y,14));
    if(hit>=0) drag=hit;
    else if(pts.length>=2){pts=[]; drag=null;}
    else{const st=toST(q.x,q.y);
         pts.push({s:((st.s%1)+1)%1, t:((st.t%1)+1)%1}); drag=pts.length-1;}
  }
  render(); e.preventDefault();
});
cv.addEventListener("pointermove",e=>{
  if(drag===null)return;
  const q=evXY(e);
  if(drag==="tau"){tau.x=Math.max(-1.2,Math.min(1.6,q.x)); tau.y=Math.max(0.25,Math.min(3,q.y));}
  else{const st=toST(q.x,q.y); pts[drag]={s:((st.s%1)+1)%1, t:((st.t%1)+1)%1};}
  render();
});
window.addEventListener("pointerup",()=>{drag=null;});
render();
</script>
"""


def isogeny_kernel_html() -> str:
    """S3 applet: isogenies are determined by their kernels.

    Replaces the tau/degree sliders + two plotly on_select panels. Left: the
    fundamental domain of Lambda with the d-torsion grid E[d]; drag the gold
    tau corner to reshape the lattice, pick d with the buttons, click a
    nonzero torsion point to choose the kernel generator P. Right: the
    codomain lattice Lambda' = Lambda + ZP updates instantly -- Lambda gray,
    the cosets P adds gold, and a Gaussian-reduced fundamental cell of
    Lambda' in green with 1/d the area."""
    return _HEAD + r"""
<div class="panel">
  <div class="modebar">
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin-right:4px;">kernel C = ⟨P⟩ ⊂ E[d], &nbsp; d =</span>
    <button class="seg dbtn" data-d="2">2</button>
    <button class="seg dbtn on" data-d="3">3</button>
    <button class="seg dbtn" data-d="5">5</button>
    <button class="seg dbtn" data-d="7">7</button>
  </div>
  <div class="stage" style="grid-template-columns: 1fr 1fr;">
    <div class="cell">
      <div class="cap">domain E = ℂ/Λ &nbsp;(grid = E[d]; click P, drag τ)</div>
      <canvas id="ika" width="330" height="330"></canvas>
    </div>
    <div class="cell">
      <div class="cap">codomain E/C = ℂ/Λ′, &nbsp;Λ′ = Λ + ℤP</div>
      <canvas id="ikb" width="330" height="330"></canvas>
    </div>
  </div>
  <div id="ikout" style="margin-top:8px;font-size:.93rem;"></div>
  <div class="hint" style="margin-top:4px;">red = the kernel ⟨P⟩ (labels kP) · right: Λ gray, cosets added by C gold, green = fundamental cell of Λ′ (area 1/d of Λ's)</div>
</div>
<script>
"use strict";
const INK="#d7d9dc", MUT="#9aa4ad", DOM="#4da3d8", SUP="#e0b64f", RED="#ef6f6f", GRN="#69b382";
const cva=document.getElementById("ika"), cvb=document.getElementById("ikb");
const out=document.getElementById("ikout");
let tau={x:0.30,y:1.00}, d=3, ker=[1,1], dragTau=false;

const pos=(cu,cv)=>({x:cu+cv*tau.x, y:cv*tau.y});
// shared square view box around Lambda's fundamental domain
let CX,CY,HS;
function fitView(){
  const cs=[pos(0,0),pos(1,0),pos(1,1),pos(0,1)];
  const xs=cs.map(c=>c.x), ys=cs.map(c=>c.y);
  CX=(Math.min(...xs)+Math.max(...xs))/2; CY=(Math.min(...ys)+Math.max(...ys))/2;
  HS=Math.max(Math.max(...xs)-Math.min(...xs), Math.max(...ys)-Math.min(...ys))/2+0.25;
}
const PX=(cv,x)=>(x-(CX-HS))/(2*HS)*cv.width;
const PY=(cv,y)=>((CY+HS)-y)/(2*HS)*cv.height;

function egcd(a,b){ if(!b) return [a,1,0]; const [g,x,y]=egcd(b,a%b); return [g,y,x-Math.floor(a/b)*y]; }
function lprimeCell(){                     // reduced basis of Lambda' = Lambda + ZP
  const [a0,b0]=ker; let V1,V2;
  if(a0%d!==0){ const t=((egcd(d,a0)[2])%d+d)%d; V1=[1,(t*b0)%d]; V2=[0,d]; }
  else { V1=[d,0]; V2=[0,1]; }
  let f1=pos(V1[0]/d,V1[1]/d), f2=pos(V2[0]/d,V2[1]/d);
  const dp=(u,v)=>u.x*v.x+u.y*v.y;
  for(let i=0;i<50;i++){
    if(dp(f2,f2)<dp(f1,f1)) [f1,f2]=[f2,f1];
    const m=Math.round(dp(f1,f2)/dp(f1,f1));
    if(!m) break;
    f2={x:f2.x-m*f1.x, y:f2.y-m*f1.y};
  }
  return [f1,f2];
}
function fdPath(cv,ctx,verts){
  ctx.beginPath();
  verts.forEach((v,i)=>{i?ctx.lineTo(PX(cv,v.x),PY(cv,v.y)):ctx.moveTo(PX(cv,v.x),PY(cv,v.y));});
  ctx.closePath();
}
function drawLeft(){
  const ctx=cva.getContext("2d");
  ctx.clearRect(0,0,cva.width,cva.height);
  fdPath(cva,ctx,[pos(0,0),pos(1,0),pos(1,1),pos(0,1)]);
  ctx.fillStyle="rgba(77,163,216,0.10)"; ctx.fill();
  ctx.strokeStyle="rgba(255,255,255,0.35)"; ctx.lineWidth=1.2; ctx.stroke();
  // kernel bookkeeping: k*P mod d -> k
  const km={};
  for(let k=0;k<d;k++) km[((ker[0]*k)%d)+","+((ker[1]*k)%d)]=k;
  ctx.font="11px system-ui";
  for(let a=0;a<d;a++)for(let b=0;b<d;b++){
    const p=pos(a/d,b/d), key=a+","+b, k=km[key];
    const isO=(a===0&&b===0), isP=(a===ker[0]&&b===ker[1]);
    ctx.fillStyle= isO?MUT : (k!==undefined?RED:DOM);
    ctx.beginPath(); ctx.arc(PX(cva,p.x),PY(cva,p.y), isP?7:(k!==undefined?5.5:4.5), 0,7); ctx.fill();
    if(isP){ctx.strokeStyle="#fff"; ctx.lineWidth=1.6; ctx.stroke();}
    if(k!==undefined){
      ctx.fillStyle=isO?MUT:RED;
      ctx.fillText(isO?"O":(k===1?"P":k+"P"), PX(cva,p.x)+7, PY(cva,p.y)-6);
    }
  }
  // draggable tau corner
  const t=pos(0,1);
  ctx.fillStyle=SUP; ctx.beginPath(); ctx.arc(PX(cva,t.x),PY(cva,t.y),6,0,7); ctx.fill();
  ctx.strokeStyle="#000"; ctx.lineWidth=1; ctx.stroke();
  ctx.fillStyle=INK; ctx.fillText("τ", PX(cva,t.x)-14, PY(cva,t.y)-6);
}
function drawRight(){
  const ctx=cvb.getContext("2d");
  ctx.clearRect(0,0,cvb.width,cvb.height);
  // reduced fundamental cell of Lambda'
  const [f1,f2]=lprimeCell();
  fdPath(cvb,ctx,[{x:0,y:0},f1,{x:f1.x+f2.x,y:f1.y+f2.y},f2]);
  ctx.fillStyle="rgba(105,179,130,0.18)"; ctx.fill();
  ctx.strokeStyle=GRN; ctx.lineWidth=1.8; ctx.stroke();
  // Lambda (gray) + the cosets C adds (gold)
  for(let m=-2;m<=3;m++)for(let n=-2;n<=3;n++)for(let k=0;k<d;k++){
    const p=pos(m+k*ker[0]/d, n+k*ker[1]/d);
    if(p.x<CX-HS||p.x>CX+HS||p.y<CY-HS||p.y>CY+HS) continue;
    ctx.fillStyle=k?SUP:MUT;
    ctx.beginPath(); ctx.arc(PX(cvb,p.x),PY(cvb,p.y),k?4:4.5,0,7); ctx.fill();
  }
}
function render(){
  fitView(); drawLeft(); drawRight();
  out.innerHTML=`τ = ${tau.x.toFixed(2)} + ${tau.y.toFixed(2)}i`
    +` &nbsp;·&nbsp; <b style="color:${RED}">P</b> = (${ker[0]}, ${ker[1]})/${d}`
    +` &nbsp;·&nbsp; deg φ = #⟨P⟩ = ${d}, &nbsp; [Λ′ : Λ] = ${d}`;
}
function evXY(e){
  const b=cva.getBoundingClientRect();
  return {x:(CX-HS)+(e.clientX-b.left)/b.width*2*HS,
          y:(CY+HS)-(e.clientY-b.top)/b.height*2*HS};
}
cva.addEventListener("pointerdown",e=>{
  const q=evXY(e), t=pos(0,1);
  const dpx=PX(cva,t.x)-PX(cva,q.x), dpy=PY(cva,t.y)-PY(cva,q.y);
  cva.setPointerCapture(e.pointerId);
  if(dpx*dpx+dpy*dpy<14*14){ dragTau=true; }
  else{
    let best=null,bd=18*18;
    for(let a=0;a<d;a++)for(let b2=0;b2<d;b2++){
      if(a===0&&b2===0) continue;
      const p=pos(a/d,b2/d);
      const dx=PX(cva,p.x)-PX(cva,q.x), dy=PY(cva,p.y)-PY(cva,q.y);
      if(dx*dx+dy*dy<bd){bd=dx*dx+dy*dy; best=[a,b2];}
    }
    if(best) ker=best;
  }
  render(); e.preventDefault();
});
cva.addEventListener("pointermove",e=>{
  if(!dragTau) return;
  const q=evXY(e);
  tau.x=Math.max(-0.5,Math.min(0.5,q.x)); tau.y=Math.max(0.4,Math.min(2.0,q.y));
  render();
});
window.addEventListener("pointerup",()=>{dragTau=false;});
document.querySelectorAll(".dbtn").forEach(b=>b.addEventListener("click",()=>{
  d=+b.dataset.d; ker=[1,1%d===0?0:1];
  if(ker[0]===0&&ker[1]===0) ker=[1,0];
  document.querySelectorAll(".dbtn").forEach(z=>z.classList.toggle("on",z===b));
  render();
}));
render();
</script>
"""


def torus_folding_html() -> str:
    """S4 applet: an l-isogeny folds the torus l-to-1.

    Replaces the tau/s/t sliders + degree select + direction radio + static
    matplotlib pair. Left: the fundamental domain sliced into l slabs, the
    kernel (red) and the whole fibre of P; drag P (or tau), pick l and the
    kernel direction with buttons. Right: the folded cell of Lambda' with the
    single image phi(P), live."""
    return _HEAD + r"""
<div class="panel">
  <div class="modebar">
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin-right:4px;">ℓ =</span>
    <button class="seg lbtn" data-l="2">2</button>
    <button class="seg lbtn on" data-l="3">3</button>
    <button class="seg lbtn" data-l="4">4</button>
    <button class="seg lbtn" data-l="5">5</button>
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin:0 4px 0 14px;">kernel =</span>
    <button class="seg kbtn on" data-k="v">⟨1/ℓ⟩ vertical</button>
    <button class="seg kbtn" data-k="h">⟨τ/ℓ⟩ horizontal</button>
  </div>
  <div class="stage" style="grid-template-columns: 1fr 1fr;">
    <div class="cell">
      <div class="cap">E = ℂ/Λ, sliced into ℓ slabs — drag P (and τ)</div>
      <canvas id="fda" width="330" height="330"></canvas>
    </div>
    <div class="cell">
      <div class="cap">E′ = ℂ/Λ′ — the slabs folded onto one cell</div>
      <canvas id="fdb" width="330" height="330"></canvas>
    </div>
  </div>
  <div id="fdout" style="margin-top:8px;font-size:.93rem;"></div>
  <div class="hint" style="margin-top:4px;">white = the fibre of P (its ℓ kernel-translates, one per slab) · red = the kernel · all ℓ fibre points land on the single image φ(P)</div>
</div>
<script>
"use strict";
const INK="#d7d9dc", MUT="#9aa4ad", DOM="#4da3d8", SUP="#e0b64f", RED="#ef6f6f", GRN="#69b382";
const cva=document.getElementById("fda"), cvb=document.getElementById("fdb");
const out=document.getElementById("fdout");
let tau={x:0.25,y:1.0}, ell=3, vert=true, P={s:0.30,t:0.55};
let drag=null;                                  // 'tau' | 'P' | null

const pos=(cu,cv)=>({x:cu+cv*tau.x, y:cv*tau.y});
let CX,CY,HS;
function fitView(){
  const cs=[pos(0,0),pos(1,0),pos(1,1),pos(0,1)];
  const xs=cs.map(c=>c.x), ys=cs.map(c=>c.y);
  CX=(Math.min(...xs)+Math.max(...xs))/2; CY=(Math.min(...ys)+Math.max(...ys))/2;
  HS=Math.max(Math.max(...xs)-Math.min(...xs), Math.max(...ys)-Math.min(...ys))/2+0.25;
}
const PX=(cv,x)=>(x-(CX-HS))/(2*HS)*cv.width;
const PY=(cv,y)=>((CY+HS)-y)/(2*HS)*cv.height;
function poly(cv,ctx,verts){
  ctx.beginPath();
  verts.forEach((v,i)=>{i?ctx.lineTo(PX(cv,v.x),PY(cv,v.y)):ctx.moveTo(PX(cv,v.x),PY(cv,v.y));});
  ctx.closePath();
}
const slabCol=k=>`hsla(${215-135*k/Math.max(1,ell-1)},55%,50%,0.32)`;
function dot(cv,ctx,x,y,r,fill,ring){
  ctx.fillStyle=fill; ctx.beginPath(); ctx.arc(PX(cv,x),PY(cv,y),r,0,7); ctx.fill();
  if(ring){ctx.strokeStyle=ring; ctx.lineWidth=1.6; ctx.stroke();}
}
function drawLeft(){
  const ctx=cva.getContext("2d");
  ctx.clearRect(0,0,cva.width,cva.height);
  for(let k=0;k<ell;k++){
    const c = vert
      ? [pos(k/ell,0),pos((k+1)/ell,0),pos((k+1)/ell,1),pos(k/ell,1)]
      : [pos(0,k/ell),pos(1,k/ell),pos(1,(k+1)/ell),pos(0,(k+1)/ell)];
    poly(cva,ctx,c); ctx.fillStyle=slabCol(k); ctx.fill();
    ctx.strokeStyle="rgba(255,255,255,0.25)"; ctx.lineWidth=0.8; ctx.stroke();
  }
  // kernel points (red) and the fibre of P (white)
  for(let k=0;k<ell;k++){
    const K = vert?pos(k/ell,0):pos(0,k/ell);
    dot(cva,ctx,K.x,K.y,5.5,RED,"#fff");
  }
  ctx.font="12px system-ui";
  for(let k=0;k<ell;k++){
    const F = vert?pos((P.s+k/ell)%1,P.t):pos(P.s,(P.t+k/ell)%1);
    dot(cva,ctx,F.x,F.y,k?4.5:7,"#fff",k?null:DOM);
    if(!k){ctx.fillStyle=INK; ctx.fillText("P",PX(cva,F.x)+8,PY(cva,F.y)-6);}
  }
  const t=pos(0,1);
  dot(cva,ctx,t.x,t.y,6,SUP,"#000");
  ctx.fillStyle=INK; ctx.fillText("τ",PX(cva,t.x)-14,PY(cva,t.y)-6);
}
function drawRight(){
  const ctx=cvb.getContext("2d");
  ctx.clearRect(0,0,cvb.width,cvb.height);
  const u = vert?pos(1/ell,0):pos(1,0);
  const v = vert?pos(0,1):pos(0,1/ell);
  poly(cvb,ctx,[{x:0,y:0},u,{x:u.x+v.x,y:u.y+v.y},v]);
  ctx.fillStyle="rgba(105,179,130,0.18)"; ctx.fill();
  ctx.strokeStyle=GRN; ctx.lineWidth=1.8; ctx.stroke();
  dot(cvb,ctx,0,0,5.5,RED,"#fff");
  const si = vert?(P.s%(1/ell)):P.s, ti = vert?P.t:(P.t%(1/ell));
  const I = pos(si,ti);
  ctx.font="12px system-ui";
  dot(cvb,ctx,I.x,I.y,7,"#fff",SUP);
  ctx.fillStyle=INK; ctx.fillText("φ(P)",PX(cvb,I.x)+8,PY(cvb,I.y)-6);
}
function render(){
  fitView(); drawLeft(); drawRight();
  const lam = vert?`(1/${ell})ℤ + τℤ`:`ℤ + (τ/${ell})ℤ`;
  out.innerHTML=`τ = ${tau.x.toFixed(2)} + ${tau.y.toFixed(2)}i`
    +` &nbsp;·&nbsp; P: (s,t) = (${P.s.toFixed(2)}, ${P.t.toFixed(2)})`
    +` &nbsp;·&nbsp; Λ′ = ${lam} &nbsp;·&nbsp; φ is ${ell}-to-1`;
}
function evXY(e){
  const b=cva.getBoundingClientRect();
  return {x:(CX-HS)+(e.clientX-b.left)/b.width*2*HS,
          y:(CY+HS)-(e.clientY-b.top)/b.height*2*HS};
}
cva.addEventListener("pointerdown",e=>{
  const q=evXY(e), t=pos(0,1);
  cva.setPointerCapture(e.pointerId);
  const dx=PX(cva,t.x)-PX(cva,q.x), dy=PY(cva,t.y)-PY(cva,q.y);
  drag = (dx*dx+dy*dy<14*14) ? "tau" : "P";
  if(drag==="P"){const tt=q.y/tau.y, ss=q.x-tt*tau.x;
    P={s:((ss%1)+1)%1, t:((tt%1)+1)%1};}
  render(); e.preventDefault();
});
cva.addEventListener("pointermove",e=>{
  if(!drag) return;
  const q=evXY(e);
  if(drag==="tau"){tau.x=Math.max(-0.5,Math.min(0.5,q.x)); tau.y=Math.max(0.4,Math.min(2.0,q.y));}
  else{const tt=q.y/tau.y, ss=q.x-tt*tau.x;
    P={s:((ss%1)+1)%1, t:((tt%1)+1)%1};}
  render();
});
window.addEventListener("pointerup",()=>{drag=null;});
document.querySelectorAll(".lbtn").forEach(b=>b.addEventListener("click",()=>{
  ell=+b.dataset.l;
  document.querySelectorAll(".lbtn").forEach(z=>z.classList.toggle("on",z===b));
  render();
}));
document.querySelectorAll(".kbtn").forEach(b=>b.addEventListener("click",()=>{
  vert=(b.dataset.k==="v");
  document.querySelectorAll(".kbtn").forEach(z=>z.classList.toggle("on",z===b));
  render();
}));
render();
</script>
"""


def hasse_count_html() -> str:
    """S2 applet: count points, read off the trace, live in the Hasse interval.

    Replaces the p selectbox + f/g number inputs + static matplotlib Hasse
    line. Pick p with buttons, step f and g; the point count runs in JS
    (Euler-criterion character sum) and the trace lands on the interval.
    UPGRADE over the original: the gray dots are clickable -- choosing an
    admissible trace a searches for a curve realizing it (Deuring/Waterhouse:
    every |a| <= 2 sqrt p occurs), so the applet runs in both directions."""
    return _HEAD + r"""
<div class="panel">
  <div class="modebar" style="flex-wrap:wrap;">
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin-right:4px;">p =</span>
    <button class="seg pbtn" data-p="5">5</button>
    <button class="seg pbtn" data-p="11">11</button>
    <button class="seg pbtn on" data-p="23">23</button>
    <button class="seg pbtn" data-p="47">47</button>
    <button class="seg pbtn" data-p="101">101</button>
    <button class="seg pbtn" data-p="199">199</button>
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin:0 2px 0 16px;">f =</span>
    <button class="seg" id="fdn">−</button><span id="fval" style="align-self:center;min-width:26px;text-align:center;">1</span><button class="seg" id="fup">+</button>
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin:0 2px 0 10px;">g =</span>
    <button class="seg" id="gdn">−</button><span id="gval" style="align-self:center;min-width:26px;text-align:center;">1</span><button class="seg" id="gup">+</button>
  </div>
  <div class="stage" style="grid-template-columns: 1fr;">
    <div class="cell">
      <div class="cap">the Hasse interval |a| ≤ 2√p — every dot is an isogeny class</div>
      <canvas id="hc" width="640" height="130"></canvas>
      <div id="hcout" style="margin-top:8px;font-size:.93rem;line-height:1.7;"></div>
      <div class="hint" style="margin-top:4px;">step f, g to change the curve — or <b>click a gray dot</b> and a curve with that trace is found for you (every admissible a occurs: Deuring)</div>
    </div>
  </div>
</div>
<script>
"use strict";
const INK="#d7d9dc", MUT="#9aa4ad", DOM="#4da3d8", SUP="#e0b64f", RED="#ef6f6f";
const cv=document.getElementById("hc"), ctx=cv.getContext("2d");
const out=document.getElementById("hcout");
let p=23, f=1, g=1;

function powmod(b,e,m){let r=1n;b=BigInt(b)%BigInt(m);let E=BigInt(e);const M=BigInt(m);
  while(E>0n){if(E&1n)r=r*b%M; b=b*b%M; E>>=1n;} return Number(r);}
const chi=(v,pp)=>{v=((v%pp)+pp)%pp; if(v===0)return 0;
  return powmod(v,(pp-1)/2,pp)===1?1:-1;};
function trace(ff,gg,pp){                    // a = -sum chi(x^3+fx+g); #E = p+1-a
  let s=0;
  for(let x=0;x<pp;x++) s+=chi(x*x*x+ff*x+gg,pp);
  return -s;
}
const disc=(ff,gg,pp)=>((-16*(4*ff*ff*ff+27*gg*gg))%pp+pp)%pp;
function findCurve(a0){                      // search a curve with trace a0
  for(let tries=0;tries<6000;tries++){
    const ff=Math.floor(Math.random()*p), gg=Math.floor(Math.random()*p);
    if(disc(ff,gg,p)===0) continue;
    if(trace(ff,gg,p)===a0){f=ff; g=gg; return true;}
  }
  return false;
}
const amax=()=>Math.floor(2*Math.sqrt(p));
const AX=a=>cv.width/2 + a/(2*Math.sqrt(p)+1.5)*(cv.width/2-14);
const AY=cv.height/2+8;

function render(){
  ctx.clearRect(0,0,cv.width,cv.height);
  const bound=2*Math.sqrt(p), sing=disc(f,g,p)===0;
  const a=sing?null:trace(f,g,p);
  // Hasse bounds + axis
  ctx.strokeStyle="rgba(255,255,255,0.2)"; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(AX(-bound-1),AY); ctx.lineTo(AX(bound+1),AY); ctx.stroke();
  ctx.strokeStyle=DOM; ctx.setLineDash([5,4]); ctx.lineWidth=1.3;
  for(const s of [-bound,bound]){ctx.beginPath();ctx.moveTo(AX(s),AY-34);ctx.lineTo(AX(s),AY+26);ctx.stroke();}
  ctx.setLineDash([]);
  ctx.fillStyle=DOM; ctx.font="11px system-ui"; ctx.textAlign="center";
  ctx.fillText(`−2√p ≈ ${(-bound).toFixed(1)}`, AX(-bound), AY+42);
  ctx.fillText(`+2√p ≈ ${bound.toFixed(1)}`, AX(bound), AY+42);
  // admissible traces
  for(let v=-amax();v<=amax();v++){
    ctx.fillStyle=(a!==null&&v===a)?RED:"rgba(255,255,255,0.4)";
    ctx.beginPath(); ctx.arc(AX(v),AY,(a!==null&&v===a)?7:4,0,7); ctx.fill();
  }
  if(a!==null){ctx.fillStyle=RED; ctx.font="bold 13px system-ui";
    ctx.fillText(`a = ${a}`, AX(a), AY-18);}
  ctx.textAlign="left";
  // readout
  if(sing){
    out.innerHTML=`<span style="color:${SUP}">y² = x³ + ${f}x + ${g} is singular mod ${p}</span> — step f or g (or click a dot).`;
  } else {
    const n=p+1-a, d=a*a-4*p;
    out.innerHTML=`E: y² = x³ + ${f}x + ${g} over 𝔽<sub>${p}</sub>`
      +` &nbsp;·&nbsp; #E = ${n} &nbsp;·&nbsp; <b style="color:${RED}">a = p + 1 − #E = ${a}</b>`
      +` &nbsp;·&nbsp; d = a² − 4p = ${d}`
      +` &nbsp;·&nbsp; ${a%p===0?`<b style="color:${SUP}">supersingular</b> (a ≡ 0 mod p)`:"ordinary"}`;
  }
}
cv.addEventListener("pointerdown",e=>{
  const b=cv.getBoundingClientRect();
  const mx=(e.clientX-b.left)*cv.width/b.width, my=(e.clientY-b.top)*cv.height/b.height;
  if(Math.abs(my-AY)>22) return;
  let best=null,bd=14*14;
  for(let v=-amax();v<=amax();v++){
    const dx=AX(v)-mx, dy=AY-my;
    if(dx*dx+dy*dy<bd){bd=dx*dx+dy*dy; best=v;}
  }
  if(best!==null){ if(findCurve(best)) render(); }
});
const clamp=v=>((v%p)+p)%p;
document.getElementById("fup").onclick=()=>{f=clamp(f+1);sync();};
document.getElementById("fdn").onclick=()=>{f=clamp(f-1);sync();};
document.getElementById("gup").onclick=()=>{g=clamp(g+1);sync();};
document.getElementById("gdn").onclick=()=>{g=clamp(g-1);sync();};
function sync(){document.getElementById("fval").textContent=f;
  document.getElementById("gval").textContent=g; render();}
document.querySelectorAll(".pbtn").forEach(b=>b.addEventListener("click",()=>{
  p=+b.dataset.p; f=Math.min(f,p-1); g=Math.min(g,p-1);
  document.querySelectorAll(".pbtn").forEach(z=>z.classList.toggle("on",z===b));
  sync();
}));
sync();
</script>
"""


def frobenius_lift_html() -> str:
    """S2 applet: Frobenius lifted -- multiplication by alpha on Z + alpha Z.

    Replaces the p selectbox + trace select_slider + static matplotlib lattice
    figure. Left: the (a, p) picker -- for each small prime a row of admissible
    traces under the Hasse parabola; click one. Right: the lattice
    Lambda = Z + alpha Z with the fundamental cell (blue) and its image under
    multiplication by alpha (gold, dashed) -- rotation by arg alpha, scaling by
    |alpha| = sqrt p; both image generators land on lattice points."""
    return _HEAD + r"""
<div class="panel">
  <div class="stage" style="grid-template-columns: 300px 1fr;">
    <div class="cell">
      <div class="cap">pick (a, p) with a² &lt; 4p — click a dot</div>
      <canvas id="fla" width="300" height="330"></canvas>
    </div>
    <div class="cell">
      <div class="cap">Λ = ℤ + αℤ and multiplication by α</div>
      <canvas id="flb" width="340" height="330"></canvas>
    </div>
  </div>
  <div id="flout" style="margin-top:8px;font-size:.93rem;line-height:1.7;"></div>
  <div class="hint" style="margin-top:4px;">blue cell: spanned by 1, α · gold dashed cell: its image under ×α — rotated by arg α, scaled by |α| = √p, with α·1 and α·α on lattice points: αΛ ⊆ Λ</div>
</div>
<script>
"use strict";
const INK="#d7d9dc", MUT="#9aa4ad", DOM="#4da3d8", SUP="#e0b64f", RED="#ef6f6f";
const PRIMES=[5,7,11,13,17,19,23];
const cva=document.getElementById("fla"), cvb=document.getElementById("flb");
const out=document.getElementById("flout");
let p=7, a=1;

// ---- left: the (a, p) picker --------------------------------------------
const AMAXALL=Math.floor(2*Math.sqrt(23));
const LX=v=>56 + (v+AMAXALL)/(2*AMAXALL)*(cva.width-70);
const LY=i=>36+i*(cva.height-64)/(PRIMES.length-1);
function drawPicker(){
  const ctx=cva.getContext("2d");
  ctx.clearRect(0,0,cva.width,cva.height);
  ctx.font="11px system-ui";
  PRIMES.forEach((pp,i)=>{
    const am=Math.floor(2*Math.sqrt(pp));
    ctx.fillStyle=MUT; ctx.textAlign="left";
    ctx.fillText("p = "+pp, 4, LY(i)+4);
    for(let v=-am;v<=am;v++){
      if(v===0) continue;
      const on=(pp===p&&v===a);
      ctx.fillStyle=on?RED:"rgba(255,255,255,0.4)";
      ctx.beginPath(); ctx.arc(LX(v),LY(i),on?6:3.5,0,7); ctx.fill();
    }
  });
  // Hasse parabola a^2 = 4p (through the row endpoints)
  ctx.strokeStyle="rgba(77,163,216,0.5)"; ctx.setLineDash([4,4]); ctx.lineWidth=1.2;
  for(const s of [-1,1]){
    ctx.beginPath();
    PRIMES.forEach((pp,i)=>{
      const x=LX(s*2*Math.sqrt(pp)), y=LY(i);
      i?ctx.lineTo(x,y):ctx.moveTo(x,y);
    });
    ctx.stroke();
  }
  ctx.setLineDash([]); ctx.textAlign="left";
}
cva.addEventListener("pointerdown",e=>{
  const b=cva.getBoundingClientRect();
  const mx=(e.clientX-b.left)*cva.width/b.width, my=(e.clientY-b.top)*cva.height/b.height;
  let best=null,bd=13*13;
  PRIMES.forEach((pp,i)=>{
    const am=Math.floor(2*Math.sqrt(pp));
    for(let v=-am;v<=am;v++){
      if(v===0)continue;
      const dx=LX(v)-mx, dy=LY(i)-my;
      if(dx*dx+dy*dy<bd){bd=dx*dx+dy*dy; best=[pp,v];}
    }
  });
  if(best){p=best[0]; a=best[1]; render();}
});

// ---- right: the lattice and the two cells --------------------------------
function drawLattice(){
  const ctx=cvb.getContext("2d");
  ctx.clearRect(0,0,cvb.width,cvb.height);
  const ax=a/2, ay=Math.sqrt(4*p-a*a)/2;              // alpha
  const i1={x:ax,y:ay};                                // alpha*1
  const i2={x:a*ax-p, y:a*ay};                         // alpha^2 = a*alpha - p
  const keys=[{x:0,y:0},{x:1,y:0},i1,{x:1+ax,y:ay},i2,{x:i1.x+i2.x,y:i1.y+i2.y}];
  const pad=1.0;
  let x0=Math.min(...keys.map(k=>k.x))-pad, x1=Math.max(...keys.map(k=>k.x))+pad,
      y0=Math.min(...keys.map(k=>k.y))-pad, y1=Math.max(...keys.map(k=>k.y))+pad;
  const need=(x1-x0)/(y1-y0), have=cvb.width/cvb.height;
  if(need<have){const c=(x0+x1)/2,w=(y1-y0)*have; x0=c-w/2; x1=c+w/2;}
  else{const c=(y0+y1)/2,h=(x1-x0)/have; y0=c-h/2; y1=c+h/2;}
  const PX=x=>(x-x0)/(x1-x0)*cvb.width, PY=y=>(y1-y)/(y1-y0)*cvb.height;
  // lattice m + n*alpha
  ctx.fillStyle="rgba(255,255,255,0.3)";
  for(let n=Math.floor(y0/ay)-1;n<=Math.ceil(y1/ay)+1;n++)
    for(let m=Math.floor(x0-n*ax)-1;m<=Math.ceil(x1-n*ax)+1;m++){
      const x=m+n*ax, y=n*ay;
      if(x>x0&&x<x1&&y>y0&&y<y1){ctx.beginPath();ctx.arc(PX(x),PY(y),2.2,0,7);ctx.fill();}
    }
  // axes
  ctx.strokeStyle="rgba(255,255,255,0.12)"; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(PX(x0),PY(0)); ctx.lineTo(PX(x1),PY(0));
  ctx.moveTo(PX(0),PY(y0)); ctx.lineTo(PX(0),PY(y1)); ctx.stroke();
  const poly=(vs,fill,stroke,dash)=>{
    ctx.beginPath();
    vs.forEach((v,i)=>{i?ctx.lineTo(PX(v.x),PY(v.y)):ctx.moveTo(PX(v.x),PY(v.y));});
    ctx.closePath();
    if(dash)ctx.setLineDash([5,4]);
    ctx.fillStyle=fill; ctx.fill(); ctx.strokeStyle=stroke; ctx.lineWidth=1.6; ctx.stroke();
    ctx.setLineDash([]);
  };
  poly([{x:0,y:0},{x:1,y:0},{x:1+ax,y:ay},i1],"rgba(77,163,216,0.16)",DOM,false);
  poly([{x:0,y:0},i1,{x:i1.x+i2.x,y:i1.y+i2.y},i2],"rgba(224,182,79,0.12)",SUP,true);
  const arrow=(v,col,dash)=>{
    ctx.strokeStyle=col; ctx.fillStyle=col; ctx.lineWidth=1.8;
    if(dash)ctx.setLineDash([5,4]);
    ctx.beginPath(); ctx.moveTo(PX(0),PY(0)); ctx.lineTo(PX(v.x),PY(v.y)); ctx.stroke();
    ctx.setLineDash([]);
    const an=Math.atan2(PY(v.y)-PY(0),PX(v.x)-PX(0));
    ctx.beginPath(); ctx.moveTo(PX(v.x),PY(v.y));
    ctx.lineTo(PX(v.x)-9*Math.cos(an-0.42),PY(v.y)-9*Math.sin(an-0.42));
    ctx.lineTo(PX(v.x)-9*Math.cos(an+0.42),PY(v.y)-9*Math.sin(an+0.42));
    ctx.closePath(); ctx.fill();
  };
  ctx.font="12px system-ui";
  arrow({x:1,y:0},"#fff",false); ctx.fillStyle="#fff"; ctx.fillText("1",PX(1)+5,PY(0)+14);
  arrow(i1,RED,false); ctx.fillStyle=RED; ctx.fillText("α = α·1",PX(i1.x)+7,PY(i1.y)-6);
  arrow(i2,SUP,true); ctx.fillStyle=SUP; ctx.fillText("α·α = aα − p",PX(i2.x)+7,PY(i2.y)-6);
}
function render(){
  drawPicker(); drawLattice();
  const ay=Math.sqrt(4*p-a*a)/2, d=a*a-4*p;
  out.innerHTML=`p = ${p}, a = ${a} &nbsp;·&nbsp; d = a² − 4p = ${d}`
    +` &nbsp;·&nbsp; α = ${(a/2).toFixed(2)} + ${ay.toFixed(2)}i, |α| = √${p} ≈ ${Math.sqrt(p).toFixed(2)}`
    +` &nbsp;·&nbsp; on the basis (1, α): [α] = ( 0&nbsp;&nbsp;−${p}&nbsp;; 1&nbsp;&nbsp;${a} ) — integer entries ⇒ αΛ ⊆ Λ`;
}
render();
</script>
"""


def cubic_family_html() -> str:
    """S1 applet: the family y^2 = x^3 + fx + g over the (f, g)-plane.

    Replaces the f/g sliders + matplotlib curve. Left: the (f, g)-plane with
    the discriminant locus 4f^3 + 27g^2 = 0 drawn (the cuspidal curve
    separating one-component from two-component curves); drag the point --
    it SNAPS onto the locus nearby, so the singular curves (node, and the
    cusp at f = g = 0) are actually reachable. Right: the real curve, live,
    with the roots of the cubic marked."""
    return _HEAD + r"""
<div class="panel">
  <div class="stage" style="grid-template-columns: 300px 1fr;">
    <div class="cell">
      <div class="cap">the (f, g)-plane — drag; Δ = 0 locus drawn</div>
      <canvas id="cfa" width="300" height="330"></canvas>
    </div>
    <div class="cell">
      <div class="cap">y² = x³ + fx + g over ℝ</div>
      <canvas id="cfb" width="340" height="330"></canvas>
    </div>
  </div>
  <div id="cfout" style="margin-top:8px;font-size:.93rem;"></div>
  <div class="hint" style="margin-top:4px;">inside the blue cusp region: Δ &gt; 0, two components · outside: Δ &lt; 0, one component · the marker snaps onto Δ = 0, where the curve degenerates (node; cusp at f = g = 0)</div>
</div>
<script>
"use strict";
const INK="#d7d9dc", MUT="#9aa4ad", DOM="#4da3d8", SUP="#e0b64f", RED="#ef6f6f";
const cva=document.getElementById("cfa"), cvb=document.getElementById("cfb");
const out=document.getElementById("cfout");
const F0=-5,F1=5,G0=-5,G1=5;
let f=-1.0, g=1.0, dragging=false;

const PXa=x=>(x-F0)/(F1-F0)*cva.width, PYa=y=>(G1-y)/(G1-G0)*cva.height;
const gLoc=ff=>Math.sqrt(Math.max(0,-4*ff*ff*ff/27));   // |g| on the locus
function drawPlane(){
  const ctx=cva.getContext("2d");
  ctx.clearRect(0,0,cva.width,cva.height);
  // region Delta > 0 (between the branches), tinted
  ctx.fillStyle="rgba(77,163,216,0.12)";
  ctx.beginPath();
  for(let i=0;i<=120;i++){const ff=F0+i*(0-F0)/120; const px=PXa(ff),py=PYa(gLoc(ff));
    i?ctx.lineTo(px,py):ctx.moveTo(px,py);}
  for(let i=120;i>=0;i--){const ff=F0+i*(0-F0)/120;
    ctx.lineTo(PXa(ff),PYa(-gLoc(ff)));}
  ctx.closePath(); ctx.fill();
  // axes
  ctx.strokeStyle="rgba(255,255,255,0.12)"; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(PXa(F0),PYa(0)); ctx.lineTo(PXa(F1),PYa(0));
  ctx.moveTo(PXa(0),PYa(G0)); ctx.lineTo(PXa(0),PYa(G1)); ctx.stroke();
  // the discriminant locus
  ctx.strokeStyle=SUP; ctx.lineWidth=1.8;
  for(const s of [1,-1]){
    ctx.beginPath();
    for(let i=0;i<=160;i++){const ff=F0+i*(0-F0)/160;
      const px=PXa(ff),py=PYa(s*gLoc(ff));
      i?ctx.lineTo(px,py):ctx.moveTo(px,py);}
    ctx.stroke();
  }
  ctx.font="11px system-ui"; ctx.fillStyle=MUT;
  ctx.fillText("f", PXa(F1)-12, PYa(0)-6);
  ctx.fillText("g", PXa(0)+6, PYa(G1)+14);
  ctx.fillStyle=DOM; ctx.fillText("Δ > 0", PXa(-4.2), PYa(0.4));
  ctx.fillStyle=MUT; ctx.fillText("Δ < 0", PXa(2.2), PYa(3.6));
  // marker
  const disc=-16*(4*f*f*f+27*g*g), sing=Math.abs(disc)<1e-9;
  ctx.fillStyle=sing?SUP:"#fff";
  ctx.beginPath(); ctx.arc(PXa(f),PYa(g),6,0,7); ctx.fill();
  ctx.strokeStyle=sing?RED:DOM; ctx.lineWidth=2; ctx.stroke();
}
function drawCurve(){
  const ctx=cvb.getContext("2d"), X0=-3.3,X1=3.3,Y0=-4.5,Y1=4.5;
  const PX=x=>(x-X0)/(X1-X0)*cvb.width, PY=y=>(Y1-y)/(Y1-Y0)*cvb.height;
  ctx.clearRect(0,0,cvb.width,cvb.height);
  ctx.strokeStyle="rgba(255,255,255,0.12)"; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(PX(X0),PY(0)); ctx.lineTo(PX(X1),PY(0));
  ctx.moveTo(PX(0),PY(Y0)); ctx.lineTo(PX(0),PY(Y1)); ctx.stroke();
  const NS=900;
  ctx.strokeStyle=DOM; ctx.lineWidth=2;
  for(const s of [1,-1]){
    ctx.beginPath(); let pen=false;
    for(let i=0;i<=NS;i++){
      const x=X0+(X1-X0)*i/NS, y2=x*x*x+f*x+g;
      if(y2>=0){const py=PY(s*Math.sqrt(y2));
        pen?ctx.lineTo(PX(x),py):ctx.moveTo(PX(x),py); pen=true;}
      else pen=false;
    }
    ctx.stroke();
  }
  // roots of the cubic on the x-axis
  ctx.fillStyle=SUP;
  let prev=X0*X0*X0+f*X0+g;
  for(let i=1;i<=NS;i++){
    const x=X0+(X1-X0)*i/NS, v=x*x*x+f*x+g;
    if(prev*v<=0&&(prev!==0||v!==0)){
      let lo=x-(X1-X0)/NS, hi=x;
      for(let k=0;k<30;k++){const m=(lo+hi)/2;
        ((lo*lo*lo+f*lo+g)*(m*m*m+f*m+g)<=0)?hi=m:lo=m;}
      ctx.beginPath(); ctx.arc(PX((lo+hi)/2),PY(0),3.5,0,7); ctx.fill();
    }
    prev=v;
  }
}
function render(){
  drawPlane(); drawCurve();
  const disc=-16*(4*f*f*f+27*g*g);
  let sh;
  if(Math.abs(disc)<1e-9) sh=`<b style="color:${RED}">singular</b> — ${(Math.abs(f)<1e-9&&Math.abs(g)<1e-9)?"a cusp":"a node"} (Δ = 0)`;
  else sh=disc>0?"two components (Δ > 0)":"one component (Δ < 0)";
  out.innerHTML=`y² = x³ ${f<0?"−":"+"} ${Math.abs(f).toFixed(1)}x ${g<0?"−":"+"} ${Math.abs(g).toFixed(1)}`
    +` &nbsp;·&nbsp; Δ = −16(4f³ + 27g²) = ${disc.toFixed(1)} &nbsp;·&nbsp; ${sh}`;
}
function set(e){
  const b=cva.getBoundingClientRect();
  let ff=F0+(e.clientX-b.left)/b.width*(F1-F0),
      gg=G1-(e.clientY-b.top)/b.height*(G1-G0);
  ff=Math.max(F0,Math.min(F1,ff)); gg=Math.max(G0,Math.min(G1,gg));
  // snap to the discriminant locus when close (in pixels)
  if(ff<=0){
    const gl=gLoc(ff), s=gg>=0?1:-1;
    if(Math.abs(PYa(gg)-PYa(s*gl))<7 && Math.abs(gg)<G1) gg=s*gl;
    if(Math.abs(ff)<0.15&&Math.abs(gg)<0.15&&Math.abs(PXa(ff)-PXa(0))<6) {ff=0;gg=0;}
  }
  f=ff; g=gg;
}
cva.addEventListener("pointerdown",e=>{dragging=true; cva.setPointerCapture(e.pointerId); set(e); render(); e.preventDefault();});
cva.addEventListener("pointermove",e=>{if(dragging){set(e); render();}});
window.addEventListener("pointerup",()=>{dragging=false;});
render();
</script>
"""


def velu_html() -> str:
    """S4 applet: Velu's formulas over F_p, numerically.

    Replaces the p selectbox + f/g number inputs + degree select_slider +
    x(P) select_slider + y(P) radio + static matplotlib pair. Controls are
    buttons/steppers; the source point P is chosen by CLICKING a point of
    E(F_p) on the left grid (kernel red, P green), and phi(P) lights up on
    the codomain grid. When the current curve has no order-l kernel the
    applet searches out one that does. Math is a JS port of the page's
    _v_kernel/_v_run (Velu sums, 3-division polynomial for l=3)."""
    return _HEAD + r"""
<div class="panel">
  <div class="modebar" style="flex-wrap:wrap;">
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin-right:4px;">p =</span>
    <button class="seg vpbtn" data-p="11">11</button>
    <button class="seg vpbtn on" data-p="23">23</button>
    <button class="seg vpbtn" data-p="31">31</button>
    <button class="seg vpbtn" data-p="47">47</button>
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin:0 4px 0 14px;">ℓ =</span>
    <button class="seg vlbtn on" data-l="2">2</button>
    <button class="seg vlbtn" data-l="3">3</button>
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin:0 2px 0 14px;">f =</span>
    <button class="seg" id="vfdn">−</button><span id="vfval" style="align-self:center;min-width:26px;text-align:center;">1</span><button class="seg" id="vfup">+</button>
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin:0 2px 0 10px;">g =</span>
    <button class="seg" id="vgdn">−</button><span id="vgval" style="align-self:center;min-width:26px;text-align:center;">1</span><button class="seg" id="vgup">+</button>
  </div>
  <div class="stage" style="grid-template-columns: 1fr 1fr;">
    <div class="cell">
      <div class="cap">E(𝔽ₚ) — kernel red; <b>click a point P</b></div>
      <canvas id="vla" width="330" height="330"></canvas>
    </div>
    <div class="cell">
      <div class="cap">(E/C)(𝔽ₚ) — Vélu's codomain, with φ(P)</div>
      <canvas id="vlb" width="330" height="330"></canvas>
    </div>
  </div>
  <div id="vlout" style="margin-top:8px;font-size:.93rem;line-height:1.7;"></div>
  <div class="hint" style="margin-top:4px;">both grids use symmetric representatives (−p/2, p/2] · #E = #(E/C): isogenous curves have equal point counts · if the curve has no order-ℓ kernel over 𝔽ₚ, one is searched out for you</div>
</div>
<script>
"use strict";
const INK="#d7d9dc", MUT="#9aa4ad", DOM="#4da3d8", SUP="#e0b64f", RED="#ef6f6f", GRN="#69b382";
const cva=document.getElementById("vla"), cvb=document.getElementById("vlb");
const out=document.getElementById("vlout");
let p=23, ell=2, f=1, g=1, P=null;
let state=null;   // {S, tt, f2, g2, Epts, E2pts, kerSet}

const mod=(a,m)=>((a%m)+m)%m;
function powmod(b,e,m){let r=1n;b=BigInt(mod(b,m));let E=BigInt(e);const M=BigInt(m);
  while(E>0n){if(E&1n)r=r*b%M; b=b*b%M; E>>=1n;} return Number(r);}
const inv=(a,m)=>powmod(a,m-2,m);
const sym=x=>{const r=mod(x,p); return 2*r<p?r:r-p;};
const rhs=(x,ff,gg)=>mod(x*x*x+ff*x+gg,p);
function curvePts(ff,gg){
  const sq={};                                  // sqrt table
  for(let y=0;y<p;y++){const v=mod(y*y,p); (sq[v]=sq[v]||[]).push(y);}
  const pts=[];
  for(let x=0;x<p;x++){const r=rhs(x,ff,gg);
    for(const y of (sq[r]||[])) pts.push([x,y]);}
  return pts;
}
function kernel(ff,gg){
  if(ell===2){
    for(let x=0;x<p;x++) if(rhs(x,ff,gg)===0) return {S:[[x,0]], tt:new Set([0])};
    return null;
  }
  for(let x=0;x<p;x++){
    if(mod(3*x*x*x*x+6*ff*x*x+12*gg*x-ff*ff,p)!==0) continue;
    const r=rhs(x,ff,gg);
    for(let y=1;y<p;y++) if(mod(y*y,p)===r) return {S:[[x,y]], tt:new Set()};
  }
  return null;
}
function velu(ff,gg,S,tt){
  let v=0,w=0; const data=[];
  S.forEach((q,i)=>{
    const gx=mod(3*q[0]*q[0]+ff,p), uq=mod(4*q[1]*q[1],p);
    const vq=tt.has(i)?gx:mod(2*gx,p);
    v=mod(v+vq,p); w=mod(w+uq+q[0]*vq,p);
    data.push([q[0],q[1],uq,vq]);
  });
  const f2=mod(ff-5*v,p), g2=mod(gg-7*w,p);
  const phi=(x,y)=>{
    let X=mod(x,p), corr=0;
    for(const [xq,yq,uq,vq] of data){
      const iv=inv(mod(x-xq,p),p);
      X=mod(X+vq*iv+uq*iv*iv,p);
      corr=mod(corr+vq*iv*iv+2*uq*powmod(iv,3,p),p);
    }
    return [X, mod(y*mod(1-corr,p),p)];
  };
  return {f2,g2,phi};
}
function rebuild(searchIfMissing){
  const disc=mod(-16*(4*f*f*f+27*g*g),p);
  let ker=disc===0?null:kernel(f,g);
  if(!ker&&searchIfMissing){
    for(let t=0;t<4000&&!ker;t++){
      const ff=Math.floor(Math.random()*p), gg=Math.floor(Math.random()*p);
      if(mod(-16*(4*ff*ff*ff+27*gg*gg),p)===0) continue;
      const k=kernel(ff,gg);
      if(k){f=ff; g=gg; ker=k;}
    }
    syncFG();
  }
  if(!ker){state=null; P=null; return;}
  const {S,tt}=ker, {f2,g2,phi}=velu(f,g,S,tt);
  const kerSet=new Set();
  for(const [xq,yq] of S){
    if(ell===2) kerSet.add(xq+","+0);
    else {kerSet.add(xq+","+yq); kerSet.add(xq+","+mod(-yq,p));}
  }
  state={S,tt,f2,g2,phi,kerSet,Epts:curvePts(f,g),E2pts:curvePts(f2,g2)};
  if(P&&(state.kerSet.has(P[0]+","+P[1])||rhs(P[0],f,g)!==mod(P[1]*P[1],p))) P=null;
  if(!P) P=state.Epts.find(q=>!state.kerSet.has(q[0]+","+q[1]))||null;
}
const h=()=>Math.floor(p/2);
const GX=(cv,x)=>(sym(x)+h()+0.5)/(p)*cv.width;
const GY=(cv,y)=>(h()-sym(y)+0.5)/(p)*cv.height;
function drawGrid(cv,pts,ffs,extra){
  const ctx=cv.getContext("2d");
  ctx.clearRect(0,0,cv.width,cv.height);
  ctx.fillStyle="rgba(255,255,255,0.10)";
  for(let x=0;x<p;x++)for(let y=0;y<p;y++){
    ctx.beginPath(); ctx.arc(GX(cv,x),GY(cv,y),1.2,0,7); ctx.fill();
  }
  ctx.fillStyle=DOM;
  for(const [x,y] of pts){
    ctx.beginPath(); ctx.arc(GX(cv,x),GY(cv,y),Math.max(2.6,60/p),0,7); ctx.fill();
  }
  if(extra) extra(ctx);
}
function render(){
  if(!state){
    drawGrid(cva,[],f); drawGrid(cvb,[],0);
    out.innerHTML=`<span style="color:${SUP}">y² = x³ + ${f}x + ${g} has no order-${ell} kernel over 𝔽<sub>${p}</sub> (or is singular)</span> — step f/g, or switch ℓ/p to auto-search.`;
    return;
  }
  const st=state;
  drawGrid(cva,st.Epts,f,ctx=>{
    ctx.font="12px system-ui";
    for(const key of st.kerSet){
      const [x,y]=key.split(",").map(Number);
      ctx.fillStyle=RED; ctx.beginPath();
      ctx.arc(GX(cva,x),GY(cva,y),Math.max(4,80/p),0,7); ctx.fill();
      ctx.strokeStyle="#fff"; ctx.lineWidth=1; ctx.stroke();
    }
    if(P){
      ctx.fillStyle=GRN; ctx.beginPath();
      ctx.arc(GX(cva,P[0]),GY(cva,P[1]),Math.max(5,90/p),0,7); ctx.fill();
      ctx.strokeStyle="#fff"; ctx.lineWidth=1.4; ctx.stroke();
      ctx.fillStyle=GRN; ctx.fillText("P",GX(cva,P[0])+8,GY(cva,P[1])-6);
    }
  });
  const IM=P?st.phi(P[0],P[1]):null;
  drawGrid(cvb,st.E2pts,st.f2,ctx=>{
    if(IM){
      ctx.font="12px system-ui";
      ctx.fillStyle=GRN; ctx.beginPath();
      ctx.arc(GX(cvb,IM[0]),GY(cvb,IM[1]),Math.max(5,90/p),0,7); ctx.fill();
      ctx.strokeStyle="#fff"; ctx.lineWidth=1.4; ctx.stroke();
      ctx.fillStyle=GRN; ctx.fillText("φ(P)",GX(cvb,IM[0])+8,GY(cvb,IM[1])-6);
    }
  });
  const q0=st.S[0];
  out.innerHTML=`E: y² = x³ + ${f}x + ${g} &nbsp;·&nbsp; kernel gen (${sym(q0[0])}, ${sym(q0[1])})`
    +` &nbsp;·&nbsp; <b style="color:${SUP}">E/C: Y² = X³ + ${st.f2}X + ${st.g2}</b> (mod ${p})`
    +(P?` &nbsp;·&nbsp; <b style="color:${GRN}">P = (${sym(P[0])}, ${sym(P[1])}) ↦ φ(P) = (${sym(IM[0])}, ${sym(IM[1])})</b>`:"")
    +` &nbsp;·&nbsp; #E = #(E/C) = ${st.Epts.length+1}`;
}
cva.addEventListener("pointerdown",e=>{
  if(!state) return;
  const b=cva.getBoundingClientRect();
  const mx=(e.clientX-b.left)*cva.width/b.width, my=(e.clientY-b.top)*cva.height/b.height;
  let best=null,bd=15*15;
  for(const q of state.Epts){
    if(state.kerSet.has(q[0]+","+q[1])) continue;
    const dx=GX(cva,q[0])-mx, dy=GY(cva,q[1])-my;
    if(dx*dx+dy*dy<bd){bd=dx*dx+dy*dy; best=q;}
  }
  if(best){P=best; render();}
});
function syncFG(){document.getElementById("vfval").textContent=f;
  document.getElementById("vgval").textContent=g;}
const step=(df,dg)=>{f=mod(f+df,p); g=mod(g+dg,p); syncFG(); P=null; rebuild(false); render();};
document.getElementById("vfup").onclick=()=>step(1,0);
document.getElementById("vfdn").onclick=()=>step(-1,0);
document.getElementById("vgup").onclick=()=>step(0,1);
document.getElementById("vgdn").onclick=()=>step(0,-1);
document.querySelectorAll(".vpbtn").forEach(b=>b.addEventListener("click",()=>{
  p=+b.dataset.p; f=mod(f,p); g=mod(g,p); P=null;
  document.querySelectorAll(".vpbtn").forEach(z=>z.classList.toggle("on",z===b));
  syncFG(); rebuild(true); render();
}));
document.querySelectorAll(".vlbtn").forEach(b=>b.addEventListener("click",()=>{
  ell=+b.dataset.l; P=null;
  document.querySelectorAll(".vlbtn").forEach(z=>z.classList.toggle("on",z===b));
  rebuild(true); render();
}));
rebuild(true); render();
</script>
"""


def group_law_real_html() -> str:
    """S1 applet: the chord-tangent group law on a real curve, live.

    Replaces the f/g sliders + plotly click-samples + Clear/Show line/Compute
    sum buttons. Click near the curve to place P, then Q (both snap onto the
    curve and stay draggable); the whole construction -- chord or tangent,
    third point R, reflection to P+Q -- updates LIVE as you drag. Drag Q onto
    P to see doubling; click the O marker to use the identity; drag Q to P's
    mirror image to see inverses (sum = O, vertical line)."""
    return _HEAD + r"""
<div class="panel">
  <div class="modebar" style="flex-wrap:wrap;">
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin:0 2px;">f =</span>
    <button class="seg" id="rfdn">−</button><span id="rfval" style="align-self:center;min-width:34px;text-align:center;">−1.0</span><button class="seg" id="rfup">+</button>
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin:0 2px 0 10px;">g =</span>
    <button class="seg" id="rgdn">−</button><span id="rgval" style="align-self:center;min-width:34px;text-align:center;">1.0</span><button class="seg" id="rgup">+</button>
  </div>
  <div class="stage" style="grid-template-columns: 1fr;">
    <div class="cell">
      <div class="cap">y² = x³ + fx + g — click to place P, Q; drag them along the curve</div>
      <canvas id="gr" width="640" height="430"></canvas>
      <div id="grout" style="margin-top:8px;font-size:.93rem;line-height:1.7;"></div>
      <div class="hint" style="margin-top:4px;">the line through P and Q meets the curve in R; reflecting R gives P + Q · drag Q onto P for the tangent (2P) · drag Q to P's mirror for inverses (sum = 𝒪) · click 𝒪 to use the identity · click empty space to clear</div>
    </div>
  </div>
</div>
<script>
"use strict";
const INK="#d7d9dc", MUT="#9aa4ad", DOM="#4da3d8", SUP="#e0b64f", RED="#ef6f6f", GRN="#69b382", PUR="#b58fd8", OLI="#a9b665";
const XR=3.3, YR=4.5, OXY=[0,4.1];
const cv=document.getElementById("gr"), ctx=cv.getContext("2d");
const out=document.getElementById("grout");
let f=-1.0, g=1.0;
let pts=[];                 // 0..2 entries: "O" or [x,y] on the curve
let drag=null;
const PX=x=>(x+XR)/(2*XR)*cv.width, PY=y=>(YR-y)/(2*YR)*cv.height;
const y2of=x=>x*x*x+f*x+g;

function snap(mx,my){                      // nearest curve point to pixel (mx,my)
  let best=null,bd=1e18;
  for(let i=0;i<=800;i++){
    const x=-XR+2*XR*i/800, v=y2of(x);
    if(v<0) continue;
    const y=Math.sqrt(v);
    for(const s of [1,-1]){
      const dx=PX(x)-mx, dy=PY(s*y)-my, d=dx*dx+dy*dy;
      if(d<bd){bd=d; best=[x,s*y];}
    }
  }
  return {pt:best, d2:bd};
}
function sum2(P,Q){                        // -> {S, R, msg}; S = "O" | [x,y]
  if(P==="O"&&Q==="O") return {S:"O",R:null,msg:"𝒪 + 𝒪 = 𝒪."};
  if(P==="O") return {S:Q,R:null,msg:"𝒪 + Q = Q (𝒪 is the identity)."};
  if(Q==="O") return {S:P,R:null,msg:"P + 𝒪 = P (𝒪 is the identity)."};
  const dbl=Math.abs(P[0]-Q[0])<1e-9&&Math.abs(P[1]-Q[1])<1e-9;
  if(dbl&&Math.abs(P[1])<1e-6) return {S:"O",R:null,msg:"P is 2-torsion, so 2P = 𝒪."};
  if(!dbl&&Math.abs(P[0]-Q[0])<1e-9) return {S:"O",R:null,msg:"P and Q are inverses, so P + Q = 𝒪."};
  const m=dbl?(3*P[0]*P[0]+f)/(2*P[1]):(Q[1]-P[1])/(Q[0]-P[0]);
  const x3=m*m-P[0]-Q[0], yR=m*(x3-P[0])+P[1];
  return {S:[x3,-yR], R:[x3,yR], msg:null, m};
}
const fmt=P=>P==="O"?"𝒪":`(${P[0].toFixed(2)}, ${P[1].toFixed(2)})`;
function draw(){
  ctx.clearRect(0,0,cv.width,cv.height);
  ctx.strokeStyle="rgba(255,255,255,0.12)"; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(PX(-XR),PY(0)); ctx.lineTo(PX(XR),PY(0));
  ctx.moveTo(PX(0),PY(-YR)); ctx.lineTo(PX(0),PY(YR)); ctx.stroke();
  // the curve
  ctx.strokeStyle=DOM; ctx.lineWidth=2;
  for(const s of [1,-1]){
    ctx.beginPath(); let pen=false;
    for(let i=0;i<=900;i++){
      const x=-XR+2*XR*i/900, v=y2of(x);
      if(v>=0){const py=PY(s*Math.sqrt(v));
        pen?ctx.lineTo(PX(x),py):ctx.moveTo(PX(x),py); pen=true;}
      else pen=false;
    }
    ctx.stroke();
  }
  ctx.font="12px system-ui";
  const both=pts.length===2, res=both?sum2(pts[0],pts[1]):null;
  const sumIsO=res&&res.S==="O";
  // O marker
  ctx.fillStyle=(pts.includes("O"))?(pts[0]==="O"?RED:GRN):(sumIsO?SUP:OLI);
  ctx.beginPath(); ctx.arc(PX(OXY[0]),PY(OXY[1]),sumIsO?8:6,0,7); ctx.fill();
  ctx.fillStyle=OLI; ctx.fillText("𝒪", PX(OXY[0])+9, PY(OXY[1])+4);
  if(both){
    const P=pts[0], Q=pts[1];
    if(P!=="O"&&Q!=="O"){
      const dbl=Math.abs(P[0]-Q[0])<1e-9&&Math.abs(P[1]-Q[1])<1e-9;
      ctx.setLineDash([6,4]); ctx.strokeStyle=SUP; ctx.lineWidth=1.4;
      if(res.msg&&res.S==="O"){         // vertical line to O
        ctx.beginPath(); ctx.moveTo(PX(P[0]),PY(-YR)); ctx.lineTo(PX(P[0]),PY(YR)); ctx.stroke();
      } else if(!res.msg){
        const m=res.m, y0=m*(-XR-P[0])+P[1], y1=m*(XR-P[0])+P[1];
        ctx.beginPath(); ctx.moveTo(PX(-XR),PY(y0)); ctx.lineTo(PX(XR),PY(y1)); ctx.stroke();
      }
      ctx.setLineDash([]);
    }
    if(res&&!res.msg&&res.R){
      // reflection R -> sum
      ctx.setLineDash([2,4]); ctx.strokeStyle="rgba(255,255,255,0.5)"; ctx.lineWidth=1;
      ctx.beginPath(); ctx.moveTo(PX(res.R[0]),PY(res.R[1]));
      ctx.lineTo(PX(res.S[0]),PY(res.S[1])); ctx.stroke(); ctx.setLineDash([]);
      ctx.fillStyle=PUR; ctx.beginPath(); ctx.arc(PX(res.R[0]),PY(res.R[1]),5.5,0,7); ctx.fill();
      ctx.fillText("R", PX(res.R[0])+8, PY(res.R[1])-6);
      const dbl=Math.abs(pts[0][0]-pts[1][0])<1e-9&&Math.abs(pts[0][1]-pts[1][1])<1e-9;
      ctx.fillStyle=SUP; ctx.beginPath(); ctx.arc(PX(res.S[0]),PY(res.S[1]),7,0,7); ctx.fill();
      ctx.strokeStyle="#000"; ctx.lineWidth=1; ctx.stroke();
      ctx.fillText(dbl?"2P":"P + Q", PX(res.S[0])+9, PY(res.S[1])-7);
    }
  }
  pts.forEach((P,i)=>{
    if(P==="O") return;
    ctx.fillStyle=i?GRN:RED; ctx.beginPath(); ctx.arc(PX(P[0]),PY(P[1]),6.5,0,7); ctx.fill();
    ctx.strokeStyle="#fff"; ctx.lineWidth=1.4; ctx.stroke();
    ctx.fillStyle=i?GRN:RED; ctx.fillText(i?"Q":"P", PX(P[0])+9, PY(P[1])-7);
  });
}
function render(){
  draw();
  const disc=-16*(4*f*f*f+27*g*g);
  let s=`y² = x³ ${f<0?"−":"+"} ${Math.abs(f).toFixed(1)}x ${g<0?"−":"+"} ${Math.abs(g).toFixed(1)} &nbsp;·&nbsp; Δ = ${disc.toFixed(1)}`;
  if(Math.abs(disc)<1e-6) s+=` &nbsp;<span style="color:${SUP}">singular — step f or g</span>`;
  pts.forEach((P,i)=>{s+=` &nbsp;·&nbsp; <b style="color:${i?GRN:RED}">${i?"Q":"P"}</b> = ${fmt(P)}`;});
  if(pts.length===2){
    const r=sum2(pts[0],pts[1]);
    const dbl=pts[0]!=="O"&&pts[1]!=="O"&&Math.abs(pts[0][0]-pts[1][0])<1e-9&&Math.abs(pts[0][1]-pts[1][1])<1e-9;
    s+=` &nbsp;·&nbsp; <b style="color:${SUP}">${dbl?"2P":"P + Q"} = ${fmt(r.S)}</b>`;
    if(r.msg) s+=` <span style="color:${MUT}">— ${r.msg}</span>`;
  }
  out.innerHTML=s;
}
cv.addEventListener("pointerdown",e=>{
  const b=cv.getBoundingClientRect(), mx=(e.clientX-b.left)*cv.width/b.width,
        my=(e.clientY-b.top)*cv.height/b.height;
  cv.setPointerCapture(e.pointerId);
  // existing point hit?
  for(let i=0;i<pts.length;i++){
    const P=pts[i]; if(P==="O") continue;
    const dx=PX(P[0])-mx, dy=PY(P[1])-my;
    if(dx*dx+dy*dy<15*15){drag=i; render(); e.preventDefault(); return;}
  }
  // O marker?
  {const dx=PX(OXY[0])-mx, dy=PY(OXY[1])-my;
   if(dx*dx+dy*dy<14*14){ if(pts.length>=2)pts=[]; pts.push("O"); render(); e.preventDefault(); return;}}
  const sn=snap(mx,my);
  if(sn.pt&&sn.d2<28*28){
    if(pts.length>=2) pts=[];
    pts.push(sn.pt); drag=pts.length-1;
  } else if(pts.length) pts=[];
  render(); e.preventDefault();
});
cv.addEventListener("pointermove",e=>{
  if(drag===null||pts[drag]==="O")return;
  const b=cv.getBoundingClientRect(), mx=(e.clientX-b.left)*cv.width/b.width,
        my=(e.clientY-b.top)*cv.height/b.height;
  const sn=snap(mx,my);
  if(sn.pt){
    // snap onto the OTHER point when close (doubling / clean coincidence)
    const o=pts[1-drag];
    if(o&&o!=="O"&&Math.abs(PX(sn.pt[0])-PX(o[0]))<6&&Math.abs(PY(sn.pt[1])-PY(o[1]))<6)
      pts[drag]=[o[0],o[1]];
    else pts[drag]=sn.pt;
    render();
  }
});
window.addEventListener("pointerup",()=>{drag=null;});
const clampfg=v=>Math.max(-5,Math.min(5,Math.round(v*10)/10));
function syncFG(){document.getElementById("rfval").textContent=f.toFixed(1);
  document.getElementById("rgval").textContent=g.toFixed(1); pts=pts.filter(P=>P==="O"); render();}
document.getElementById("rfup").onclick=()=>{f=clampfg(f+0.5);syncFG();};
document.getElementById("rfdn").onclick=()=>{f=clampfg(f-0.5);syncFG();};
document.getElementById("rgup").onclick=()=>{g=clampfg(g+0.5);syncFG();};
document.getElementById("rgdn").onclick=()=>{g=clampfg(g-0.5);syncFG();};
render();
</script>
"""


def group_law_fp_html() -> str:
    """S1 applet: the same chord-tangent law over F_p, live.

    Replaces the p selectbox + f/g number inputs + plotly click-grid +
    Clear/Show line/Compute sum buttons. Click two curve points (same point
    twice = doubling); the 'line' through them -- all its F_p points -- plus
    R and the reflected sum appear instantly. O sits beyond the top-right
    corner and is clickable; vertical lines (inverses, 2-torsion tangents)
    pass through it."""
    return _HEAD + r"""
<div class="panel">
  <div class="modebar" style="flex-wrap:wrap;">
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin-right:4px;">p =</span>
    <button class="seg gpbtn" data-p="7">7</button>
    <button class="seg gpbtn" data-p="11">11</button>
    <button class="seg gpbtn on" data-p="17">17</button>
    <button class="seg gpbtn" data-p="23">23</button>
    <button class="seg gpbtn" data-p="31">31</button>
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin:0 2px 0 14px;">f =</span>
    <button class="seg" id="qfdn">−</button><span id="qfval" style="align-self:center;min-width:26px;text-align:center;">0</span><button class="seg" id="qfup">+</button>
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin:0 2px 0 10px;">g =</span>
    <button class="seg" id="qgdn">−</button><span id="qgval" style="align-self:center;min-width:26px;text-align:center;">1</span><button class="seg" id="qgup">+</button>
  </div>
  <div class="stage" style="grid-template-columns: 1fr;">
    <div class="cell">
      <div class="cap">E(𝔽ₚ) — click P, then Q (same point twice = 2P)</div>
      <canvas id="gq" width="480" height="480"></canvas>
      <div id="gqout" style="margin-top:8px;font-size:.93rem;line-height:1.7;"></div>
      <div class="hint" style="margin-top:4px;">gold dots: the 𝔽ₚ-points of the line through P and Q (mod p it wraps!) · purple: the third intersection R · gold ring: P + Q = R reflected · vertical lines run to 𝒪 · a third click clears</div>
    </div>
  </div>
</div>
<script>
"use strict";
const INK="#d7d9dc", MUT="#9aa4ad", DOM="#4da3d8", SUP="#e0b64f", RED="#ef6f6f", GRN="#69b382", PUR="#b58fd8", OLI="#a9b665";
const cv=document.getElementById("gq"), ctx=cv.getContext("2d");
const out=document.getElementById("gqout");
let p=17, f=0, g=1;
let sel=[];                              // 0..2 entries: "O" or [x,y] (raw 0..p-1)
const mod=(a,m)=>((a%m)+m)%m;
function powmod(b,e,m){let r=1n;b=BigInt(mod(b,m));let E=BigInt(e);const M=BigInt(m);
  while(E>0n){if(E&1n)r=r*b%M; b=b*b%M; E>>=1n;} return Number(r);}
const inv=(a,m)=>powmod(a,m-2,m);
const sym=x=>{const r=mod(x,p); return 2*r<p?r:r-p;};
const rhs=x=>mod(x*x*x+f*x+g,p);
function curvePts(){
  const sq={};
  for(let y=0;y<p;y++){const v=mod(y*y,p);(sq[v]=sq[v]||[]).push(y);}
  const pts=[];
  for(let x=0;x<p;x++) for(const y of (sq[rhs(x)]||[])) pts.push([x,y]);
  return pts;
}
function linePts(P,Q){
  if(P==="O"&&Q==="O") return [];
  if(P==="O"||Q==="O"){const F=P==="O"?Q:P;
    return Array.from({length:p},(_,y)=>[F[0],y]);}
  const same=mod(P[0]-Q[0],p)===0&&mod(P[1]-Q[1],p)===0;
  if(same){
    if(mod(P[1],p)===0) return Array.from({length:p},(_,y)=>[P[0],y]);
    const m=mod((3*P[0]*P[0]+f)*inv(mod(2*P[1],p),p),p);
    const c=mod(P[1]-m*P[0],p);
    return Array.from({length:p},(_,x)=>[x,mod(m*x+c,p)]);
  }
  if(mod(P[0]-Q[0],p)===0) return Array.from({length:p},(_,y)=>[P[0],y]);
  const m=mod((Q[1]-P[1])*inv(mod(Q[0]-P[0],p),p),p);
  const c=mod(P[1]-m*P[0],p);
  return Array.from({length:p},(_,x)=>[x,mod(m*x+c,p)]);
}
function sum2(P,Q){
  if(P==="O"&&Q==="O") return {S:"O",R:null,msg:"𝒪 + 𝒪 = 𝒪."};
  if(P==="O") return {S:Q,R:null,msg:"𝒪 + Q = Q (𝒪 is the identity)."};
  if(Q==="O") return {S:P,R:null,msg:"P + 𝒪 = P (𝒪 is the identity)."};
  const same=mod(P[0]-Q[0],p)===0&&mod(P[1]-Q[1],p)===0;
  if(same&&mod(P[1],p)===0) return {S:"O",R:null,msg:"P = Q is 2-torsion, so 2P = 𝒪."};
  if(!same&&mod(P[0]-Q[0],p)===0) return {S:"O",R:null,msg:"P and Q are inverses, so P + Q = 𝒪."};
  const m=same?mod((3*P[0]*P[0]+f)*inv(mod(2*P[1],p),p),p)
              :mod((Q[1]-P[1])*inv(mod(Q[0]-P[0],p),p),p);
  const x3=mod(m*m-P[0]-Q[0],p), yR=mod(m*(x3-P[0])+P[1],p);
  return {S:[x3,mod(-yR,p)], R:[x3,yR], msg:null};
}
const h=()=>Math.floor(p/2);
const M=34;                              // margin for the O corner
const GX=x=>M+ (sym(x)+h())/(p-1)*(cv.width-2*M);
const GY=y=>cv.height-M-(sym(y)+h())/(p-1)*(cv.height-2*M);
const OPX=()=>cv.width-14, OPY=()=>14;
const fmt=P=>P==="O"?"𝒪":`(${sym(P[0])}, ${sym(P[1])})`;
function render(){
  ctx.clearRect(0,0,cv.width,cv.height);
  const disc=mod(-16*(4*f*f*f+27*g*g),p);
  const pts=disc===0?[]:curvePts();
  ctx.fillStyle="rgba(255,255,255,0.10)";
  for(let x=0;x<p;x++)for(let y=0;y<p;y++){
    ctx.beginPath(); ctx.arc(GX(x),GY(y),1.3,0,7); ctx.fill();}
  const both=sel.length===2, res=both?sum2(sel[0],sel[1]):null;
  if(both){
    ctx.fillStyle="rgba(224,182,79,0.55)";
    for(const [x,y] of linePts(sel[0],sel[1])){
      ctx.beginPath(); ctx.arc(GX(x),GY(y),Math.max(2.4,42/p),0,7); ctx.fill();}
  }
  ctx.fillStyle=DOM;
  for(const [x,y] of pts){
    ctx.beginPath(); ctx.arc(GX(x),GY(y),Math.max(3,66/p),0,7); ctx.fill();}
  ctx.font="12px system-ui";
  if(res&&!res.msg){
    ctx.fillStyle=PUR; ctx.beginPath();
    ctx.arc(GX(res.R[0]),GY(res.R[1]),Math.max(4.5,80/p),0,7); ctx.fill();
    ctx.fillText("R",GX(res.R[0])+8,GY(res.R[1])-6);
    ctx.fillStyle=SUP; ctx.beginPath();
    ctx.arc(GX(res.S[0]),GY(res.S[1]),Math.max(5.5,90/p),0,7); ctx.fill();
    ctx.strokeStyle="#000"; ctx.lineWidth=1; ctx.stroke();
    ctx.fillText(fmt(sel[0])===fmt(sel[1])?"2P":"P + Q",GX(res.S[0])+9,GY(res.S[1])-7);
  }
  sel.forEach((P,i)=>{
    if(P==="O")return;
    ctx.fillStyle=i?GRN:RED; ctx.beginPath();
    ctx.arc(GX(P[0]),GY(P[1]),Math.max(5.5,90/p),0,7); ctx.fill();
    ctx.strokeStyle="#fff"; ctx.lineWidth=1.4; ctx.stroke();
    ctx.fillStyle=i?GRN:RED; ctx.fillText(i?"Q":"P",GX(P[0])+9,GY(P[1])-7);
  });
  // O marker top-right
  const sumIsO=res&&res.S==="O";
  ctx.fillStyle=sel.includes("O")?(sel[0]==="O"?RED:GRN):(sumIsO?SUP:OLI);
  ctx.beginPath(); ctx.arc(OPX(),OPY(),sumIsO?8:6,0,7); ctx.fill();
  ctx.fillStyle=OLI; ctx.fillText("𝒪",OPX()-22,OPY()+5);
  let s=`E: y² = x³ + ${f}x + ${g} (mod ${p})`;
  if(disc===0) s+=` &nbsp;<span style="color:${SUP}">singular — step f or g</span>`;
  else s+=` &nbsp;·&nbsp; #E = ${pts.length+1}`;
  sel.forEach((P,i)=>{s+=` &nbsp;·&nbsp; <b style="color:${i?GRN:RED}">${i?"Q":"P"}</b> = ${fmt(P)}`;});
  if(res){
    const dbl=sel[0]!=="O"&&sel[1]!=="O"&&mod(sel[0][0]-sel[1][0],p)===0&&mod(sel[0][1]-sel[1][1],p)===0;
    s+=` &nbsp;·&nbsp; <b style="color:${SUP}">${dbl?"2P":"P + Q"} = ${fmt(res.S)}</b>`;
    if(res.msg) s+=` <span style="color:${MUT}">— ${res.msg}</span>`;
  }
  out.innerHTML=s;
}
cv.addEventListener("pointerdown",e=>{
  const b=cv.getBoundingClientRect(), mx=(e.clientX-b.left)*cv.width/b.width,
        my=(e.clientY-b.top)*cv.height/b.height;
  {const dx=OPX()-mx, dy=OPY()-my;
   if(dx*dx+dy*dy<15*15){ if(sel.length>=2)sel=[]; sel.push("O"); render(); return;}}
  const disc=mod(-16*(4*f*f*f+27*g*g),p);
  if(disc===0) return;
  let best=null,bd=16*16;
  for(const q of curvePts()){
    const dx=GX(q[0])-mx, dy=GY(q[1])-my;
    if(dx*dx+dy*dy<bd){bd=dx*dx+dy*dy; best=q;}
  }
  if(best){ if(sel.length>=2)sel=[]; sel.push(best); }
  else sel=[];
  render();
});
function syncFG(){document.getElementById("qfval").textContent=f;
  document.getElementById("qgval").textContent=g; sel=[]; render();}
document.getElementById("qfup").onclick=()=>{f=mod(f+1,p);syncFG();};
document.getElementById("qfdn").onclick=()=>{f=mod(f-1,p);syncFG();};
document.getElementById("qgup").onclick=()=>{g=mod(g+1,p);syncFG();};
document.getElementById("qgdn").onclick=()=>{g=mod(g-1,p);syncFG();};
document.querySelectorAll(".gpbtn").forEach(b=>b.addEventListener("click",()=>{
  p=+b.dataset.p; f=mod(f,p); g=mod(g,p); sel=[];
  document.querySelectorAll(".gpbtn").forEach(z=>z.classList.toggle("on",z===b));
  syncFG();
}));
render();
</script>
"""


# Shared JS: tiny marching-squares contour tracer for the static implicit
# pictures (real curve galleries, quadratic-form level sets, cubic pairs).
_MS_JS = r"""
function contour(ctx,cv,x0,x1,y0,y1,fn,col,lw){
  const N=110, vals=new Float64Array((N+1)*(N+1));
  for(let iy=0;iy<=N;iy++)for(let ix=0;ix<=N;ix++)
    vals[iy*(N+1)+ix]=fn(x0+(x1-x0)*ix/N, y1-(y1-y0)*iy/N);
  const PX=x=>(x-x0)/(x1-x0)*cv.width, PY=y=>(y1-y)/(y1-y0)*cv.height;
  const gx=ix=>x0+(x1-x0)*ix/N, gy=iy=>y1-(y1-y0)*iy/N;
  ctx.strokeStyle=col; ctx.lineWidth=lw; ctx.beginPath();
  for(let iy=0;iy<N;iy++)for(let ix=0;ix<N;ix++){
    const a=vals[iy*(N+1)+ix], b=vals[iy*(N+1)+ix+1],
          c=vals[(iy+1)*(N+1)+ix+1], d=vals[(iy+1)*(N+1)+ix];
    const pts=[];
    if(a*b<0) pts.push([gx(ix)+(gx(ix+1)-gx(ix))*a/(a-b), gy(iy)]);
    if(b*c<0) pts.push([gx(ix+1), gy(iy)+(gy(iy+1)-gy(iy))*b/(b-c)]);
    if(d*c<0) pts.push([gx(ix)+(gx(ix+1)-gx(ix))*d/(d-c), gy(iy+1)]);
    if(a*d<0) pts.push([gx(ix), gy(iy)+(gy(iy+1)-gy(iy))*a/(a-d)]);
    for(let k=0;k+1<pts.length;k+=2){
      ctx.moveTo(PX(pts[k][0]),PY(pts[k][1]));
      ctx.lineTo(PX(pts[k+1][0]),PY(pts[k+1][1]));
    }
  }
  ctx.stroke();
}
function axes(ctx,cv,x0,x1,y0,y1){
  const PX=x=>(x-x0)/(x1-x0)*cv.width, PY=y=>(y1-y)/(y1-y0)*cv.height;
  ctx.strokeStyle="rgba(255,255,255,0.14)"; ctx.lineWidth=1; ctx.beginPath();
  ctx.moveTo(PX(x0),PY(0)); ctx.lineTo(PX(x1),PY(0));
  ctx.moveTo(PX(0),PY(y0)); ctx.lineTo(PX(0),PY(y1)); ctx.stroke();
}
"""


def real_examples_html() -> str:
    """S1 static gallery (replaces fig_ex): five familiar algebraic curves."""
    return _HEAD + r"""
<div class="panel">
  <div class="stage" style="grid-template-columns: repeat(5, 1fr); gap:8px;">
    <div class="cell"><div class="cap">Line<br>y = x</div><canvas class="exc" data-i="0" width="120" height="120"></canvas></div>
    <div class="cell"><div class="cap">Parabola<br>y = x²</div><canvas class="exc" data-i="1" width="120" height="120"></canvas></div>
    <div class="cell"><div class="cap">Hyperbola<br>xy = 1</div><canvas class="exc" data-i="2" width="120" height="120"></canvas></div>
    <div class="cell"><div class="cap">Circle<br>x²+y² = 1</div><canvas class="exc" data-i="3" width="120" height="120"></canvas></div>
    <div class="cell"><div class="cap">Ellipse<br>x²/4+y² = 1</div><canvas class="exc" data-i="4" width="120" height="120"></canvas></div>
  </div>
</div>
<script>
"use strict";
""" + _MS_JS + r"""
const FNS=[(x,y)=>y-x,(x,y)=>y-x*x,(x,y)=>x*y-1,(x,y)=>x*x+y*y-1,(x,y)=>x*x/4+y*y-1];
document.querySelectorAll(".exc").forEach(cv=>{
  const ctx=cv.getContext("2d");
  axes(ctx,cv,-2.2,2.2,-2.2,2.2);
  contour(ctx,cv,-2.2,2.2,-2.2,2.2,FNS[+cv.dataset.i],"#4da3d8",2);
});
</script>
"""


def fp_examples_html() -> str:
    """S1 static gallery (replaces fig17): the same five curves over F_17."""
    return _HEAD + r"""
<div class="panel">
  <div class="stage" style="grid-template-columns: repeat(5, 1fr); gap:8px;">
    <div class="cell"><div class="cap">Line<br>y = x</div><canvas class="fpc" data-i="0" width="120" height="120"></canvas></div>
    <div class="cell"><div class="cap">Parabola<br>y = x²</div><canvas class="fpc" data-i="1" width="120" height="120"></canvas></div>
    <div class="cell"><div class="cap">Hyperbola<br>xy = 1</div><canvas class="fpc" data-i="2" width="120" height="120"></canvas></div>
    <div class="cell"><div class="cap">Circle<br>x²+y² = 1</div><canvas class="fpc" data-i="3" width="120" height="120"></canvas></div>
    <div class="cell"><div class="cap">Ellipse<br>x²+4y² = 4</div><canvas class="fpc" data-i="4" width="120" height="120"></canvas></div>
  </div>
</div>
<script>
"use strict";
const p=17, h=8;
const FNS=[(x,y)=>y-x,(x,y)=>y-x*x,(x,y)=>x*y-1,(x,y)=>x*x+y*y-1,(x,y)=>x*x+4*y*y-4];
const mod=(a)=>((a%p)+p)%p;
document.querySelectorAll(".fpc").forEach(cv=>{
  const ctx=cv.getContext("2d"), fn=FNS[+cv.dataset.i];
  const G=v=>(v+h+0.5)/p*cv.width;
  ctx.fillStyle="rgba(255,255,255,0.14)";
  for(let x=-h;x<=h;x++)for(let y=-h;y<=h;y++){
    ctx.beginPath(); ctx.arc(G(x),cv.height-G(y),1,0,7); ctx.fill();}
  ctx.fillStyle="#4da3d8";
  for(let x=-h;x<=h;x++)for(let y=-h;y<=h;y++)
    if(mod(fn(x,y))===0){
      ctx.beginPath(); ctx.arc(G(x),cv.height-G(y),2.4,0,7); ctx.fill();}
});
</script>
"""


def quadratic_form_html() -> str:
    """S1 applet (replaces fig0 + a/b/c number inputs): the ellipse q = 1."""
    return _HEAD + r"""
<div class="panel">
  <div class="modebar" style="flex-wrap:wrap;">
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin:0 2px;">a =</span>
    <button class="seg" id="qadn">−</button><span id="qaval" style="align-self:center;min-width:32px;text-align:center;">2.0</span><button class="seg" id="qaup">+</button>
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin:0 2px 0 10px;">b =</span>
    <button class="seg" id="qbdn">−</button><span id="qbval" style="align-self:center;min-width:32px;text-align:center;">1.0</span><button class="seg" id="qbup">+</button>
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin:0 2px 0 10px;">c =</span>
    <button class="seg" id="qcdn">−</button><span id="qcval" style="align-self:center;min-width:32px;text-align:center;">2.0</span><button class="seg" id="qcup">+</button>
  </div>
  <div class="stage" style="grid-template-columns: 1fr;">
    <div class="cell">
      <div class="cap">level sets of q(x, y) = ax² + bxy + cy²; the ellipse q = 1 in bold</div>
      <canvas id="qfc" width="420" height="420"></canvas>
      <div id="qfout" style="margin-top:8px;font-size:.93rem;"></div>
    </div>
  </div>
</div>
<script>
"use strict";
""" + _MS_JS + r"""
const cv=document.getElementById("qfc"), ctx=cv.getContext("2d");
const out=document.getElementById("qfout");
let a=2.0, b=1.0, c=2.0;
function render(){
  ctx.clearRect(0,0,cv.width,cv.height);
  axes(ctx,cv,-2.5,2.5,-2.5,2.5);
  const disc=4*a*c-b*b, ok=a>0&&disc>0;
  const q=(x,y)=>a*x*x+b*x*y+c*y*y;
  if(ok){
    for(const [lv,al] of [[0.25,0.25],[0.5,0.32],[2,0.32],[4,0.25]])
      contour(ctx,cv,-2.5,2.5,-2.5,2.5,(x,y)=>q(x,y)-lv,`rgba(77,163,216,${al})`,1);
    contour(ctx,cv,-2.5,2.5,-2.5,2.5,(x,y)=>q(x,y)-1,"#4da3d8",2.5);
  }
  out.innerHTML=`q(x,y) = ${a.toFixed(1)}x² ${b<0?"−":"+"} ${Math.abs(b).toFixed(1)}xy ${c<0?"−":"+"} ${Math.abs(c).toFixed(1)}y²`
    +` &nbsp;·&nbsp; 4ac − b² = ${disc.toFixed(2)} &nbsp;·&nbsp; `
    +(ok?`<span style="color:#69b382">positive definite ✓ — q = 1 is an ellipse</span>`
        :`<span style="color:#e0b64f">not positive definite (need a &gt; 0 and 4ac − b² &gt; 0)</span>`);
}
const mk=(id,fn)=>document.getElementById(id).onclick=fn;
const sync=()=>{for(const [i,v] of [["qaval",a],["qbval",b],["qcval",c]])
  document.getElementById(i).textContent=v.toFixed(1); render();};
mk("qaup",()=>{a+=0.5;sync();}); mk("qadn",()=>{a-=0.5;sync();});
mk("qbup",()=>{b+=0.5;sync();}); mk("qbdn",()=>{b-=0.5;sync();});
mk("qcup",()=>{c+=0.5;sync();}); mk("qcdn",()=>{c-=0.5;sync();});
render();
</script>
"""


def cubic_pair_html(panels) -> str:
    """S1 static pair (replaces fig_sing / fig_smth): two cubics side by side.

    panels: list of two dicts {f, g, title, sing (optional [x,y])}."""
    import json as _j
    return _HEAD + r"""
<div class="panel">
  <div class="stage" style="grid-template-columns: 1fr 1fr;">
    <div class="cell"><div class="cap" id="cpt0"></div><canvas class="cpc" data-i="0" width="280" height="230"></canvas></div>
    <div class="cell"><div class="cap" id="cpt1"></div><canvas class="cpc" data-i="1" width="280" height="230"></canvas></div>
  </div>
</div>
<script>
"use strict";
""" + _MS_JS + r"""
const PANELS=""" + _j.dumps(panels) + r""";
PANELS.forEach((P,i)=>{document.getElementById("cpt"+i).textContent=P.title;});
document.querySelectorAll(".cpc").forEach(cv=>{
  const P=PANELS[+cv.dataset.i], ctx=cv.getContext("2d");
  const x0=P.win[0],x1=P.win[1],y0=P.win[2],y1=P.win[3];
  axes(ctx,cv,x0,x1,y0,y1);
  contour(ctx,cv,x0,x1,y0,y1,(x,y)=>y*y-(x*x*x+P.f*x+P.g),"#4da3d8",2);
  if(P.sing){
    const PX=x=>(x-x0)/(x1-x0)*cv.width, PY=y=>(y1-y)/(y1-y0)*cv.height;
    ctx.fillStyle="#ef6f6f"; ctx.beginPath();
    ctx.arc(PX(P.sing[0]),PY(P.sing[1]),4.5,0,7); ctx.fill();
  }
});
</script>
"""


def endo_lattice_html() -> str:
    """S2 applet (replaces fig_endo + a/b inputs): endomorphism or not.

    Step alpha = a + bi; left panel Z[i] (every alpha works), right panel
    Z[2i] (images with odd imaginary part fall OUTSIDE the lattice, red)."""
    return _HEAD + r"""
<div class="panel">
  <div class="modebar">
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin:0 2px;">α = a + bi, &nbsp; a =</span>
    <button class="seg" id="eadn">−</button><span id="eaval" style="align-self:center;min-width:24px;text-align:center;">1</span><button class="seg" id="eaup">+</button>
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin:0 2px 0 10px;">b =</span>
    <button class="seg" id="ebdn">−</button><span id="ebval" style="align-self:center;min-width:24px;text-align:center;">1</span><button class="seg" id="ebup">+</button>
  </div>
  <div class="stage" style="grid-template-columns: 1fr 1fr;">
    <div class="cell"><div class="cap">ℤ[i] — images of ×α (all land in the lattice)</div>
      <canvas id="eza" width="320" height="320"></canvas></div>
    <div class="cell"><div class="cap">ℤ[2i] — images of ×α</div>
      <canvas id="ezb" width="320" height="320"></canvas></div>
  </div>
  <div id="ezout" style="margin-top:8px;font-size:.93rem;"></div>
  <div class="hint" style="margin-top:4px;">gray: the lattice · blue: images α·z that land IN the lattice · red: images that fall OUTSIDE — α is an endomorphism exactly when there is no red (for ℤ[2i]: when b is even)</div>
</div>
<script>
"use strict";
const cva=document.getElementById("eza"), cvb=document.getElementById("ezb");
const out=document.getElementById("ezout");
let a=1, b=1;
const SRC=3;
function render(){
  // images
  const ziImg=[], zgood=[], zbad=[];
  for(let m=-SRC;m<=SRC;m++)for(let n=-SRC;n<=SRC;n++){
    ziImg.push([a*m-b*n, a*n+b*m]);
    if(Math.abs(n)<=1){                    // Z[2i] source patch: m + 2ni
      const rx=a*m-b*2*n, ry=b*m+a*2*n;
      (ry%2===0?zgood:zbad).push([rx,ry]);
    }
  }
  const R=Math.max(5,Math.min(14,Math.max(...[...ziImg,...zgood,...zbad].flat().map(Math.abs))+1));
  const draw=(cv,isZ2,good,bad)=>{
    const ctx=cv.getContext("2d");
    ctx.clearRect(0,0,cv.width,cv.height);
    const PX=x=>(x+R)/(2*R)*cv.width, PY=y=>(R-y)/(2*R)*cv.height;
    ctx.strokeStyle="rgba(255,255,255,0.12)"; ctx.lineWidth=1; ctx.beginPath();
    ctx.moveTo(PX(-R),PY(0)); ctx.lineTo(PX(R),PY(0));
    ctx.moveTo(PX(0),PY(-R)); ctx.lineTo(PX(0),PY(R)); ctx.stroke();
    ctx.fillStyle="rgba(255,255,255,0.28)";
    for(let m=-R;m<=R;m++)for(let n=-R;n<=R;n++){
      if(isZ2&&n%2!==0) continue;
      ctx.beginPath(); ctx.arc(PX(m),PY(n),1.8,0,7); ctx.fill();
    }
    ctx.fillStyle="#4da3d8";
    for(const [x,y] of good){if(Math.abs(x)<=R&&Math.abs(y)<=R){
      ctx.beginPath(); ctx.arc(PX(x),PY(y),3.6,0,7); ctx.fill();}}
    ctx.fillStyle="#ef6f6f";
    for(const [x,y] of bad){if(Math.abs(x)<=R&&Math.abs(y)<=R){
      ctx.beginPath(); ctx.arc(PX(x),PY(y),3.6,0,7); ctx.fill();}}
  };
  draw(cva,false,ziImg,[]);
  draw(cvb,true,zgood,zbad);
  const ast=`${a} ${b<0?"−":"+"} ${Math.abs(b)}i`;
  out.innerHTML=`α = ${ast} &nbsp;·&nbsp; ℤ[i]: always an endomorphism &nbsp;·&nbsp; `
    +(zbad.length? `<b style="color:#ef6f6f">ℤ[2i]: NOT an endomorphism</b> (b = ${b} odd — red images escape the lattice)`
                 : `<b style="color:#69b382">ℤ[2i]: an endomorphism ✓</b> (b = ${b} even)`);
}
const mk=(id,fn)=>document.getElementById(id).onclick=fn;
const sync=()=>{document.getElementById("eaval").textContent=a;
  document.getElementById("ebval").textContent=b; render();};
mk("eaup",()=>{a++;sync();}); mk("eadn",()=>{a--;sync();});
mk("ebup",()=>{b++;sync();}); mk("ebdn",()=>{b--;sync();});
render();
</script>
"""


def mult_frobenius_html() -> str:
    """S3 applet (replaces fig_m + p selectbox + n/k sliders): F_{p^n}^* on
    the unit circle with the full Frobenius action.

    p and n are buttons (n rebuilt per p so the point count stays <= 160);
    the orbit is chosen by CLICKING a point (was: the k slider). Points are
    coloured by their field of definition; fixed points carry loops."""
    return _HEAD + r"""
<div class="panel">
  <div class="modebar" style="flex-wrap:wrap;">
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin-right:4px;">p =</span>
    <button class="seg mpbtn" data-p="2">2</button>
    <button class="seg mpbtn on" data-p="3">3</button>
    <button class="seg mpbtn" data-p="5">5</button>
    <button class="seg mpbtn" data-p="7">7</button>
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin:0 4px 0 14px;">n =</span>
    <span id="nbtns" style="display:flex;gap:6px;"></span>
  </div>
  <div class="stage" style="grid-template-columns: 1fr;">
    <div class="cell">
      <div class="cap">𝔽<sub>pⁿ</sub>* on the unit circle — <b>click a point</b> to follow its Frobenius orbit</div>
      <canvas id="mfc" width="460" height="460"></canvas>
      <div id="mfout" style="margin-top:8px;font-size:.93rem;line-height:1.7;"></div>
      <div class="hint" style="margin-top:4px;">every point gets an arrow x ↦ xᵖ; loops mark the fixed points 𝔽ₚ* · orbit length = degree of the field the point generates</div>
    </div>
  </div>
</div>
<script>
"use strict";
const cv=document.getElementById("mfc"), ctx=cv.getContext("2d");
const out=document.getElementById("mfout");
let p=3, n=2, khi=1;
const nmax=pp=>{let m=1; while(Math.pow(pp,m+1)-1<=160)m++; return m;};
const M=()=>Math.pow(p,n)-1;
const CX=cv.width/2, CY=cv.height/2, RAD=cv.width/2-46;
const XY=k=>{const t=2*Math.PI*k/M(); return [CX+RAD*Math.cos(t), CY-RAD*Math.sin(t)];};
function orbit(k0){const seen=[k0]; let k=(p*k0)%M();
  while(k!==k0){seen.push(k); k=(p*k)%M();} return seen;}
const PAL={1:"#ef6f6f",2:"#4da3d8",3:"#69b382",4:"#e0b64f",5:"#b58fd8",6:"#5fc4c4",7:"#c9a35f"};
function render(){
  const m=M();
  if(khi>=m) khi=Math.min(1,m-1);
  ctx.clearRect(0,0,cv.width,cv.height);
  ctx.strokeStyle="rgba(255,255,255,0.18)"; ctx.lineWidth=1;
  ctx.beginPath(); ctx.arc(CX,CY,RAD,0,7); ctx.stroke();
  const deg={}; for(let k=0;k<m;k++) deg[k]=orbit(k).length;
  const hot=new Set(orbit(khi));
  // arrows
  for(let k=0;k<m;k++){
    const kk=(p*k)%m, isHot=hot.has(k);
    ctx.strokeStyle=isHot?"#fff":"rgba(255,255,255,0.22)";
    ctx.fillStyle=ctx.strokeStyle;
    ctx.lineWidth=isHot?1.8:0.9;
    if(kk===k){
      const t=2*Math.PI*k/m, r=13;
      ctx.beginPath();
      ctx.arc(CX+(RAD+r+3)*Math.cos(t), CY-(RAD+r+3)*Math.sin(t), r, 0, 7);
      ctx.stroke();
    } else {
      const A=XY(k), B=XY(kk);
      const mx=(A[0]+B[0])/2, my=(A[1]+B[1])/2;
      const dx=B[0]-A[0], dy=B[1]-A[1], L=Math.hypot(dx,dy);
      const cx2=mx-0.18*dy, cy2=my+0.18*dx;              // arc3 rad=0.18-ish
      ctx.beginPath(); ctx.moveTo(A[0],A[1]);
      ctx.quadraticCurveTo(cx2,cy2,B[0],B[1]); ctx.stroke();
      const an=Math.atan2(B[1]-cy2,B[0]-cx2);
      ctx.beginPath(); ctx.moveTo(B[0],B[1]);
      ctx.lineTo(B[0]-8*Math.cos(an-0.42),B[1]-8*Math.sin(an-0.42));
      ctx.lineTo(B[0]-8*Math.cos(an+0.42),B[1]-8*Math.sin(an+0.42));
      ctx.closePath(); ctx.fill();
    }
  }
  // points
  for(let k=0;k<m;k++){
    const [x,y]=XY(k), d=deg[k];
    ctx.fillStyle=PAL[d]||"#9aa4ad";
    ctx.beginPath(); ctx.arc(x,y,d===1?7:5.5,0,7); ctx.fill();
    ctx.strokeStyle=hot.has(k)?"#fff":"rgba(0,0,0,0.5)";
    ctx.lineWidth=hot.has(k)?2:1; ctx.stroke();
  }
  // legend
  const divs=[]; for(let d=1;d<=n;d++) if(n%d===0) divs.push(d);
  ctx.font="12px system-ui"; ctx.textAlign="left";
  divs.forEach((d,i)=>{
    ctx.fillStyle=PAL[d]||"#9aa4ad";
    ctx.beginPath(); ctx.arc(14,18+i*20,5,0,7); ctx.fill();
    ctx.fillStyle="#9aa4ad";
    ctx.fillText(d===1?`𝔽${p}* (deg 1, fixed)`:`deg ${d} (𝔽${p}^${d})`, 26, 22+i*20);
  });
  const dsel=deg[khi];
  out.innerHTML=`#𝔽<sub>${p}<sup>${n}</sup></sub>* = ${p}<sup>${n}</sup> − 1 = ${m}`
    +` &nbsp;·&nbsp; highlighted orbit: k = ${khi}, length ${dsel} — ζ<sup>${khi}</sup> generates 𝔽<sub>${p}<sup>${dsel}</sup></sub>`
    +(dsel===1?" (a fixed point: it lies in 𝔽<sub>p</sub>*)":`; applying F ${dsel} times returns it to the start`);
}
cv.addEventListener("pointerdown",e=>{
  const b=cv.getBoundingClientRect();
  const mx=(e.clientX-b.left)*cv.width/b.width, my=(e.clientY-b.top)*cv.height/b.height;
  let best=null,bd=16*16;
  for(let k=0;k<M();k++){
    const [x,y]=XY(k), dx=x-mx, dy=y-my;
    if(dx*dx+dy*dy<bd){bd=dx*dx+dy*dy; best=k;}
  }
  if(best!==null){khi=best; render();}
});
function buildN(){
  const span=document.getElementById("nbtns");
  span.innerHTML="";
  for(let v=1;v<=nmax(p);v++){
    const btn=document.createElement("button");
    btn.className="seg"+(v===n?" on":"");
    btn.textContent=v;
    btn.onclick=()=>{n=v; khi=Math.min(1,M()-1); buildN(); render();};
    span.appendChild(btn);
  }
}
document.querySelectorAll(".mpbtn").forEach(b=>b.addEventListener("click",()=>{
  p=+b.dataset.p; n=Math.min(2,nmax(p)); khi=1;
  document.querySelectorAll(".mpbtn").forEach(z=>z.classList.toggle("on",z===b));
  buildN(); render();
}));
buildN(); render();
</script>
"""


def f5_lattice_html() -> str:
    """S3 applet (replaces fig_ec + a select_slider + F_25 checkbox): the
    curves y^2 = x^3 + ax over F_5 next to their lattice models.

    a buttons pick the curve; the right panel shows the fixed points of
    multiplication by alpha on C/Z[i] (the lattice model), with an F_25
    toggle. The identity 0 is visible at the corner -- the point the
    classical picture cannot show."""
    return _HEAD + r"""
<div class="panel">
  <div class="modebar">
    <span style="align-self:center;color:var(--muted);font-size:.9rem;margin-right:4px;">y² = x³ + a·x (mod 5), &nbsp; a =</span>
    <button class="seg abtn" data-a="1">1</button>
    <button class="seg abtn" data-a="2">2</button>
    <button class="seg abtn on" data-a="3">3</button>
    <button class="seg abtn" data-a="4">4</button>
    <button class="seg" id="t25" style="margin-left:14px;">show 𝔽₂₅ points</button>
  </div>
  <div class="stage" style="grid-template-columns: 1fr 1fr;">
    <div class="cell"><div class="cap">classical: the 𝔽₅ grid</div>
      <canvas id="fla5" width="300" height="300"></canvas></div>
    <div class="cell"><div class="cap">lattice: fixed points of ×α on ℂ/ℤ[i]</div>
      <canvas id="flb5" width="300" height="300"></canvas></div>
  </div>
  <div id="fl5out" style="margin-top:8px;font-size:.93rem;line-height:1.7;"></div>
  <div class="hint" style="margin-top:4px;">the red 0 at the corner is the identity — visible in the lattice model, at infinity in the classical one · gold: the extra 𝔽₂₅ points, which the classical 𝔽₅ grid cannot show at all</div>
</div>
<script>
"use strict";
const ALPHAS={1:[1,2],2:[2,-1],3:[-2,1],4:[-1,-2]};
const cva=document.getElementById("fla5"), cvb=document.getElementById("flb5");
const out=document.getElementById("fl5out");
let a=3, show25=false;
const p=5, h=2;
function fixedPts(br,bi){
  const N=br*br+bi*bi, seen=new Set();
  for(let u=0;u<N;u++)for(let v=0;v<N;v++)
    seen.add((((u*br+v*bi)%N+N)%N)+","+(((v*br-u*bi)%N+N)%N));
  return {pts:[...seen].map(s=>s.split(",").map(Number).map(t=>t/N)), N};
}
function render(){
  const al=ALPHAS[a], b1=[al[0]-1,al[1]];
  const f5=fixedPts(b1[0],b1[1]);
  const al2=[al[0]*al[0]-al[1]*al[1], 2*al[0]*al[1]];
  const b2=[al2[0]-1,al2[1]], f25=fixedPts(b2[0],b2[1]);
  // classical panel
  {
    const ctx=cva.getContext("2d");
    ctx.clearRect(0,0,cva.width,cva.height);
    const G=v=>(v+h+0.5)/p*cva.width;
    ctx.fillStyle="rgba(255,255,255,0.15)";
    for(let x=-h;x<=h;x++)for(let y=-h;y<=h;y++){
      ctx.beginPath(); ctx.arc(G(x),cva.height-G(y),3,0,7); ctx.fill();}
    ctx.fillStyle="#4da3d8";
    let naff=0;
    for(let x=0;x<p;x++)for(let y=0;y<p;y++)
      if((y*y)%p===(((x*x*x+a*x)%p)+p)%p){
        naff++;
        const sx=2*x<p?x:x-p, sy=2*y<p?y:y-p;
        ctx.beginPath(); ctx.arc(G(sx),cva.height-G(sy),6,0,7); ctx.fill();}
  }
  // lattice panel
  {
    const ctx=cvb.getContext("2d");
    ctx.clearRect(0,0,cvb.width,cvb.height);
    const M=26, PX=v=>M+v*(cvb.width-2*M), PY=v=>cvb.height-M-v*(cvb.height-2*M);
    ctx.fillStyle="rgba(255,255,255,0.05)";
    ctx.fillRect(PX(0),PY(1),PX(1)-PX(0),PY(0)-PY(1));
    ctx.strokeStyle="rgba(255,255,255,0.25)"; ctx.lineWidth=1;
    ctx.strokeRect(PX(0),PY(1),PX(1)-PX(0),PY(0)-PY(1));
    if(show25){
      const f5set=new Set(f5.pts.map(q=>q.join(",")));
      ctx.fillStyle="#e0b64f";
      for(const q of f25.pts) if(!f5set.has(q.join(","))){
        ctx.beginPath(); ctx.arc(PX(q[0]),PY(q[1]),4,0,7); ctx.fill();}
    }
    ctx.fillStyle="#4da3d8";
    for(const q of f5.pts){
      ctx.beginPath(); ctx.arc(PX(q[0]),PY(q[1]),5.5,0,7); ctx.fill();}
    ctx.fillStyle="#ef6f6f";
    ctx.beginPath(); ctx.arc(PX(0),PY(0),7,0,7); ctx.fill();
    ctx.font="13px system-ui"; ctx.fillText("0", PX(0)+9, PY(0)-6);
  }
  const fmt=z=>`${z[0]} ${z[1]<0?"−":"+"} ${Math.abs(z[1])===1?"":Math.abs(z[1])}i`;
  out.innerHTML=`α = ${fmt(ALPHAS[a])} &nbsp;·&nbsp; #E(𝔽₅) = |α − 1|² = ${f5.N}`
    +(show25?` &nbsp;·&nbsp; <span style="color:#e0b64f">#E(𝔽₂₅) = |α² − 1|² = ${f25.N}</span>`:"");
}
document.querySelectorAll(".abtn").forEach(b=>b.addEventListener("click",()=>{
  a=+b.dataset.a;
  document.querySelectorAll(".abtn").forEach(z=>z.classList.toggle("on",z===b));
  render();
}));
document.getElementById("t25").addEventListener("click",e=>{
  show25=!show25; e.target.classList.toggle("on",show25); render();
});
render();
</script>
"""




def frobenius_flow_html() -> str:
    """S3 animation: the continuous Frobenius flow, torus and galaxy views.

    Pinned to y^2 = x^3 + 3x over F_5 -- lattice Lambda = Z[i], Frobenius the
    lift z -> alpha*z with alpha = -2 + i (|alpha| = sqrt 5, arg alpha ~ 153.4
    deg). The extension points E(F_5^n) = ker(alpha^n - 1) are drawn together
    in one picture. The flow interpolates Frobenius continuously,
    z(phi) = alpha^phi * z: every point spirals (scale |alpha|^phi, rotate
    phi*arg alpha) and at phi = 1 lands exactly on alpha*z, so the whole set
    maps to itself and the loop closes seamlessly.

    Two renderings share the one flow:
      - TORUS: the fundamental cell of C/Lambda, drawn *centred on 0* so
        Frobenius reads as a rotation-and-scaling about the middle of the frame.
      - GALAXY: the multiplicative model C*/q^Z via w = exp(2 pi i z), on a
        rebased basis (tau' = 0.4 + 0.2 i, |q| ~ 0.285) so the annulus is
        legible; straight/spiral paths upstairs become log-spiral arcs.

    Each point keeps a fixed colour taken from a continuous wheel (hue = its
    angle about 0 on the centred torus), a faint shadow marks where it started,
    and a fading spiral trail shows the path it has taken this beat. Because
    Frobenius is multiplication by the fixed alpha, it rotates *every* colour by
    the same arg alpha -- the geometric signature of CM. Click a point to follow
    one Frobenius orbit: a white marker walks it, returning after (degree) beats.
    Starts paused."""
    return _HEAD + r"""
<div class="panel">
  <div class="modebar" style="flex-wrap:wrap;align-items:center;gap:8px;">
    <button class="seg" id="ffplay">&#9654; play</button>
    <span style="color:var(--muted);font-size:.9rem;margin-left:8px;">view</span>
    <button class="seg vwbtn on" data-v="torus">torus</button>
    <button class="seg vwbtn" data-v="galaxy">galaxy</button>
    <span style="color:var(--muted);font-size:.9rem;margin-left:8px;">up to</span>
    <button class="seg lvbtn" data-l="1">&#120125;&#8325;</button>
    <button class="seg lvbtn on" data-l="2">&#120125;&#8322;&#8325;</button>
    <button class="seg lvbtn" data-l="3">&#120125;&#8321;&#8322;&#8325;</button>
    <span style="color:var(--muted);font-size:.9rem;margin-left:8px;">speed</span>
    <input type="range" id="ffspeed" min="1" max="8" value="6" style="width:96px;align-self:center;">
  </div>
  <div class="stage" style="grid-template-columns:1fr;">
    <div class="cell">
      <div class="cap">𝔼(𝔽₅ⁿ) under the Frobenius flow z ↦ αᵠ·z, α = −2+i &nbsp;·&nbsp; colour = angle about 0 &nbsp;·&nbsp; <b>click a point</b> to follow its orbit</div>
      <canvas id="ffc" width="460" height="460" style="max-width:440px;margin:0 auto;"></canvas>
    </div>
  </div>
  <div id="ffout" class="info" style="margin-top:10px;"></div>
  <div class="hint" style="margin-top:4px;">faint dot = where each point rests (φ = 0); trail = the spiral it has swept this beat; at φ = 1 every point has reached its Frobenius image αz, landing on the resting spot of another — its colour arg α ≈ 153° behind, the same turn for all points</div>
</div>
<script>
"use strict";
const cv=document.getElementById("ffc"), ctx=cv.getContext("2d");
const out=document.getElementById("ffout");
const MUT="#9aa4ad";
const P=5, ALPHA=[-2,1];                        // Frobenius multiplier alpha = -2 + i
const TAU0=[0,1];                               // tau0 = i  (Lambda = Z + i Z)
const THETA=Math.atan2(ALPHA[1],ALPHA[0]);      // arg alpha ~ 2.678 rad
const RHO=Math.hypot(ALPHA[0],ALPHA[1]);        // |alpha| = sqrt 5
// galaxy (rebased) basis: tau' = 0.4 + 0.2 i, denominator c*tau0 + d = -2 + i
const TAUP=[0.4,0.2], DEN=[-2,1], QABS=Math.exp(-2*Math.PI*TAUP[1]);
const amul=(z,w)=>[z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]];
const apow=phi=>[Math.pow(RHO,phi)*Math.cos(THETA*phi), Math.pow(RHO,phi)*Math.sin(THETA*phi)];
const r01=x=>x-Math.floor(x), r05=x=>x-Math.round(x);   // reduce to [0,1) / [-.5,.5)

// ---- point set: union of ker(alpha^n - 1), n = 1..level -------------------
let level=2, pts=[], deg=[], hue=[], perm=[], PERIOD=1, baseTab=[[]];
function fixedPts(b){
  const br=b[0], bi=b[1], N=br*br+bi*bi, seen=new Set(), a=[];
  for(let u=0;u<N;u++)for(let v=0;v<N;v++){
    const x=((u*br+v*bi)%N+N)%N, y=((v*br-u*bi)%N+N)%N, key=x+","+y;
    if(!seen.has(key)){ seen.add(key); a.push([x/N,y/N]); }
  }
  return a;
}
const keyOf=q=>Math.round(q[0]*720000)+","+Math.round(q[1]*720000);
function build(){
  pts=[]; deg=[]; hue=[]; const idx=new Map();
  let an=[1,0];
  for(let n=1;n<=level;n++){
    an=amul(an,ALPHA);
    for(const q of fixedPts([an[0]-1,an[1]])){
      const k=keyOf(q);
      if(!idx.has(k)){ idx.set(k,pts.length); pts.push(q); deg.push(n); }
    }
  }
  // continuous colour wheel: hue = angle of the resting point about 0 (centred cell)
  hue=pts.map(q=>{ const uc=r05(q[0]), vc=r05(q[1]);
    return (uc===0&&vc===0)?-1 : (Math.atan2(vc,uc)*180/Math.PI+360)%360; });  // -1 = identity
  // Frobenius as an index permutation z -> alpha*z
  perm=pts.map(q=>{ const w=amul(ALPHA,q); return idx.get(keyOf([r01(w[0]),r01(w[1])])); });
  // loop period = lcm of the field degrees present: Frobenius^PERIOD = identity,
  // so after PERIOD beats every point (and its colour) is home -> a seamless loop.
  const glcm=(x,y)=>{const g=(a,b)=>b?g(b,a%b):a; return x/g(x,y)*y;};
  PERIOD=[...new Set(deg)].reduce((a,d)=>glcm(a,d),1);
  // baseTab[n][i] = the point that dot i sits on at beat n as it walks its orbit
  baseTab=[pts.map((_,i)=>i)];
  for(let n=1;n<PERIOD;n++) baseTab.push(baseTab[n-1].map(j=>perm[j]));
  sel=null; orbit=[]; orbitSet=new Set();
}

// ---- projections: cover coords z (in basis 1, tau0=i) -> screen -----------
const W=cv.width, H=cv.height, CX=W/2, CY=H/2, MG=40;
const S=W-2*MG, R=W/2-MG;                        // torus half-cell scale / galaxy radius
// each projection returns [X, Y, cellA, cellB]: the integer cell tells the
// trail loop exactly when a segment crosses a wrap boundary (skip it).
function toTorus(z){ const ca=Math.round(z[0]), cb=Math.round(z[1]);
  return [CX + (z[0]-ca)*S, CY - (z[1]-cb)*S, ca, cb]; }
function toGalaxy(z){
  const [zr,zi]=z, [dr,di]=DEN, n2=dr*dr+di*di;   // z' = z / (c tau0 + d)
  const wr=(zr*dr+zi*di)/n2, wi=(zi*dr-zr*di)/n2;
  const vs=wi/TAUP[1], vp=r01(vs), up=r01(wr - vs*TAUP[0]);   // reduce each (1, tau') coord independently
  const rr=R*Math.exp(-2*Math.PI*vp*TAUP[1]), th=2*Math.PI*(up+vp*TAUP[0]);
  return [CX + rr*Math.cos(th), CY - rr*Math.sin(th), 0, 0];
}
let view="torus";
const project=z=> view==="torus" ? toTorus(z) : toGalaxy(z);
const fill=(i,l,a)=> hue[i]<0 ? `rgba(235,235,235,${a})` : `hsla(${hue[i]},68%,${l}%,${a})`;

// ---- strip coords for continuous (unwrapped) trails -----------------------
// (u, v): fractional coords in the fundamental domain BEFORE reduction. Trails
// are drawn relative to the reduced head, so a point that leaves one edge keeps
// going (spiralling into the galaxy core / out past the rim) instead of resetting.
function stripOf(z){
  if(view==="torus") return [z[0], z[1]];
  const [zr,zi]=z, [dr,di]=DEN, n2=dr*dr+di*di;   // z' = z / (c tau0 + d)
  const wr=(zr*dr+zi*di)/n2, wi=(zi*dr-zr*di)/n2, vs=wi/TAUP[1];
  return [wr - vs*TAUP[0], vs];
}
function headRed(u,v){                              // reduce so it coincides with the dot
  if(view==="torus") return [r05(u), r05(v)];
  return [r01(u), r01(v)];
}
function stripToScreen(u,v){
  if(view==="torus") return [CX + u*S, CY - v*S];
  const rr=R*Math.exp(-2*Math.PI*v*TAUP[1]), th=2*Math.PI*(u+v*TAUP[0]);
  return [CX + rr*Math.cos(th), CY - rr*Math.sin(th)];
}

// ---- frames ---------------------------------------------------------------
function frame(){
  ctx.strokeStyle="rgba(255,255,255,0.11)"; ctx.lineWidth=1;
  if(view==="torus"){
    for(const [dx,dy] of [[-1,0],[1,0],[0,-1],[0,1]])
      ctx.strokeRect(CX+(dx-0.5)*S, CY-(dy+0.5)*S, S, S);
    ctx.strokeStyle="rgba(255,255,255,0.28)"; ctx.lineWidth=1.4;
    ctx.strokeRect(CX-0.5*S, CY-0.5*S, S, S);
  } else {
    ctx.beginPath(); ctx.arc(CX,CY,R,0,7); ctx.stroke();          // |w| = 1
    ctx.beginPath(); ctx.arc(CX,CY,R*QABS,0,7); ctx.stroke();     // |w| = |q|
    ctx.beginPath();                                              // spiral edge u' = 0
    for(let vv=0;vv<=1.0001;vv+=0.004){
      const rr=R*Math.exp(-2*Math.PI*vv*TAUP[1]), th=2*Math.PI*vv*TAUP[0];
      const x=CX+rr*Math.cos(th), y=CY-rr*Math.sin(th);
      vv?ctx.lineTo(x,y):ctx.moveTo(x,y);
    }
    ctx.stroke();
  }
}

// ---- animation ------------------------------------------------------------
let tau=0, playing=false, last=null, speed=0.06, sel=null, orbit=[], orbitSet=new Set();
const TRAILMAX=52, TRAILSTEP=0.010, TRAILSTEPS=22;  // comet tail capped by SCREEN length, so fast points don't balloon
function render(){
  const phi=tau-Math.floor(tau), beat=Math.floor(tau);
  const bt=baseTab[((beat%PERIOD)+PERIOD)%PERIOD];   // each dot's orbit position this beat
  const af=apow(phi);
  ctx.clearRect(0,0,W,H);
  frame();
  // shadows: resting positions (phi = 0)
  for(let i=0;i<pts.length;i++){
    const s=project(pts[i]);
    ctx.beginPath(); ctx.arc(s[0],s[1],3.4,0,7); ctx.fillStyle=fill(i,52,0.30); ctx.fill();
  }
  // trails: walk backward from the head (the dot) along the flow, drawn as a
  // continuous displacement from the reduced head -- so the tail crosses a
  // domain edge without resetting -- and stopped once it has run TRAILMAX
  // pixels, so a fast point's tail is no longer than a slow one's. Fades to 0.
  if(playing || phi>0.002){
    for(let i=0;i<pts.length;i++){
      const base=pts[bt[i]];                                // walk this dot's orbit
      const zh=amul(apow(phi),base), sh=stripOf(zh);        // head, unreduced strip
      const hr=headRed(sh[0],sh[1]), head=stripToScreen(hr[0],hr[1]);  // = the dot
      const path=[head]; let acc=0, prev=head;
      for(let k=1;k<=TRAILSTEPS && acc<TRAILMAX;k++){
        const ss=stripOf(amul(apow(phi-k*TRAILSTEP),base));
        const cur=stripToScreen(hr[0]+(ss[0]-sh[0]), hr[1]+(ss[1]-sh[1]));
        acc+=Math.hypot(cur[0]-prev[0],cur[1]-prev[1]); path.push(cur); prev=cur;
      }
      for(let k=1;k<path.length;k++){
        const f=1-(k-1)/path.length;                        // head bright -> tail gone
        ctx.beginPath(); ctx.moveTo(path[k-1][0],path[k-1][1]); ctx.lineTo(path[k][0],path[k][1]);
        ctx.lineWidth=0.6+1.7*f; ctx.strokeStyle=fill(i,62,0.5*f*f); ctx.stroke();
      }
    }
  }
  // dots: current flowed position (each dot walks its Frobenius orbit)
  const cur=[];
  for(let i=0;i<pts.length;i++){
    const s=project(amul(af,pts[bt[i]])); cur.push(s);
    const inOrb=orbitSet.has(i), isId=hue[i]<0;
    ctx.beginPath(); ctx.arc(s[0],s[1], isId?6:4.6, 0,7); ctx.fillStyle=fill(i,58,0.98); ctx.fill();
    if(inOrb){ ctx.lineWidth=2; ctx.strokeStyle="#fff"; ctx.stroke(); }
    else if(isId){ ctx.lineWidth=1.4; ctx.strokeStyle="rgba(0,0,0,0.5)"; ctx.stroke(); }
  }
  if(sel!==null){                                 // the selected dot itself walks its orbit
    const s=cur[sel];
    ctx.beginPath(); ctx.arc(s[0],s[1],9,0,7); ctx.lineWidth=2.4; ctx.strokeStyle="#fff"; ctx.stroke();
  }
  info(phi);
}
function info(phi){
  let s=`α = −2+i &nbsp;·&nbsp; |α| = √5, arg α ≈ 153.4° &nbsp;·&nbsp; ${view} view &nbsp;·&nbsp; ${pts.length} points &nbsp;·&nbsp; loops every ${PERIOD} beat${PERIOD>1?"s":""} = LCM of the degrees`;
  if(sel!==null){
    const L=orbit.length;
    s+=`<br><span style="color:#fff">selected a degree-${deg[sel]} point — it walks a Frobenius orbit of length ${L}`
      +(L===1?" (an 𝔽₅ point: Frobenius fixes it)":`, returning home after ${L} beats`)+"</span>";
  } else {
    s+=`<br><span style="color:var(--muted)">each colour rides its Frobenius orbit and returns after (its degree) beats; the whole picture repeats every ${PERIOD}. Click a point to trace one orbit.</span>`;
  }
  out.innerHTML=s;
}

// ---- interaction ----------------------------------------------------------
cv.addEventListener("pointerdown",e=>{
  const b=cv.getBoundingClientRect();
  const mx=(e.clientX-b.left)*W/b.width, my=(e.clientY-b.top)*H/b.height;
  const beat=Math.floor(tau), bt=baseTab[((beat%PERIOD)+PERIOD)%PERIOD];
  const af=apow(tau-beat);
  let best=null,bd=16*16;
  for(let i=0;i<pts.length;i++){
    const s=project(amul(af,pts[bt[i]])), dx=s[0]-mx, dy=s[1]-my;
    if(dx*dx+dy*dy<bd){ bd=dx*dx+dy*dy; best=i; }
  }
  if(best!==null){
    sel=best; orbit=[]; orbitSet=new Set(); let k=best;
    do{ orbit.push(k); orbitSet.add(k); k=perm[k]; }while(k!==best && orbit.length<600);
    if(!playing) render();
  }
});
document.getElementById("ffplay").addEventListener("click",e=>{
  playing=!playing; e.target.innerHTML=playing?"&#10073;&#10073; pause":"&#9654; play";
  if(playing){ last=null; requestAnimationFrame(loop); }
});
document.querySelectorAll(".vwbtn").forEach(bt=>bt.addEventListener("click",()=>{
  view=bt.dataset.v;
  document.querySelectorAll(".vwbtn").forEach(z=>z.classList.toggle("on",z===bt));
  if(!playing) render();
}));
document.querySelectorAll(".lvbtn").forEach(bt=>bt.addEventListener("click",()=>{
  level=+bt.dataset.l;
  document.querySelectorAll(".lvbtn").forEach(z=>z.classList.toggle("on",z===bt));
  build(); if(!playing) render();
}));
document.getElementById("ffspeed").addEventListener("input",e=>{ speed=+e.target.value/100; });

function loop(ms){
  if(last===null) last=ms;
  tau+=(ms-last)/1000*speed; last=ms;
  if(tau>=PERIOD) tau-=PERIOD;          // seamless: Frobenius^PERIOD = identity
  render();
  if(playing) requestAnimationFrame(loop);
}
build(); render();
</script>
"""
