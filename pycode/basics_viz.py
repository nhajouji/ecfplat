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
