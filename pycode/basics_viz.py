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
