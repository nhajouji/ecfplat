"""Canvas widgets for the Explorer page.

The flagship widget is ``isogeny_graph_html``: a self-contained canvas
renderer for the l-isogeny graph of an isogeny class. It consumes the JSON
descriptors produced by ``ecqf.class_graph_descriptor`` (one per prime l) and
lays the graph out structurally:

  * crater cycles (horizontal edges) on inner circles,
  * vertical trees hanging outward in recursive, leaf-weighted angular
    sectors (the volcano),
  * disjoint components (the l ∤ cond case: cycles grouped by endomorphism
    ring) packed side by side.

Interaction: drag to pan, wheel to zoom, hover for a readout, click to
select. When ``link_base`` is given, the selected node's readout offers a
drill-down link with ``&node={i}`` appended, navigated via ``navParent``
(Streamlit's iframe sandbox blocks target=_parent, so the link is clicked
from the parent document instead) — this is the query-param navigation used
by the drill-down Explorer.

Embedded via ``st.components.v1.html`` like the Background applets; styling
matches the shared panel look (``modular_viz._HEAD``).
"""

import json

from modular_viz import _HEAD

# Shared JS: primitive reduced forms of a discriminant + total lattice-class
# count over all conductor levels (sum of h(d0 e^2) over e | cond). Small
# enough to recompute client-side, so the Hasse picker needs no server data.
# Streamlit's component iframe is sandboxed WITHOUT allow-top-navigation, so
# target=_parent anchors are silently blocked. It does allow same-origin
# scripting, so we navigate by clicking a link created in the parent document
# (relative ?query URLs resolve against the parent there). Fallback: new tab.
_NAV_JS = r"""
window.navParent = function(url){
  try{
    const doc = window.parent.document;
    const a = doc.createElement("a");
    a.href = url; doc.body.appendChild(a); a.click(); a.remove();
  }catch(e){ window.open(url, "_blank"); }
};
"""

_FORMS_JS = r"""
const gcd=(x,y)=>{x=Math.abs(x);y=Math.abs(y);while(y)[x,y]=[y,x%y];return x;};
function reducedForms(D){
  const forms=[]; const B0=((D%2)+2)%2;
  for(let B=B0;B*B<=-D/3+1e-9;B+=2){
    const AC=(B*B-D)/4;
    for(let A=Math.max(B,1);A*A<=AC;A++){
      if(AC%A) continue;
      const C=AC/A;
      if(gcd(gcd(A,B),C)!==1) continue;
      forms.push([A,B,C]);
      if(B>0&&B<A&&A<C) forms.push([A,-B,C]);
    }
  }
  return forms;
}
function condOf(d){          // d = d0 * c^2, d0 fundamental: c = max valid e
  let c=1;
  for(let e=2;e*e*3<=-d;e++){
    if(d%(e*e)===0){ const r=(((d/(e*e))%4)+4)%4;
      if(r===0||r===1) c=e; }
  }
  return c;
}
function classCount(d){      // total lattice classes across conductor levels
  const c=condOf(d), d0=d/(c*c);
  let n=0;
  for(let e=1;e<=c;e++) if(c%e===0) n+=reducedForms(d0*e*e).length;
  return n;
}
"""


def isogeny_graph_html(descs: dict, link_base: str = None,
                       link_label: str = "open curve view →",
                       height_px: int = 560) -> str:
    """Canvas l-isogeny-graph widget.

    descs: {l: class_graph_descriptor(cls, l)} — one entry per offered prime;
    pills switch between them. link_base: query-param prefix for the selected
    node's drill-down link (e.g. "?a=6&p=101&view=curve"); None hides it.
    """
    if not descs:
        raise ValueError("need at least one descriptor")
    ls = sorted(descs)
    payload = {
        "ls": ls,
        "descs": {str(l): descs[l] for l in ls},
        "linkBase": link_base,
        "linkLabel": link_label,
        "H": height_px,
    }
    return '<meta charset="utf-8">' + _HEAD + """
<div class="panel" style="max-width:860px;">
  <div class="modebar" style="flex-wrap:wrap; align-items:center;">
    <span style="font-size:.85rem; color:var(--muted); margin-right:2px;">isogeny degree</span>
    <span id="lbar" style="display:flex; gap:8px;"></span>
    <span style="flex:1;"></span>
    <button class="seg" id="fitBtn">reset view</button>
  </div>
  <canvas id="volc" width="820" height="__H__" style="touch-action:none; cursor:grab;"></canvas>
  <div class="info" id="vInfo" style="min-height:1.4em;">&nbsp;</div>
  <div class="hint" id="vHint">drag to pan · scroll to zoom · hover a node · click to select</div>
  <div class="hint" id="vLegend" style="margin-top:6px;"></div>
</div>
<script>
(() => {
"use strict";
__NAV_JS__
const DATA = __DATA__;
const ACCENT="#4da3d8", GOLD="#e0b64f", RED="#ef6f6f", INK="#d7d9dc", MUT="#9aa4ad";
const LVL_COLORS=[GOLD, ACCENT, "#8fb4c7", "#6f8a9c", "#5a7080", "#4a5c68"];
const H_EDGE="#3f7fa8", V_EDGE="#565b61";

const cv=document.getElementById("volc"), ctx=cv.getContext("2d");
const info=document.getElementById("vInfo"), legend=document.getElementById("vLegend");

let curL=DATA.ls[0];
const layouts={};                 // l -> computed layout
let view={k:1, tx:0, ty:0};       // world -> canvas transform
let hoverI=-1, selI=-1;

const vval=(l,n)=>{let v=0; while(n%l===0){n/=l; v++;} return v;};

/* ── structural layout ──────────────────────────────────────────────── */
function computeLayout(desc){
  const N=desc.nodes.length;
  const adjH=Array.from({length:N},()=>[]), adjV=Array.from({length:N},()=>[]);
  for(const e of desc.edges){
    if(e.kind==="h"){ adjH[e.s].push(e.t); adjH[e.t].push(e.s); }
    else            { adjV[e.s].push(e.t); adjV[e.t].push(e.s); }
  }
  const lev=desc.nodes.map(nd=>vval(desc.l, nd.endo_cond));
  // connected components
  const comp=new Array(N).fill(-1); let nc=0;
  for(let i=0;i<N;i++){
    if(comp[i]>=0) continue;
    const st=[i]; comp[i]=nc;
    while(st.length){ const u=st.pop();
      for(const w of adjH[u].concat(adjV[u])) if(comp[w]<0){ comp[w]=nc; st.push(w); } }
    nc++;
  }
  const pos=new Array(N);
  // subtree leaf-weights for angular allocation (v-edges only, downward)
  const kidsOf=i=>adjV[i].filter(w=>lev[w]===lev[i]+1);
  const wt=new Array(N).fill(0);
  const weigh=i=>{ const ks=kidsOf(i); if(!ks.length){ wt[i]=1; return 1; }
                   let s=0; for(const k of ks) s+=weigh(k); wt[i]=s; return s; };
  const comps=[];
  for(let c=0;c<nc;c++){
    const nodes=[]; for(let i=0;i<N;i++) if(comp[i]===c) nodes.push(i);
    const minLev=Math.min(...nodes.map(i=>lev[i]));
    const maxLev=Math.max(...nodes.map(i=>lev[i]));
    const depth=maxLev-minLev;
    let crater=nodes.filter(i=>lev[i]===minLev);
    // order the crater along its h-cycle so the rim draws as a polygon
    if(crater.length>2){
      const inCr=new Set(crater), cyc=[crater[0]]; let prev=-1;
      while(cyc.length<crater.length){
        const cur=cyc[cyc.length-1];
        const nxt=adjH[cur].find(w=>w!==prev&&inCr.has(w)&&!cyc.includes(w));
        if(nxt===undefined) break;
        cyc.push(nxt); prev=cur;
      }
      for(const i of crater) if(!cyc.includes(i)) cyc.push(i);
      crater=cyc;
    }
    for(const i of crater) weigh(i);
    const K=crater.length;
    const R0=K===1?0:Math.max(44, 15*K/Math.PI);
    const DR=64;
    const rOf=l=>R0+(l-minLev)*DR;
    const place=(i,th,wedge)=>{
      const r=rOf(lev[i]);
      pos[i]={x:r*Math.cos(th), y:r*Math.sin(th)};
      const ks=kidsOf(i);
      if(!ks.length) return;
      const tot=ks.reduce((s,k)=>s+wt[k],0);
      let a=th-wedge/2;
      for(const k of ks){
        const w=wedge*wt[k]/tot;
        place(k, a+w/2, w);
        a+=w;
      }
    };
    crater.forEach((i,k)=>place(i, -Math.PI/2+2*Math.PI*k/K, 2*Math.PI/Math.max(K,1)));
    const Rout=rOf(maxLev)+30;
    comps.push({nodes, Rout, disc:desc.nodes[nodes[0]].endo_disc,
                cond:desc.nodes[nodes[0]].endo_cond, depth});
  }
  // pack components: rows, aiming at a wide-ish aspect
  comps.sort((a,b)=>b.Rout-a.Rout || a.cond-b.cond || a.disc-b.disc);
  const totalW=comps.reduce((s,c)=>s+2*c.Rout,0);
  const targetW=Math.max(2*comps[0].Rout, Math.min(totalW, Math.sqrt(totalW*2.4*comps[0].Rout)*1.35));
  let x=0, y=0, rowH=0;
  for(const c of comps){
    const w=2*c.Rout;
    if(x>0 && x+w>targetW){ x=0; y+=rowH+34; rowH=0; }
    c.cx=x+c.Rout; c.cy=y+c.Rout;
    x+=w+18; rowH=Math.max(rowH,2*c.Rout);
    for(const i of c.nodes){ pos[i]={x:pos[i].x+c.cx, y:pos[i].y+c.cy}; }
  }
  // world bounding box
  let x0=1e9,y0=1e9,x1=-1e9,y1=-1e9;
  for(const p of pos){ x0=Math.min(x0,p.x); y0=Math.min(y0,p.y);
                       x1=Math.max(x1,p.x); y1=Math.max(y1,p.y); }
  const maxLevAll=Math.max(...lev), minLevAll=Math.min(...lev);
  return {pos, lev, comps, bbox:[x0,y0,x1,y1], N,
          multiComp:comps.length>1, minLev:minLevAll, maxLev:maxLevAll};
}

function getLayout(){ if(!layouts[curL]) layouts[curL]=computeLayout(DATA.descs[curL]); return layouts[curL]; }

/* ── view transform ─────────────────────────────────────────────────── */
function fitView(){
  const {bbox}=getLayout();
  const [x0,y0,x1,y1]=bbox, pad=42;
  const w=Math.max(x1-x0,1), h=Math.max(y1-y0,1);
  const k=Math.min((cv.width-2*pad)/w, (cv.height-2*pad)/h, 2.2);
  view={k, tx:cv.width/2-k*(x0+x1)/2, ty:cv.height/2-k*(y0+y1)/2};
}
const W2C=p=>[view.k*p.x+view.tx, view.k*p.y+view.ty];

function nodeRadius(N){ return N<=14?11:N<=40?8.5:N<=120?6:4.2; }

/* ── drawing ────────────────────────────────────────────────────────── */
function draw(){
  const desc=DATA.descs[curL], L=getLayout();
  ctx.clearRect(0,0,cv.width,cv.height);
  // component group labels — only informative when the graph splits into
  // rings of genuinely different endomorphism discriminants
  const discs=new Set(L.comps.map(c=>c.disc));
  if(L.multiComp && discs.size>1){
    ctx.fillStyle="#6c7681"; ctx.font="12px 'Avenir Next',Helvetica,sans-serif";
    ctx.textAlign="center"; ctx.textBaseline="top";
    for(const c of L.comps){
      if(c.nodes.length<2) continue;
      const [cx,cy]=W2C({x:c.cx, y:c.cy+c.Rout-16});
      ctx.fillText("disc "+c.disc, cx, cy);
    }
  }
  // edges
  for(const e of desc.edges){
    const [ax,ay]=W2C(L.pos[e.s]), [bx,by]=W2C(L.pos[e.t]);
    ctx.strokeStyle=e.kind==="h"?H_EDGE:V_EDGE;
    ctx.lineWidth=e.kind==="h"?2.2:1.6;
    ctx.beginPath(); ctx.moveTo(ax,ay); ctx.lineTo(bx,by); ctx.stroke();
  }
  // nodes
  const rN=nodeRadius(L.N);
  desc.nodes.forEach((nd,i)=>{
    const [x,y]=W2C(L.pos[i]);
    ctx.fillStyle=LVL_COLORS[Math.min(L.lev[i], LVL_COLORS.length-1)];
    ctx.beginPath(); ctx.arc(x,y,rN,0,7); ctx.fill();
    if(i===selI||i===hoverI){
      ctx.strokeStyle=i===selI?RED:INK; ctx.lineWidth=2.4;
      ctx.beginPath(); ctx.arc(x,y,rN+3.2,0,7); ctx.stroke();
    }
  });
  // labels when the graph is small
  if(L.N<=40){
    ctx.font="500 12.5px 'Avenir Next',Helvetica,sans-serif";
    ctx.fillStyle=INK; ctx.textAlign="center"; ctx.textBaseline="bottom";
    desc.nodes.forEach((nd,i)=>{
      const [x,y]=W2C(L.pos[i]);
      ctx.fillText(nodeLabel(nd), x, y-rN-4);
    });
  }
}

function nodeLabel(nd){
  if(nd.curves&&nd.curves.length){
    const js=nd.curves.map(c=>"j="+c.j);
    return js.length>2?js.slice(0,2).join(" ")+" …":js.join("  ");
  }
  return "("+nd.qf.join(",")+")";
}

function nodeReadout(nd, i){
  const tau="τ ≈ "+nd.tau[0].toFixed(3)+" + "+nd.tau[1].toFixed(3)+"i";
  let s="⟨"+nd.qf.join(", ")+"⟩ · End ring disc "+nd.endo_disc+" (cond "+nd.endo_cond+") · "+tau;
  if(nd.curves&&nd.curves.length){
    const c=nd.curves[0];
    s+=" · j = "+nd.curves.map(k=>k.j).join(", ");
    if(c.model&&c.model.length===2) s+=" · y² = x³ + "+c.model[0]+"x + "+c.model[1];
  }
  if(DATA.linkBase!=null && i!=null)
    s+=' &nbsp;<a style="color:'+ACCENT+';cursor:pointer;" onclick="navParent(&quot;'
      +DATA.linkBase+"&node="+i+'&quot;)">'+DATA.linkLabel+"</a>";
  return s;
}

function refreshInfo(){
  const desc=DATA.descs[curL];
  if(selI>=0) info.innerHTML=nodeReadout(desc.nodes[selI], selI);
  else if(hoverI>=0) info.innerHTML=nodeReadout(desc.nodes[hoverI], null);
  else{
    let s="ℓ = "+curL+" · "+desc.nodes.length+" classes · disc "+desc.disc+" (cond "+desc.cond+")";
    if(!desc.edges.length) s+=" · no ℓ-isogenies here (ℓ is inert)";
    info.innerHTML=s;
  }
}

function refreshLegend(){
  const L=getLayout();
  if(L.maxLev===L.minLev){ legend.innerHTML=""; return; }
  let s="";
  for(let l=L.minLev;l<=L.maxLev;l++){
    const col=LVL_COLORS[Math.min(l,LVL_COLORS.length-1)];
    s+='<span style="color:'+col+';">●</span> '+(l===L.minLev?"crater (surface)":"level "+(l-L.minLev))+" &nbsp; ";
  }
  legend.innerHTML=s;
}

/* ── interaction ────────────────────────────────────────────────────── */
function hitTest(px,py){
  const L=getLayout(), rN=nodeRadius(L.N)+4;
  for(let i=L.N-1;i>=0;i--){
    const [x,y]=W2C(L.pos[i]);
    if((px-x)*(px-x)+(py-y)*(py-y)<=rN*rN) return i;
  }
  return -1;
}
const evPos=ev=>{ const r=cv.getBoundingClientRect();
  return [(ev.clientX-r.left)*cv.width/r.width, (ev.clientY-r.top)*cv.height/r.height]; };

let panning=false, moved=false, px0=0, py0=0;
cv.addEventListener("pointerdown",ev=>{
  const [x,y]=evPos(ev);
  panning=true; moved=false; px0=x; py0=y;
  cv.setPointerCapture(ev.pointerId); cv.style.cursor="grabbing";
});
cv.addEventListener("pointermove",ev=>{
  const [x,y]=evPos(ev);
  if(panning){
    if(Math.abs(x-px0)+Math.abs(y-py0)>3) moved=true;
    if(moved){ view.tx+=x-px0; view.ty+=y-py0; px0=x; py0=y; draw(); }
  } else {
    const h=hitTest(x,y);
    if(h!==hoverI){ hoverI=h; cv.style.cursor=h>=0?"pointer":"grab"; draw(); refreshInfo(); }
  }
});
cv.addEventListener("pointerup",ev=>{
  cv.style.cursor="grab";
  if(panning&&!moved){
    const [x,y]=evPos(ev);
    const h=hitTest(x,y);
    selI=(h===selI)?-1:h;
    draw(); refreshInfo();
  }
  panning=false;
});
cv.addEventListener("wheel",ev=>{
  ev.preventDefault();
  const [x,y]=evPos(ev);
  const f=Math.exp(-ev.deltaY*0.0016);
  const k2=Math.min(9,Math.max(0.15,view.k*f)), g=k2/view.k;
  view={k:k2, tx:x-g*(x-view.tx), ty:y-g*(y-view.ty)};
  draw();
},{passive:false});
document.getElementById("fitBtn").addEventListener("click",()=>{ fitView(); draw(); });

/* ── ℓ pills ────────────────────────────────────────────────────────── */
const lbar=document.getElementById("lbar");
DATA.ls.forEach(l=>{
  const b=document.createElement("button");
  b.className="seg"+(l===curL?" on":""); b.textContent="ℓ = "+l;
  b.addEventListener("click",()=>{
    curL=l; hoverI=-1; selI=-1;
    document.querySelectorAll("#lbar .seg").forEach(bb=>bb.classList.toggle("on",bb===b));
    fitView(); draw(); refreshInfo(); refreshLegend();
  });
  lbar.appendChild(b);
});

fitView(); draw(); refreshInfo(); refreshLegend();
})();
</script>
""".replace("__DATA__", json.dumps(payload)).replace("__NAV_JS__", _NAV_JS).replace("__H__", str(height_px))


def hasse_picker_html(p: int, store_p_max: int = 1021, height_px: int = 300) -> str:
    """The clickable Hasse interval for a fixed prime p.

    One stem per admissible trace a (a^2 < 4p, p ∤ a, plus a = 0 for the
    supersingular class); stem height = total number of lattice classes of
    a^2 - 4p across all conductor levels (computed client-side). Clicking a
    stem navigates the parent page to ?a=..&p=.. — the class view."""
    payload = {"p": p, "storeMax": store_p_max, "H": height_px}
    return '<meta charset="utf-8">' + _HEAD + """
<div class="panel" style="max-width:860px;">
  <canvas id="hp" width="820" height="__H__" style="touch-action:none;"></canvas>
  <div class="info" id="hpInfo" style="min-height:1.4em;">&nbsp;</div>
  <div class="hint">each stem is an isogeny class over 𝔽_p · height = number of curves in it · click one to open it</div>
</div>
<script>
(() => {
"use strict";
__NAV_JS__
const DATA = __DATA__;
__FORMS_JS__
const ACCENT="#4da3d8", GOLD="#e0b64f", RED="#ef6f6f", INK="#d7d9dc", MUT="#9aa4ad";
const cv=document.getElementById("hp"), ctx=cv.getContext("2d");
const info=document.getElementById("hpInfo");
const p=DATA.p, am=Math.floor(2*Math.sqrt(p)-1e-9);

// one entry per admissible trace
const traces=[];
for(let a=-am;a<=am;a++){
  if(a!==0 && a%p===0) continue;
  const d=a*a-4*p;
  traces.push({a, d, n:classCount(d)});
}
const nMax=Math.max(...traces.map(t=>t.n));

const padL=46, padR=18, padB=40, padT=26;
const X=a=>padL+(a+am)/(2*am)*(cv.width-padL-padR);
const Y=n=>cv.height-padB-n/nMax*(cv.height-padB-padT);
let hover=-1;

function draw(){
  ctx.clearRect(0,0,cv.width,cv.height);
  // baseline + axis labels
  ctx.strokeStyle="rgba(255,255,255,0.15)"; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(padL-14,Y(0)); ctx.lineTo(cv.width-padR+6,Y(0)); ctx.stroke();
  ctx.fillStyle=MUT; ctx.font="12px system-ui"; ctx.textAlign="center";
  ctx.fillText("a = −⌊2√p⌋", X(-am), cv.height-14);
  ctx.fillText("trace of Frobenius", X(0), cv.height-14);
  ctx.fillText("a = ⌊2√p⌋", X(am), cv.height-14);
  ctx.save(); ctx.translate(13,(padT+Y(0))/2); ctx.rotate(-Math.PI/2);
  ctx.fillText("# curves", 0, 0); ctx.restore();
  traces.forEach((t,i)=>{
    const x=X(t.a), y=Y(t.n), on=i===hover;
    const col=t.a===0?GOLD:ACCENT;
    ctx.strokeStyle=on?RED:col; ctx.globalAlpha=on?1:0.55; ctx.lineWidth=on?3:2;
    ctx.beginPath(); ctx.moveTo(x,Y(0)); ctx.lineTo(x,y); ctx.stroke();
    ctx.globalAlpha=1;
    ctx.fillStyle=on?RED:col;
    ctx.beginPath(); ctx.arc(x,y,on?6:4,0,7); ctx.fill();
  });
}
function hit(px,py){
  let best=-1, bd=14*14;
  traces.forEach((t,i)=>{
    const dx=X(t.a)-px;
    const yy=Math.min(Math.max(py,Y(t.n)),Y(0));   // clamp to the stem segment
    const q=dx*dx+(py-yy)*(py-yy);
    if(q<bd){bd=q;best=i;}
  });
  return best;
}
const evPos=ev=>{const r=cv.getBoundingClientRect();
  return [(ev.clientX-r.left)*cv.width/r.width,(ev.clientY-r.top)*cv.height/r.height];};
cv.addEventListener("pointermove",ev=>{
  const [x,y]=evPos(ev); const h=hit(x,y);
  if(h!==hover){ hover=h; cv.style.cursor=h>=0?"pointer":"default"; draw();
    if(h>=0){ const t=traces[h];
      info.textContent="a = "+t.a+" · disc = a² − 4p = "+t.d+" · "+t.n+" curve"+(t.n>1?"s":"")
        +(t.a===0?" · supersingular":"")+(p>DATA.storeMax?" · structure only (p beyond curve tables)":"");
    } else info.innerHTML="&nbsp;";
  }
});
cv.addEventListener("pointerleave",()=>{hover=-1;draw();info.innerHTML="&nbsp;";});
cv.addEventListener("pointerdown",ev=>{
  const [x,y]=evPos(ev); const h=hit(x,y);
  if(h>=0) navParent("?a="+traces[h].a+"&p="+p);
});
draw();
})();
</script>
""".replace("__DATA__", json.dumps(payload)).replace("__NAV_JS__", _NAV_JS) \
   .replace("__FORMS_JS__", _FORMS_JS) \
   .replace("__H__", str(height_px))


def fd_points_html(points: list, link_base: str = None, height_px: int = 440) -> str:
    """The class's CM points in the fundamental domain — canvas scatter.

    points: [{'x','y','color','label','sub'}] — τ coordinates, the shared
    per-row palette colour, a short label (j or the form) and a hover line.
    Two views toggled by segs: the standard FD strip and the unit-disc image
    under τ ↦ −1/τ (always bounded; the trivial class lands at the origin).
    Clicking a point drills straight into the curve view via link_base."""
    payload = {"pts": points, "linkBase": link_base, "H": height_px}
    return '<meta charset="utf-8">' + _HEAD + """
<div class="panel" style="max-width:640px;">
  <div class="modebar">
    <button class="seg on" id="fdSegDisc">unit disc (τ ↦ −1/τ)</button>
    <button class="seg" id="fdSegStd">fundamental domain</button>
  </div>
  <canvas id="fdc" width="600" height="__H__" style="touch-action:none;"></canvas>
  <div class="info" id="fdInfo" style="min-height:1.4em;">&nbsp;</div>
  <div class="hint" id="fdHint">every lattice class is a point τ · hover for details__CLICKHINT__</div>
</div>
<script>
(() => {
"use strict";
__NAV_JS__
const DATA = __DATA__;
const ACCENT="#4da3d8", GOLD="#e0b64f", RED="#ef6f6f", INK="#d7d9dc", MUT="#9aa4ad";
const cv=document.getElementById("fdc"), ctx=cv.getContext("2d");
const info=document.getElementById("fdInfo");
let mode="disc", hover=-1;

// screen positions per mode
function screenPts(){
  const pts=DATA.pts;
  if(mode==="disc"){
    // w = -1/τ maps the FD into the unit disc
    const cx=cv.width/2, cy=cv.height/2, R=Math.min(cx,cy)-28;
    return pts.map(p=>{
      const n=p.x*p.x+p.y*p.y;
      return [cx+R*(-p.x/n), cy+R*(-p.y/n)];   // canvas y down = flip of Im(w)>0? draw symmetric
    });
  }
  const ymax=Math.max(2.1, ...DATA.pts.map(p=>p.y))+0.4;
  const S=Math.min((cv.width-40)/1.6, (cv.height-46)/ymax);
  const ox=cv.width/2, oy=cv.height-30;
  return pts.map(p=>[ox+S*p.x, oy-S*p.y]);
}

function drawFrame(){
  ctx.strokeStyle="rgba(255,255,255,0.22)"; ctx.lineWidth=1.4;
  if(mode==="disc"){
    const cx=cv.width/2, cy=cv.height/2, R=Math.min(cx,cy)-28;
    ctx.beginPath(); ctx.arc(cx,cy,R,0,7); ctx.stroke();
    ctx.fillStyle=MUT; ctx.font="11px system-ui"; ctx.textAlign="center";
    ctx.fillText("|−1/τ| = 1", cx, cy-R-8);
    ctx.beginPath(); ctx.arc(cx,cy,2,0,7); ctx.fillStyle=MUT; ctx.fill();
    ctx.fillText("0 (trivial class τ → ∞)", cx, cy+14);
  } else {
    const ymax=Math.max(2.1, ...DATA.pts.map(p=>p.y))+0.4;
    const S=Math.min((cv.width-40)/1.6, (cv.height-46)/ymax);
    const ox=cv.width/2, oy=cv.height-30;
    // FD boundary: |Re τ| = 1/2 walls + unit-circle floor
    ctx.beginPath();
    ctx.moveTo(ox-S/2, oy-S*Math.sqrt(3)/2);
    ctx.lineTo(ox-S/2, 24);
    ctx.moveTo(ox+S/2, oy-S*Math.sqrt(3)/2);
    ctx.lineTo(ox+S/2, 24);
    ctx.stroke();
    ctx.beginPath();
    for(let t=120;t>=60;t-=2){const r=t*Math.PI/180;
      const x=ox+S*Math.cos(r), y=oy-S*Math.sin(r);
      t===120?ctx.moveTo(x,y):ctx.lineTo(x,y);}
    ctx.stroke();
    ctx.fillStyle=MUT; ctx.font="11px system-ui"; ctx.textAlign="center";
    ctx.fillText("−½", ox-S/2, oy+14); ctx.fillText("½", ox+S/2, oy+14);
    ctx.fillText("i", ox-10, oy-S-4);
  }
}

function draw(){
  ctx.clearRect(0,0,cv.width,cv.height);
  drawFrame();
  const sp=screenPts();
  DATA.pts.forEach((p,i)=>{
    const [x,y]=sp[i], on=i===hover;
    ctx.fillStyle=p.color||ACCENT;
    ctx.beginPath(); ctx.arc(x,y,on?8:5.5,0,7); ctx.fill();
    if(on){ ctx.strokeStyle=RED; ctx.lineWidth=2.2;
      ctx.beginPath(); ctx.arc(x,y,10,0,7); ctx.stroke(); }
  });
  if(DATA.pts.length<=40){
    ctx.font="500 11.5px system-ui"; ctx.textAlign="center"; ctx.textBaseline="bottom";
    ctx.fillStyle=INK;
    DATA.pts.forEach((p,i)=>{ const [x,y]=sp[i]; ctx.fillText(p.label, x, y-9); });
  }
}
const evPos=ev=>{const r=cv.getBoundingClientRect();
  return [(ev.clientX-r.left)*cv.width/r.width,(ev.clientY-r.top)*cv.height/r.height];};
function hit(px,py){
  const sp=screenPts(); let best=-1,bd=13*13;
  sp.forEach(([x,y],i)=>{const q=(x-px)*(x-px)+(y-py)*(y-py); if(q<bd){bd=q;best=i;}});
  return best;
}
cv.addEventListener("pointermove",ev=>{
  const [x,y]=evPos(ev); const h=hit(x,y);
  if(h!==hover){ hover=h; cv.style.cursor=h>=0&&DATA.linkBase?"pointer":"default"; draw();
    info.innerHTML=h>=0?DATA.pts[h].sub:"&nbsp;"; }
});
cv.addEventListener("pointerleave",()=>{hover=-1;draw();info.innerHTML="&nbsp;";});
cv.addEventListener("pointerdown",ev=>{
  if(!DATA.linkBase) return;
  const [x,y]=evPos(ev); const h=hit(x,y);
  if(h>=0) navParent(DATA.linkBase+"&node="+h);
});
const segD=document.getElementById("fdSegDisc"), segS=document.getElementById("fdSegStd");
segD.addEventListener("click",()=>{mode="disc"; segD.classList.add("on"); segS.classList.remove("on"); hover=-1; draw();});
segS.addEventListener("click",()=>{mode="std";  segS.classList.add("on"); segD.classList.remove("on"); hover=-1; draw();});
draw();
})();
</script>
""".replace("__DATA__", json.dumps(payload)).replace("__NAV_JS__", _NAV_JS) \
   .replace("__H__", str(height_px)) \
   .replace("__CLICKHINT__", " · click to open the curve" if link_base else "")


def curve_torus_html(qf, a: int, p: int, frobmat, n_points: int,
                     height_px: int = 380) -> str:
    """The curve view's picture: E as ℂ/Λ with E(𝔽_p) = Fix(α), plus ×α.

    Left: the fundamental cell of Λ = ⟨1, τ⟩ with the N = p+1−a fixed points
    of the Frobenius lift (the group E(𝔽_p) drawn on the torus — the CM-torus
    picture from the talk). Right: the cell and its image under ×α.

    frobmat: the integer matrix of α on the ordered basis (1, τ) — rows, with
    columns the images of the basis vectors (ECQFIsogenyClass frobmat .vec).
    It must be passed in because qf may sit at any level of the volcano: the
    (A, B, C)-only formula works just for the Frobenius-order form itself.
    Points are drawn only while N stays below the picture guard."""
    A, B, C = (int(v) for v in qf)
    (m00, m01), (m10, m11) = ((int(v) for v in row) for row in frobmat)
    payload = {"A": A, "B": B, "C": C, "a": a, "p": p,
               "M": [[m00, m01], [m10, m11]],
               "drawPts": bool(n_points), "H": height_px}
    return '<meta charset="utf-8">' + _HEAD + """
<div class="panel" style="max-width:860px;">
  <div class="stage">
    <div class="cell">
      <div class="cap">E(𝔽_p) on the torus ℂ/Λ — the fixed points of the Frobenius lift α</div>
      <canvas id="ctA" width="400" height="__H__"></canvas>
    </div>
    <div class="cell">
      <div class="cap">multiplication by α on Λ = ⟨1, τ⟩ — the blue cell and its ×α image</div>
      <canvas id="ctB" width="400" height="__H__"></canvas>
    </div>
  </div>
  <div class="info" id="ctInfo"></div>
</div>
<script>
(() => {
"use strict";
__NAV_JS__
const DATA = __DATA__;
const ACCENT="#4da3d8", GOLD="#e0b64f", RED="#ef6f6f", INK="#d7d9dc", MUT="#9aa4ad";
const mod=(x,n)=>((x%n)+n)%n;
const {A,B,C,a,p}=DATA;
const D=B*B-4*A*C, N=p+1-a;
const tx=-B/(2*A), ty=Math.sqrt(-D)/(2*A);

/* ── left: Fix(α) on the fundamental cell (talk's CM torus) ─────────── */
(function(){
  const cv=document.getElementById("ctA"), ctx=cv.getContext("2d");
  const pad=26;
  const scale=Math.min((cv.width-2*pad)/(1+Math.abs(tx)),(cv.height-2*pad)/ty);
  const ox=pad+(tx<0?-tx*scale:0), oy=cv.height-pad;
  const X=(u,v)=>ox+(u+v*tx)*scale, Y=(u,v)=>oy-v*ty*scale;
  ctx.strokeStyle="#555"; ctx.lineWidth=2;
  ctx.beginPath(); ctx.moveTo(X(0,0),Y(0,0)); ctx.lineTo(X(1,0),Y(1,0));
  ctx.lineTo(X(1,1),Y(1,1)); ctx.lineTo(X(0,1),Y(0,1)); ctx.closePath(); ctx.stroke();
  if(DATA.drawPts){
    // mult by (α − 1) on basis (1, τ): K = M − I (columns are the images);
    // the fixed points are generated by the columns of adj(K) mod N
    const M=DATA.M;
    const s=M[0][0]-1, t=M[1][0], u=M[0][1], v=M[1][1]-1;
    const g1x=mod(v,N), g1y=mod(-t,N), g2x=mod(-u,N), g2y=mod(s,N);
    const seen=new Set(), pts=[];
    outer:
    for(let i=0;i<N;i++) for(let j=0;j<N;j++){
      const x=mod(i*g1x+j*g2x,N), y=mod(i*g1y+j*g2y,N);
      const key=x*N+y;
      if(!seen.has(key)){ seen.add(key); pts.push([x,y]); if(pts.length===N) break outer; }
    }
    pts.forEach(([x,y],k)=>{
      ctx.fillStyle=`hsl(${205+130*k/N},65%,${58-18*k/N}%)`;
      ctx.beginPath(); ctx.arc(X(x/N,y/N),Y(x/N,y/N),Math.max(2.2,8-N/40),0,7); ctx.fill();
    });
  } else {
    ctx.fillStyle=MUT; ctx.font="12px system-ui"; ctx.textAlign="center";
    ctx.fillText("N = "+N+" points — beyond the drawing guard", cv.width/2, cv.height/2);
  }
})();

/* ── right: the ×α cell (adapted from the §2 lift applet) ───────────── */
(function(){
  const cv=document.getElementById("ctB"), ctx=cv.getContext("2d");
  // α on ⟨1, τ⟩ via its integer matrix: the cell corners land exactly on Λ
  const M=DATA.M, tau={x:tx,y:ty};
  const al={x:M[0][0]+M[1][0]*tx, y:M[1][0]*ty};               // α·1
  const at={x:M[0][1]+M[1][1]*tx, y:M[1][1]*ty};               // α·τ
  const keys=[{x:0,y:0},{x:1,y:0},tau,{x:1+tau.x,y:tau.y},al,at,{x:al.x+at.x,y:al.y+at.y}];
  const pad=0.7;
  let x0=Math.min(...keys.map(k=>k.x))-pad, x1=Math.max(...keys.map(k=>k.x))+pad,
      y0=Math.min(...keys.map(k=>k.y))-pad, y1=Math.max(...keys.map(k=>k.y))+pad;
  const need=(x1-x0)/(y1-y0), have=cv.width/cv.height;
  if(need<have){const c=(x0+x1)/2,w=(y1-y0)*have; x0=c-w/2; x1=c+w/2;}
  else{const c=(y0+y1)/2,h=(x1-x0)/have; y0=c-h/2; y1=c+h/2;}
  const PX=x=>(x-x0)/(x1-x0)*cv.width, PY=y=>(y1-y)/(y1-y0)*cv.height;
  // lattice m + n·τ
  ctx.fillStyle="rgba(255,255,255,0.3)";
  for(let n=Math.floor(y0/ty)-1;n<=Math.ceil(y1/ty)+1;n++)
    for(let m=Math.floor(x0-n*tx)-1;m<=Math.ceil(x1-n*tx)+1;m++){
      const x=m+n*tx, y=n*ty;
      if(x>x0&&x<x1&&y>y0&&y<y1){ctx.beginPath();ctx.arc(PX(x),PY(y),2,0,7);ctx.fill();}
    }
  const poly=(vs,fill,stroke,dash)=>{
    ctx.beginPath();
    vs.forEach((v,i)=>{i?ctx.lineTo(PX(v.x),PY(v.y)):ctx.moveTo(PX(v.x),PY(v.y));});
    ctx.closePath();
    if(dash)ctx.setLineDash([5,4]);
    ctx.fillStyle=fill; ctx.fill(); ctx.strokeStyle=stroke; ctx.lineWidth=1.6; ctx.stroke();
    ctx.setLineDash([]);
  };
  poly([{x:0,y:0},{x:1,y:0},{x:1+tau.x,y:tau.y},tau],"rgba(77,163,216,0.16)",ACCENT,false);
  poly([{x:0,y:0},al,{x:al.x+at.x,y:al.y+at.y},at],"rgba(224,182,79,0.12)",GOLD,true);
  ctx.font="12px system-ui";
  ctx.fillStyle=INK; ctx.fillText("1",PX(1)+4,PY(0)+13);
  ctx.fillStyle=ACCENT; ctx.fillText("τ",PX(tau.x)+5,PY(tau.y)-5);
  ctx.fillStyle=GOLD; ctx.fillText("α",PX(al.x)+5,PY(al.y)-5);
  ctx.fillText("ατ",PX(at.x)+5,PY(at.y)-5);
})();

document.getElementById("ctInfo").innerHTML=
  "χ(x) = x² − "+a+"x + "+p+" · #E(𝔽_p) = χ(1) = "+N
  +" · τ ≈ "+tx.toFixed(3)+" + "+ty.toFixed(3)+"i · |α| = √"+p;
})();
</script>
""".replace("__DATA__", json.dumps(payload)).replace("__NAV_JS__", _NAV_JS).replace("__H__", str(height_px))
