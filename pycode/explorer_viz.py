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
select. When ``link_base`` is given, the selected node's readout offers an
``<a target="_parent">`` link with ``&node={i}`` appended — this is the
query-param navigation used by the drill-down Explorer.

Embedded via ``st.components.v1.html`` like the Background applets; styling
matches the shared panel look (``modular_viz._HEAD``).
"""

import json

from modular_viz import _HEAD


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
    s+=' &nbsp;<a target="_parent" style="color:'+ACCENT+';" href="'+DATA.linkBase+"&node="+i+'">'+DATA.linkLabel+"</a>";
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
""".replace("__DATA__", json.dumps(payload)).replace("__H__", str(height_px))
