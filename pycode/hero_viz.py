"""Homepage hero: the Frobenius flow as a glowing galaxy banner.

A full-width, fixed-height canvas that autoplays a slow, seamless loop of the
Frobenius flow z -> alpha^phi z on the CM torus, drawn in the multiplicative
("galaxy") model with additive glow.  Pinned to y^2 = x^3 + 3x over F_5,
alpha = -2 + i; the points are E(F_5^level) = ker(alpha^level - 1) (level 2 =
F_25, the minimal 20-point constellation).  Colour is by angle, so the whole
picture returns to itself after one Frobenius beat -> a clean loop.
"""


def frobenius_hero_html(level: int = 2, speed: float = 0.05, height: int = 380) -> str:
    return r"""
<style>
  html,body{margin:0;background:#0e1013;}
  #hero{display:block;width:100%;max-width:760px;height:HEROHpx;margin:0 auto;background:#0e1013;}
</style>
<canvas id="hero"></canvas>
<script>
"use strict";
const cv=document.getElementById("hero"), ctx=cv.getContext("2d");
const dpr=Math.min(window.devicePixelRatio||1, 2);
const ALPHA=[-2,1], THETA=Math.atan2(1,-2), RHO=Math.sqrt(5);
const TAUP=[0.4,0.2], DEN=[-2,1], LEVEL=__LEVEL__, SPEED=__SPEED__;
const cmul=(a,b)=>[a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]];
const apow=p=>{const r=Math.pow(RHO,p),t=THETA*p;return [r*Math.cos(t), r*Math.sin(t)];};
const apowInt=n=>{let r=[1,0];for(let i=0;i<n;i++)r=cmul(r,ALPHA);return r;};
const r01=x=>x-Math.floor(x);

// E(F_5^level) = ker(alpha^level - 1), enumerated by box + dedup
function genField(m){
  const a=apowInt(m), b=[a[0]-1,a[1]], N=b[0]*b[0]+b[1]*b[1];
  const cr=b[0], ci=-b[1], B=Math.ceil(Math.sqrt(N))+1, seen=new Set(), pts=[];
  for(let x=-B;x<=B;x++)for(let y=-B;y<=B;y++){
    const u=r01((x*cr - y*ci)/N), v=r01((x*ci + y*cr)/N);
    const key=Math.round(u*N)+","+Math.round(v*N);
    if(!seen.has(key)){seen.add(key); pts.push([u,v]);}
  }
  return pts;
}
const G=genField(LEVEL);

function toGalaxy(q, CX, CY, R){
  const [zr,zi]=q,[dr,di]=DEN,n2=dr*dr+di*di;
  const wr=(zr*dr+zi*di)/n2, wi=(zi*dr-zr*di)/n2;
  const vs=wi/TAUP[1], vp=r01(vs), up=r01(wr-vs*TAUP[0]);
  const rr=R*Math.exp(-2*Math.PI*vp*TAUP[1]), th=2*Math.PI*(up+vp*TAUP[0]);
  return [CX+rr*Math.cos(th), CY-rr*Math.sin(th), th];
}
function hue(a){
  const h=r01(a/(2*Math.PI))*6, x=1-Math.abs(h%2-1);
  const c=[[1,x,0],[x,1,0],[0,1,x],[0,x,1],[x,0,1],[1,0,x]][Math.floor(h)%6];
  return c.map(v=>Math.round((0.28+0.72*v)*255));
}
const halo=(x,y,rad,r,g,b,a)=>{
  const grd=ctx.createRadialGradient(x,y,0,x,y,rad);
  grd.addColorStop(0,`rgba(${r},${g},${b},${a})`); grd.addColorStop(1,`rgba(${r},${g},${b},0)`);
  ctx.fillStyle=grd; ctx.beginPath(); ctx.arc(x,y,rad,0,7); ctx.fill();
};

let t0=null;
function frame(ms){
  const w=Math.round(cv.clientWidth*dpr), h=Math.round(cv.clientHeight*dpr);
  if(cv.width!==w||cv.height!==h){cv.width=w; cv.height=h;}
  const CX=cv.width/2, CY=cv.height/2, R=cv.height/2 - 26*dpr;
  if(t0===null) t0=ms;
  const phi=r01((ms-t0)*0.001*SPEED), af=apow(phi);
  ctx.globalCompositeOperation="source-over"; ctx.fillStyle="#0e1013"; ctx.fillRect(0,0,cv.width,cv.height);
  ctx.globalCompositeOperation="lighter";
  for(const g of G){
    for(let j=4;j>=1;j--){
      const s=toGalaxy(cmul(apow(phi-j*0.011),g),CX,CY,R), [r,gg,b]=hue(s[2]), f=1-j/5;
      halo(s[0],s[1], (5+7*f)*dpr, r,gg,b, 0.16*f);
    }
    const s=toGalaxy(cmul(af,g),CX,CY,R), [r,gg,b]=hue(s[2]);
    halo(s[0],s[1], 34*dpr, r,gg,b, 0.20);
    halo(s[0],s[1], 12*dpr, r,gg,b, 0.55);
    ctx.fillStyle=`rgba(${r},${gg},${b},0.95)`; ctx.beginPath(); ctx.arc(s[0],s[1],3.0*dpr,0,7); ctx.fill();
  }
  requestAnimationFrame(frame);
}
requestAnimationFrame(frame);
</script>
""".replace("HEROH", str(height)).replace("__LEVEL__", str(level)).replace("__SPEED__", repr(speed))
