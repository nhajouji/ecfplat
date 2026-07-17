"""Interactive 3-D Pinkall-torus viewer for the Background page.

Loads a precomputed geometry JSON (exported from Steve Trettel's elliptic-curve-viz
via its exact Hopf/Pinkall pipeline — points of E(F_{p^k}) on the flat torus in S³,
stereographically projected to R³ + the torus surface mesh) and renders it live with
three.js: orbitable, auto-rotating, points coloured by field degree.  Self-contained
(three.js from CDN, geometry inlined) so it drops into Streamlit `components.html`.

Regenerate the JSON with elliptic-curve-viz/scripts/export-geometry.ts (see
`ecfplat-elliptic-curve-viz-survey` notes)."""
import json
from pathlib import Path


def hopf_torus_html(json_name: str = "hopf_disc8_k2.json", height: int = 460) -> str:
    data = json.loads((Path(__file__).parent / "data" / json_name).read_text())
    return _HTML.replace("__DATA__", json.dumps(data)).replace("__H__", str(height))


_HTML = r"""
<style>
  html,body{margin:0;background:#0e1013;}
  #hopf{display:block;width:100%;height:__H__px;background:#0e1013;border-radius:12px;touch-action:none;}
</style>
<canvas id="hopf"></canvas>
<script type="importmap">
{ "imports": {
    "three": "https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js",
    "three/addons/": "https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/"
} }
</script>
<script type="module">
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const DATA = __DATA__;
const cv = document.getElementById('hopf');
const renderer = new THREE.WebGLRenderer({ canvas: cv, antialias: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x0e1013);
const camera = new THREE.PerspectiveCamera(45, 2, 0.05, 200);

const controls = new OrbitControls(camera, cv);
controls.enableDamping = true; controls.dampingFactor = 0.08;
controls.autoRotate = true; controls.autoRotateSpeed = 0.55;
controls.enablePan = false; controls.minDistance = 3; controls.maxDistance = 40;

// paper's (x, z, -y) swizzle: stand the torus up
const group = new THREE.Group();
group.rotation.x = -Math.PI / 2;
scene.add(group);

// ── the Pinkall torus surface ────────────────────────────────────────────────
const S = DATA.surface, nu = S.nu, nv = S.nv;
const pos = new Float32Array(S.verts.length * 3);
S.verts.forEach((v, i) => { pos[3*i] = v[0]; pos[3*i+1] = v[1]; pos[3*i+2] = v[2]; });
const idx = [];
for (let iu = 0; iu < nu - 1; iu++) for (let iv = 0; iv < nv - 1; iv++) {
  const a = iu*nv + iv, b = (iu+1)*nv + iv, c = (iu+1)*nv + (iv+1), d = iu*nv + (iv+1);
  idx.push(a, b, d, b, c, d);
}
const geo = new THREE.BufferGeometry();
geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
geo.setIndex(idx);
geo.computeVertexNormals();
const surf = new THREE.Mesh(geo, new THREE.MeshStandardMaterial({
  color: 0xcfd4da, roughness: 0.62, metalness: 0.0, side: THREE.DoubleSide,
  transparent: true, opacity: 0.9 }));
group.add(surf);

// ── the E(F_{p^k}) points, coloured by field degree ──────────────────────────
const N = DATA.points.length;
const inst = new THREE.InstancedMesh(
  new THREE.SphereGeometry(1, 18, 14),
  new THREE.MeshStandardMaterial({ roughness: 0.3, metalness: 0.15 }), N);
const dummy = new THREE.Object3D();
const PAL = { 1: new THREE.Color(0xef6f6f), 2: new THREE.Color(0xe0b64f),
              3: new THREE.Color(0x69b382), 6: new THREE.Color(0x8f7ae0) };
for (let i = 0; i < N; i++) {
  const p = DATA.points[i], d = DATA.degrees[i];
  dummy.position.set(p[0], p[1], p[2]);
  dummy.scale.setScalar(d === 1 ? 0.15 : 0.10);   // F_p points a bit bigger
  dummy.updateMatrix();
  inst.setMatrixAt(i, dummy.matrix);
  inst.setColorAt(i, PAL[d] || PAL[6]);
}
inst.instanceMatrix.needsUpdate = true;
if (inst.instanceColor) inst.instanceColor.needsUpdate = true;
group.add(inst);

// ── lights ───────────────────────────────────────────────────────────────────
scene.add(new THREE.AmbientLight(0xffffff, 0.55));
const key = new THREE.DirectionalLight(0xffffff, 1.0); key.position.set(4, 6, 5); scene.add(key);
const rim = new THREE.DirectionalLight(0x88aaff, 0.4); rim.position.set(-5, -3, -4); scene.add(rim);

// ── frame the camera on the geometry ─────────────────────────────────────────
const box = new THREE.Box3().setFromObject(group);
const center = box.getCenter(new THREE.Vector3());
const radius = box.getSize(new THREE.Vector3()).length() / 2;
controls.target.copy(center);
camera.position.set(center.x + radius*1.5, center.y + radius*0.95, center.z + radius*1.5);
camera.updateProjectionMatrix();

function resize() {
  const w = cv.clientWidth, h = cv.clientHeight;
  if (cv.width !== Math.round(w*renderer.getPixelRatio()) ||
      cv.height !== Math.round(h*renderer.getPixelRatio())) {
    renderer.setSize(w, h, false); camera.aspect = w / h; camera.updateProjectionMatrix();
  }
}
function loop() { resize(); controls.update(); renderer.render(scene, camera); requestAnimationFrame(loop); }
loop();
</script>
"""
