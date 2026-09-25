// three.js viewer for a gridded surface.
//
// The model hangs in a pivot whose orientation is a quaternion and the camera sits on the +z axis
// looking at the origin. Dragging premultiplies a small rotation about the screen axis perpendicular
// to the drag (no gimbal lock); shift+drag rolls about the view axis. Wheel or +/- zooms, WASDQE
// translates the model (A/D left/right, W/S up/down, Q/E toward/away), R resets, F re-frames.

import * as THREE from '../vendor/three.module.js';

const VERT = /* glsl */ `
attribute vec2 st;
attribute float val;
varying vec3 vN;
varying vec3 vPos;
varying vec2 vSt;
varying float vVal;
void main() {
  vec4 mv = modelViewMatrix * vec4(position, 1.0);
  vPos = mv.xyz;
  vN = normalMatrix * normal;
  vSt = st;
  vVal = val;
  gl_Position = projectionMatrix * mv;
}`;

const FRAG = /* glsl */ `
uniform vec3 uFront;
uniform vec3 uBack;
uniform vec3 uInk;
uniform int uMode;        // 0 two-sided, 1 sequential map of val, 2 diverging map of val
uniform vec2 uRange;      // sequential: [lo, hi] -> [0,1]; diverging: [lo, 0] -> [-1,0], [0, hi] -> [0,1]
uniform float uGrid;      // 0 = off
uniform float uGridStep;
uniform float uLineW;
uniform float uWire;      // 1 = wireframe: lines take the surface's colour
uniform float uOpacity;   // of the fill
uniform int uPass;        // 0 fill and lines, 1 lines only, 2 fill only
varying vec3 vN;
varying vec3 vPos;
varying vec2 vSt;
varying float vVal;

// viridis, polynomial fit (after mattz, CC0)
vec3 viridis(float t) {
  const vec3 c0 = vec3(0.2777, 0.0054, 0.3341);
  const vec3 c1 = vec3(0.1051, 1.4046, 1.3846);
  const vec3 c2 = vec3(-0.3309, 0.2148, 0.0951);
  const vec3 c3 = vec3(-4.6342, -5.7991, -19.3324);
  const vec3 c4 = vec3(6.2283, 14.1799, 56.6906);
  const vec3 c5 = vec3(4.7764, -13.7451, -65.3530);
  const vec3 c6 = vec3(-5.4355, 4.6459, 26.3124);
  return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}
vec3 diverging(float t) { // t in [-1, 1]: blue - paper - red
  vec3 mid = vec3(0.93, 0.91, 0.86);
  vec3 neg = vec3(0.18, 0.36, 0.62);
  vec3 pos = vec3(0.72, 0.24, 0.16);
  return t < 0.0 ? mix(mid, neg, -t) : mix(mid, pos, t);
}

void main() {
  vec3 n = normalize(vN);
  vec3 v = normalize(-vPos);
  bool nSide = dot(n, v) >= 0.0;   // are we looking at the side N points to?
  if (!nSide) n = -n;

  vec3 base = nSide ? uFront : uBack;
  if (uMode == 1) {
    base = viridis(clamp((vVal - uRange.x) / max(uRange.y - uRange.x, 1e-9), 0.0, 1.0));
    if (!nSide) base *= 0.8;
  } else if (uMode == 2) {
    // each side scaled separately: for K the positive side is bounded by 1/4, the negative side is not
    float t = vVal < 0.0 ? vVal / max(abs(uRange.x), 1e-9) : vVal / max(abs(uRange.y), 1e-9);
    base = diverging(clamp(t, -1.0, 1.0));
    if (!nSide) base *= 0.8;
  }

  // camera-fixed lights: key, fill, rim
  vec3 L1 = normalize(vec3(-0.45, 0.65, 0.62));
  vec3 L2 = normalize(vec3(0.7, -0.25, 0.45));
  float dif = 0.30 + 0.62 * max(dot(n, L1), 0.0) + 0.22 * max(dot(n, L2), 0.0);
  float rim = pow(1.0 - max(dot(n, v), 0.0), 3.0) * 0.12;
  float spec = pow(max(dot(n, normalize(L1 + v)), 0.0), 48.0) * 0.28;
  vec3 col = base * dif + vec3(spec + rim);

  float line = 0.0;
  if (uGrid > 0.5 && uPass != 2) {
    vec2 g = vSt / uGridStep;
    vec2 w = fwidth(g);
    vec2 d = abs(fract(g - 0.5) - 0.5) / max(w, vec2(1e-6));
    float dmin = min(d.x, d.y);
    // fade lines out where they would be denser than a few pixels
    float dens = clamp(1.5 - 3.0 * max(w.x, w.y), 0.0, 1.0);
    line = (1.0 - smoothstep(uLineW - 0.5, uLineW + 0.5, dmin)) * dens;
  }
  if (uPass == 1) {
    // lines alone, drawn first and writing depth, so the fill in front of them tints them
    if (line < 0.02) discard;
    gl_FragColor = vec4(uWire > 0.5 ? mix(col, uInk, 0.35) : uInk, min(1.0, line * (uWire > 0.5 ? 1.0 : 0.8)));
  } else if (uPass == 2) {
    gl_FragColor = vec4(col, uOpacity);
  } else {
    gl_FragColor = vec4(mix(col, uInk, 0.8 * line), 1.0);
  }
}`;

export class Viewer {
  constructor(canvas) {
    this.canvas = canvas;
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(32, 1, 0.01, 1000);
    this.pivot = new THREE.Group();
    this.holder = new THREE.Group(); // translated by -centre so the pivot rotates about the model centre
    this.scene.add(this.pivot);
    this.pivot.add(this.holder);

    this.uniforms = {
      uFront: { value: new THREE.Color('#6f8fb0') },
      uBack: { value: new THREE.Color('#e8c9a0') },
      uInk: { value: new THREE.Color('#2b2620') },
      uMode: { value: 0 },
      uRange: { value: new THREE.Vector2(0, 1) },
      uGrid: { value: 1 },
      uGridStep: { value: 0.25 },
      uLineW: { value: 0.6 },
      uWire: { value: 0 },
      uOpacity: { value: 1 },
    };
    this.geometry = new THREE.BufferGeometry();
    // opaque: one pass. See-through or wireframe: the lines (opaque, writing depth), then the fill
    // (blended, not writing depth), so the lines never sort wrongly against the fill.
    const mesh = (pass, opts, order) => {
      const material = new THREE.ShaderMaterial({
        vertexShader: VERT, fragmentShader: FRAG, side: THREE.DoubleSide,
        uniforms: { ...this.uniforms, uPass: { value: pass } }, ...opts,
      });
      const m = new THREE.Mesh(this.geometry, material);
      m.frustumCulled = false;
      m.renderOrder = order;
      this.holder.add(m);
      return m;
    };
    this.mesh = mesh(0, {}, 0);
    this.lines = mesh(1, { transparent: true }, 1);
    this.fill = mesh(2, { transparent: true, depthWrite: false }, 2);
    this.wire = false;
    this.opacity = 1;
    this._passes();

    this.dist = 10;           // camera distance
    this.radius = 3;          // model radius from the last framing
    this.pan = new THREE.Vector3();
    this.centre = new THREE.Vector3(); // model-space point at the pivot (set by frame)
    this.scale = 1;                    // homothety applied to the model (mean curvature H = 1/(2 scale))
    this.keys = new Set();
    this.dirty = true;
    this.gridDims = null;
    this.resetOrientation();

    this._bindPointer();
    this._bindKeys();
    new ResizeObserver(() => this.resize()).observe(canvas.parentElement);
    this.resize();
    let last = performance.now();
    const loop = (now) => {
      const dt = Math.min(0.1, (now - last) / 1000);
      last = now;
      this._applyKeys(dt);
      if (this.dirty) this.render();
      requestAnimationFrame(loop);
    };
    requestAnimationFrame(loop);
  }

  // -------------------------------------------------------------- geometry

  /** Show a surface from computeSurface. colourVal is a per-vertex Float32Array or null. */
  setSurface(res, colourVal) {
    const { nx, ny, pos, nrm, st } = res;
    const g = this.geometry;
    const key = `${nx}x${ny}`;
    const finite = (v) => Number.isFinite(pos[3 * v]) && Number.isFinite(pos[3 * v + 1]) && Number.isFinite(pos[3 * v + 2]);
    // rebuild the index when the grid or the set of bad vertices changes
    let bad = false;
    for (let v = 0; v < nx * ny && !bad; v++) bad = !finite(v);
    if (this.gridDims !== key || bad || this.hadBad) {
      const idx = [];
      for (let j = 0; j + 1 < ny; j++) {
        for (let i = 0; i + 1 < nx; i++) {
          const a = j * nx + i, b = a + 1, c = a + nx, d = c + 1;
          if (bad && !(finite(a) && finite(b) && finite(c) && finite(d))) continue;
          idx.push(a, b, d, a, d, c);
        }
      }
      g.setIndex(nx * ny > 65535 ? new THREE.Uint32BufferAttribute(idx, 1) : new THREE.Uint16BufferAttribute(idx, 1));
      this.gridDims = bad ? null : key;
      this.hadBad = bad;
    }
    if (bad) for (let e = 0; e < pos.length; e++) if (!Number.isFinite(pos[e])) pos[e] = 0;
    g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    g.setAttribute('normal', new THREE.BufferAttribute(nrm, 3));
    g.setAttribute('st', new THREE.BufferAttribute(st, 2));
    g.setAttribute('val', new THREE.BufferAttribute(colourVal || new Float32Array(nx * ny), 1));
    this.bbox = res.bbox;
    this.dirty = true;
  }

  setColourValues(vals) {
    this.geometry.setAttribute('val', new THREE.BufferAttribute(vals, 1));
    this.dirty = true;
  }

  setStyle({ mode, range, grid, gridStep, lineWidth, front, back, ink, wire, opacity }) {
    const u = this.uniforms;
    if (wire !== undefined) this.wire = wire;
    if (opacity !== undefined) this.opacity = opacity;
    if (mode !== undefined) u.uMode.value = mode;
    if (range) u.uRange.value.set(range[0], range[1]);
    if (grid !== undefined) u.uGrid.value = grid ? 1 : 0;
    if (gridStep !== undefined) u.uGridStep.value = gridStep;
    if (lineWidth !== undefined) u.uLineW.value = lineWidth;
    if (front) u.uFront.value.set(front);
    if (back) u.uBack.value.set(back);
    if (ink) u.uInk.value.set(ink);
    this._passes();
    this.dirty = true;
  }

  _passes() {
    const u = this.uniforms, see = this.opacity < 1;
    u.uWire.value = this.wire ? 1 : 0;
    u.uOpacity.value = this.opacity;
    this.mesh.visible = !this.wire && !see;
    this.lines.visible = this.wire || (see && u.uGrid.value > 0.5);
    this.fill.visible = !this.wire && see && this.opacity > 0;
  }

  setBackground(colour) {
    this.renderer.setClearColor(new THREE.Color(colour), 1);
    this.dirty = true;
  }

  // -------------------------------------------------------------- camera

  resetOrientation() {
    // a three-quarter view
    this.pivot.quaternion.setFromEuler(new THREE.Euler(0.45, -0.6, 0, 'XYZ'));
    this.pan.set(0, 0, 0);
    this.dirty = true;
  }

  /** Centre the model at the pivot and fit it in view. */
  frame() {
    const b = this.bbox;
    if (!b || !Number.isFinite(b[0])) return;
    this.centre.set((b[0] + b[3]) / 2, (b[1] + b[4]) / 2, (b[2] + b[5]) / 2);
    this.radius = Math.max(1e-3, 0.5 * this.scale * Math.hypot(b[3] - b[0], b[4] - b[1], b[5] - b[2]));
    this.holder.position.copy(this.centre).multiplyScalar(-this.scale);
    const fov = THREE.MathUtils.degToRad(this.camera.fov);
    const aspect = Math.min(1, this.camera.aspect);
    this.dist = (1.05 * this.radius) / Math.sin(fov / 2) / Math.sqrt(aspect);
    this.pan.set(0, 0, 0);
    this.dirty = true;
  }

  /** Scale the model about the pivot, keeping the camera where it is. */
  setScale(s) {
    this.scale = s;
    this.holder.scale.setScalar(s);
    this.holder.position.copy(this.centre).multiplyScalar(-s);
    this.dirty = true;
  }

  zoom(factor) {
    this.dist = THREE.MathUtils.clamp(this.dist * factor, this.radius * 0.02, this.radius * 200);
    this.dirty = true;
  }

  resize() {
    const p = this.canvas.parentElement;
    const w = p.clientWidth, h = p.clientHeight;
    if (!w || !h) return;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.dirty = true;
  }

  render() {
    const cam = this.camera;
    cam.position.set(0, 0, this.dist);
    // bounding sphere of the current model about the pivot (the surface may have grown since framing)
    let r = this.radius;
    const b = this.bbox, o = this.holder.position, k = this.scale;
    if (b && Number.isFinite(b[0])) {
      const ex = Math.max(Math.abs(k * b[0] + o.x), Math.abs(k * b[3] + o.x));
      const ey = Math.max(Math.abs(k * b[1] + o.y), Math.abs(k * b[4] + o.y));
      const ez = Math.max(Math.abs(k * b[2] + o.z), Math.abs(k * b[5] + o.z));
      r = Math.hypot(ex, ey, ez);
    }
    const d = Math.hypot(this.pan.x, this.pan.y, this.dist - this.pan.z);
    cam.near = Math.max(d - 1.05 * r, d * 1e-3);
    cam.far = d + 1.05 * r;
    cam.updateProjectionMatrix();
    this.pivot.position.copy(this.pan);
    this.renderer.render(this.scene, cam);
    this.dirty = false;
  }

  /** PNG data URL of the current view. */
  snapshot() {
    this.render();
    return this.canvas.toDataURL('image/png');
  }

  /** Wavefront OBJ of the current mesh (positions, normals, grid coordinates as texture coordinates). */
  toOBJ() {
    const g = this.geometry;
    const P = g.getAttribute('position'), N = g.getAttribute('normal'), T = g.getAttribute('st');
    const I = g.getIndex();
    if (!P || !I) return '';
    const out = ['# CMC immersion, exported from the spectral-data viewer'];
    const f = (x) => (Math.abs(x) < 1e-12 ? '0' : x.toPrecision(7));
    const k = this.scale;
    for (let v = 0; v < P.count; v++) out.push(`v ${f(k * P.getX(v))} ${f(k * P.getY(v))} ${f(k * P.getZ(v))}`);
    for (let v = 0; v < N.count; v++) out.push(`vn ${f(N.getX(v))} ${f(N.getY(v))} ${f(N.getZ(v))}`);
    for (let v = 0; v < T.count; v++) out.push(`vt ${f(T.getX(v))} ${f(T.getY(v))}`);
    for (let k = 0; k < I.count; k += 3) {
      const a = I.getX(k) + 1, b = I.getX(k + 1) + 1, c = I.getX(k + 2) + 1;
      out.push(`f ${a}/${a}/${a} ${b}/${b}/${b} ${c}/${c}/${c}`);
    }
    return out.join('\n') + '\n';
  }

  // -------------------------------------------------------------- input

  _rotateScreen(dx, dy) {
    const len = Math.hypot(dx, dy);
    if (!len) return;
    const angle = (len / this.canvas.clientHeight) * Math.PI * 1.2;
    const q = new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(dy / len, dx / len, 0), angle);
    this.pivot.quaternion.premultiply(q);
    this.dirty = true;
  }

  _roll(angle) {
    const q = new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0, 0, 1), angle);
    this.pivot.quaternion.premultiply(q);
    this.dirty = true;
  }

  _bindPointer() {
    const c = this.canvas;
    const pts = new Map();
    const centre = () => {
      const r = c.getBoundingClientRect();
      return [r.left + r.width / 2, r.top + r.height / 2];
    };
    let pinch = null;
    c.addEventListener('pointerdown', (e) => {
      c.setPointerCapture(e.pointerId);
      pts.set(e.pointerId, { x: e.clientX, y: e.clientY });
      if (pts.size === 2) {
        const [a, b] = [...pts.values()];
        pinch = { d: Math.hypot(a.x - b.x, a.y - b.y), ang: Math.atan2(b.y - a.y, b.x - a.x) };
      }
    });
    c.addEventListener('pointermove', (e) => {
      const p = pts.get(e.pointerId);
      if (!p) return;
      const nx = e.clientX, ny = e.clientY;
      if (pts.size === 1) {
        if (e.shiftKey) {
          const [cx, cy] = centre();
          const a0 = Math.atan2(p.y - cy, p.x - cx), a1 = Math.atan2(ny - cy, nx - cx);
          let d = a1 - a0;
          if (d > Math.PI) d -= 2 * Math.PI;
          if (d < -Math.PI) d += 2 * Math.PI;
          this._roll(-d);
        } else {
          this._rotateScreen(nx - p.x, ny - p.y);
        }
      }
      p.x = nx; p.y = ny;
      if (pts.size === 2 && pinch) {
        const [a, b] = [...pts.values()];
        const d = Math.hypot(a.x - b.x, a.y - b.y), ang = Math.atan2(b.y - a.y, b.x - a.x);
        if (d > 0) this.zoom(pinch.d / d);
        let da = ang - pinch.ang;
        if (da > Math.PI) da -= 2 * Math.PI;
        if (da < -Math.PI) da += 2 * Math.PI;
        this._roll(-da);
        pinch = { d, ang };
      }
    });
    const up = (e) => { pts.delete(e.pointerId); if (pts.size < 2) pinch = null; };
    c.addEventListener('pointerup', up);
    c.addEventListener('pointercancel', up);
    c.addEventListener('wheel', (e) => {
      e.preventDefault();
      const dy = e.deltaMode === 1 ? e.deltaY * 16 : e.deltaY;
      this.zoom(Math.exp(dy * 0.0012));
    }, { passive: false });
  }

  _bindKeys() {
    const typing = (e) => {
      const t = e.target;
      return t && (t.isContentEditable || (t.tagName === 'INPUT' && t.type !== 'range' && t.type !== 'checkbox')
        || t.tagName === 'TEXTAREA' || t.tagName === 'SELECT');
    };
    window.addEventListener('keydown', (e) => {
      if (typing(e) || e.ctrlKey || e.metaKey || e.altKey) return;
      const k = e.key.toLowerCase();
      if ('wasdqe'.includes(k) && k.length === 1) { this.keys.add(k); e.preventDefault(); }
      else if (k === '+' || k === '=') { this.zoom(1 / 1.15); e.preventDefault(); }
      else if (k === '-' || k === '_') { this.zoom(1.15); e.preventDefault(); }
      else if (k === 'r') { this.resetOrientation(); this.frame(); }
      else if (k === 'f') { this.frame(); }
    });
    window.addEventListener('keyup', (e) => this.keys.delete(e.key.toLowerCase()));
    window.addEventListener('blur', () => this.keys.clear());
  }

  _applyKeys(dt) {
    if (!this.keys.size) return;
    const v = 0.5 * this.dist * dt;
    const k = this.keys;
    if (k.has('a')) this.pan.x -= v;
    if (k.has('d')) this.pan.x += v;
    if (k.has('w')) this.pan.y += v;
    if (k.has('s')) this.pan.y -= v;
    if (k.has('q')) this.pan.z += v;
    if (k.has('e')) this.pan.z -= v;
    this.dirty = true;
  }
}
