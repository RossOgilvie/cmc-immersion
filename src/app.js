// Glue: state, panel controls, progressive recomputation in a worker, presets, URL state, export.

import { Viewer } from './viewer.js';
import { SpectralWidget, MAX_GENUS } from './spectral-widget.js';
import { shapeFlows, kappa0, hopfArg } from './cmc/cmc.js';
import { WhithamCurve } from './cmc/whitham.js';

const $ = (id) => document.getElementById(id);
const polar = (r, t) => [r * Math.cos(t), r * Math.sin(t)];
const PI = Math.PI;

// ------------------------------------------------------------------ state

const VIEW_DEFAULTS = {
  tau: [], z0: [0, 0], width: 12, height: 12, curv: true, phi: 0,
  gridStep: 0.25, lineW: 0.6, hmax: 0.05, res: 240,
};
// Widths of the rotational examples are exact closing periods of the parallels, pi / sqrt(det X(lam0)).
// Bubbletons: a double branch point at the k-th resonance point of the cylinder, where its monodromy is ±I.
const RES = (k) => (k - Math.sqrt(k * k - 1)) ** 2;
const PRESETS = [
  { name: 'Round cylinder (g = 0)', s: { alphas: [], theta0: 0, width: PI, height: 10, gridStep: PI / 16 } },
  { name: 'Unduloid (g = 1)', s: { alphas: [polar(0.45, 0)], theta0: PI, width: 2.90682, height: 14 } },
  { name: 'Unduloid near a sphere chain', s: { alphas: [polar(0.08, 0)], theta0: PI, width: 1.645512, height: 14, gridStep: 0.1 } },
  { name: 'Nodoid (g = 1)', s: { alphas: [polar(0.45, 0)], theta0: 0, width: 10, height: 7.663434 } },
  { name: 'Twizzler (g = 1)', s: { alphas: [polar(0.45, 0)], theta0: PI / 2, width: 8, height: 8 } },
  { name: 'Bent tube (g = 2)', s: { alphas: [polar(0.49, 1), polar(0.49, -1)], theta0: PI / 2, width: 6, height: 6 } },
  { name: 'Nearly a cylinder (g = 2)', s: { alphas: [polar(0.85, 0), polar(0.85, PI / 2)], theta0: 0.7, width: 10, height: 10 } },
  { name: 'Bubbleton, two lobes', s: { alphas: [polar(RES(2), 0), polar(RES(2), 0)], theta0: 0, width: PI, height: 3.4, gridStep: PI / 48, hmax: 0.02, res: 360 } },
  { name: 'Bubbleton, three lobes', s: { alphas: [polar(RES(3), 0), polar(RES(3), 0)], theta0: 0, width: PI, height: 1.6, gridStep: PI / 60, hmax: 0.02, res: 512 } },
  { name: 'Bubbleton, four lobes', s: { alphas: [polar(RES(4), 0), polar(RES(4), 0)], theta0: 0, width: PI, height: 1.2, gridStep: PI / 80, hmax: 0.02, res: 512 } },
  // Wente torus: McIntosh's spectral data (arXiv math/0407248, Remark 1: a = 0.1413 ± 0.1018i), refined
  // by Newton on the closing conditions. The period lattice of zeta is rhombic, <w1, w2> with w1 = 2.31698
  // real and 2 w2 - w1 = 4.47988 i; E(w1) turns by 2/3 about an axis, E(2 w2 - w1) = id, so the torus
  // lattice is <3 w1, w1 + w2>, whose fundamental domain is the brick 3 w1 x (4.47988 / 2).
  { name: 'Wente torus (three lobes)', s: { alphas: [[0.1412634686, 0.1017768953], [0.1412634686, -0.1017768953]], theta0: 0, width: 6.950948617, height: 2.239937825, gridStep: 6.950948617 / 72, res: 360 } },
  { name: 'Threefold (g = 3)', s: { alphas: [polar(0.4, 0), polar(0.4, 2 * PI / 3), polar(0.4, -2 * PI / 3)], theta0: 0.3, tau: [0], width: 6, height: 6 } },
  { name: 'Genus 4', s: { alphas: [polar(0.54, 0.38), polar(0.4, PI / 2), polar(0.58, -2.6), polar(0.63, -1.25)], theta0: 1, tau: [0.8, -0.5], width: 6, height: 6 } },
];

const state = {
  alphas: [polar(0.49, 1), polar(0.49, -1)], theta0: PI / 2,
  ...VIEW_DEFAULTS, width: 6, height: 6,
  res: 240, hmax: 0.05, colour: 'side', grid: true, H: 0.5, adapt: true,
  front: '#6f8fb0', back: '#e8c9a0', bg: '#f5f0e4',
};

const genus = () => state.alphas.length;
function syncTau() {
  const n = shapeFlows(genus()).length;
  state.tau = Array.from({ length: n }, (_, i) => state.tau[i] || 0);
}

// ------------------------------------------------------------------ URL hash

const fmt = (x, d = 4) => String(+x.toFixed(d));
function writeHash() {
  const p = new URLSearchParams();
  p.set('a', state.alphas.map(([re, im]) => `${fmt(re, 10)},${fmt(im, 10)}`).join(';'));
  p.set('l', fmt(state.theta0));
  if (state.tau.length) p.set('t', state.tau.map((x) => fmt(x)).join(','));
  p.set('c', state.z0.map((x) => fmt(x)).join(','));
  p.set('d', `${fmt(state.width, 8)},${fmt(state.height, 8)}`);
  if (!state.curv) p.set('p', fmt(state.phi));
  if (state.colour !== 'side') p.set('m', state.colour);
  p.set('g', fmt(state.gridStep));
  if (state.H !== 0.5) p.set('H', fmt(state.H));
  if (state.hmax !== 0.05) p.set('h', fmt(state.hmax));
  if (state.res !== 240) p.set('r', state.res);
  // commas and semicolons are safe in a fragment; keep the link readable
  history.replaceState(null, '', '#' + p.toString().replace(/%2C/g, ',').replace(/%3B/g, ';'));
}
function readHash() {
  const p = new URLSearchParams(location.hash.slice(1));
  if (!p.has('a')) return false;
  const nums = (s) => (s || '').split(',').filter((x) => x !== '').map(Number);
  try {
    const a = p.get('a');
    state.alphas = a ? a.split(';').map((s) => nums(s)).filter((v) => v.length === 2 && v.every(Number.isFinite)) : [];
    state.alphas = state.alphas.slice(0, MAX_GENUS);
    if (p.has('l')) state.theta0 = +p.get('l') || 0;
    state.tau = nums(p.get('t'));
    if (p.has('c')) { const c = nums(p.get('c')); if (c.length === 2) state.z0 = c; }
    if (p.has('d')) { const d = nums(p.get('d')); if (d.length === 2 && d.every((x) => x > 0)) [state.width, state.height] = d; }
    state.curv = !p.has('p');
    if (p.has('p')) state.phi = +p.get('p') || 0;
    if (p.has('m')) state.colour = p.get('m');
    if (p.has('g')) state.gridStep = +p.get('g') || state.gridStep;
    if (p.has('H')) state.H = Math.max(0.01, +p.get('H') || 0.5);
    if (p.has('h')) state.hmax = +p.get('h') || 0.05;
    if (p.has('r')) state.res = Math.min(1024, Math.max(16, +p.get('r') || 240));
  } catch (e) {
    console.warn('bad URL state', e);
  }
  return true;
}

// ------------------------------------------------------------------ computation scheduling

const viewer = new Viewer($('view'));
let worker = null, busy = false, busyFull = false, pending = null, reqId = 0, shownId = 0;
let fullTimer = null, needFrame = true, last = null;
const jobInfo = new Map(); // job id -> what the job computed

function makeWorker() {
  worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  worker.onmessage = onResult;
  worker.onerror = (e) => showStatus(`worker error: ${e.message}`, true);
}

function gridDims(res) {
  const m = Math.max(state.width, state.height);
  return [Math.max(2, Math.round((res * state.width) / m)), Math.max(2, Math.round((res * state.height) / m))];
}

// everything that determines the surface (not the mesh or the accuracy)
const surfaceKey = () => JSON.stringify([state.alphas, state.theta0, state.tau, state.z0,
  state.curv ? null : state.phi, state.width, state.height]);

function params(full) {
  const [nx, ny] = gridDims(full ? state.res : 72);
  const p = {
    alphas: state.alphas, theta0: state.theta0, tau: state.tau, z0: state.z0,
    phi: state.curv ? null : state.phi, width: state.width, height: state.height,
    nx, ny, hmax: full ? state.hmax : Math.max(state.hmax, 0.1),
  };
  // full meshes: put the rows and columns where the surface is, using the preview's metric
  if (full && state.adapt && last && last.key === surfaceKey() && last.preview) {
    p.sCoords = equidistribute(last, 's', nx);
    p.tCoords = equidistribute(last, 't', ny);
  }
  return p;
}

/**
 * Grid coordinates (relative to z0) along one axis of a result, spaced so that each interval carries
 * about the same length on the surface: density max e^u across the other axis, plus a floor.
 */
function equidistribute(res, axis, N) {
  const { nx, ny, u } = res;
  const n = axis === 's' ? nx : ny, m = axis === 's' ? ny : nx;
  const L = axis === 's' ? state.width : state.height;
  const w = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let mx = 0;
    for (let j = 0; j < m; j++) {
      const v = axis === 's' ? u[j * nx + i] : u[i * nx + j];
      if (Number.isFinite(v)) mx = Math.max(mx, Math.exp(v));
    }
    w[i] = mx;
  }
  const mean = w.reduce((a, b) => a + b, 0) / n || 1;
  const x = (i) => -L / 2 + (n > 1 ? (i * L) / (n - 1) : 0);
  const W = new Float64Array(n); // cumulative length (trapezoid)
  for (let i = 1; i < n; i++) W[i] = W[i - 1] + 0.5 * (w[i - 1] + w[i] + 0.5 * mean) * (x(i) - x(i - 1));
  const out = new Float64Array(N);
  let k = 0;
  for (let q = 0; q < N; q++) {
    const target = N > 1 ? (q * W[n - 1]) / (N - 1) : 0;
    while (k < n - 2 && W[k + 1] < target) k++;
    const f = W[k + 1] > W[k] ? (target - W[k]) / (W[k + 1] - W[k]) : 0;
    out[q] = x(k) + Math.min(1, Math.max(0, f)) * (x(k + 1) - x(k));
  }
  out[0] = -L / 2; out[N - 1] = L / 2;
  return out;
}

function request(full) {
  pending = { id: ++reqId, full };
  if (busy && busyFull && !full) { // a slow job is in the way of interactive feedback: drop it
    worker.terminate();
    makeWorker();
    busy = false;
  }
  pump();
}

function pump() {
  if (busy || !pending) return;
  const job = pending;
  pending = null;
  busy = true;
  busyFull = job.full;
  $('busy').classList.add('on');
  const p = params(job.full);
  jobInfo.set(job.id, { key: surfaceKey(), preview: !job.full });
  worker.postMessage({ id: job.id, params: p });
}

function onResult(e) {
  busy = false;
  const { id, result, error } = e.data;
  if (error) showStatus(error, true);
  else if (id > shownId) { shownId = id; Object.assign(result, jobInfo.get(id)); show(result); }
  jobInfo.delete(id);
  if (!pending) $('busy').classList.remove('on');
  pump();
}

/** Called after every state change. dragging: more changes are coming, so only preview. */
function changed(dragging = false) {
  whithamCheckAnchor();
  requestFamily();
  request(false);
  clearTimeout(fullTimer);
  fullTimer = setTimeout(() => request(true), dragging || anim ? 350 : 0);
  clearTimeout(changed.hashTimer);
  changed.hashTimer = setTimeout(writeHash, 300);
}

// ------------------------------------------------------------------ display

function colourValues(res) {
  if (state.colour === 'side') return { vals: null, range: [0, 1] };
  const N = res.u.length;
  const vals = new Float32Array(N);
  for (let v = 0; v < N; v++) {
    const u = res.u[v];
    // K = (2H)^2 (1 - e^{-4u}) / 4 after the homothety by 1/(2H)
    vals[v] = state.colour === 'K' ? state.H * state.H * (1 - Math.exp(-4 * u)) : u;
  }
  // robust range: K has a long negative tail where u is small
  const sorted = vals.filter(Number.isFinite).sort();
  const q = (f) => sorted[Math.min(sorted.length - 1, Math.floor(f * sorted.length))];
  if (!sorted.length) return { vals, range: [0, 1] };
  const H2 = state.H * state.H; // K < H^2, with equality only at umbilics
  return { vals, range: state.colour === 'K' ? [Math.min(q(0.1), -1e-3 * H2), H2] : [q(0.02), q(0.98)] };
}

function show(res) {
  last = res;
  const { vals, range } = colourValues(res);
  viewer.setSurface(res, vals);
  viewer.setStyle({ mode: { side: 0, u: 1, K: 2 }[state.colour], range });
  if (needFrame) { viewer.frame(); needFrame = false; }
  widget.set({ divisor: res.divisor });
  const g = genus();
  const argQ = hopfArg(state.alphas, state.theta0);
  const size = Math.max(res.bbox[3] - res.bbox[0], res.bbox[4] - res.bbox[1], res.bbox[5] - res.bbox[2]) / (2 * state.H);
  showStatus(
    `genus ${g} · H = ${fmt(state.H, 3)} · κ₀ = ${fmt(kappa0(state.alphas), 3)} · arg Q = ${fmt(wrapPi(argQ), 3)}` +
    ` · ${res.nx}×${res.ny} in ${res.stats.ms.toFixed(0)} ms · extent ${size.toPrecision(3)}` +
    (res.stats.detDefect > 1e-6 ? ` · det drift ${res.stats.detDefect.toExponential(1)}` : ''),
    false,
  );
}

function showStatus(text, err) {
  const s = $('status');
  s.textContent = text;
  s.classList.toggle('err', !!err);
}

const wrapPi = (x) => x - 2 * PI * Math.round(x / (2 * PI));

// ------------------------------------------------------------------ spectral widget and alpha list

const widget = new SpectralWidget($('lam'), (s, dragging) => {
  const gOld = genus();
  state.alphas = s.alphas;
  state.theta0 = s.theta0;
  if (genus() !== gOld) genusChanged();
  else refreshAlphaValues();
  changed(dragging);
});

function genusChanged() {
  whithamRecentre();
  syncTau();
  needFrame = true;
  buildAlphaList();
  buildTau();
  widget.set({ alphas: state.alphas, theta0: state.theta0, divisor: [] });
}

function buildAlphaList() {
  const box = $('alphaList');
  box.innerHTML = '';
  state.alphas.forEach((_, i) => {
    const row = document.createElement('div');
    row.className = 'alpharow';
    row.innerHTML = `<span class="name">α<sub>${i + 1}</sub></span>
      <label>r</label><input type="number" step="0.01" min="0.03" max="0.97" data-k="r">
      <label>arg</label><input type="number" step="0.01" data-k="t">
      <button type="button" title="remove">×</button>`;
    const [r, t] = row.querySelectorAll('input');
    const update = () => {
      const rr = Math.min(0.97, Math.max(0.03, +r.value || 0.03));
      state.alphas[i] = polar(rr, +t.value || 0);
      widget.set({ alphas: state.alphas });
      changed(false);
    };
    r.addEventListener('change', update);
    t.addEventListener('change', update);
    row.querySelector('button').addEventListener('click', () => {
      state.alphas.splice(i, 1);
      genusChanged();
      changed(false);
    });
    box.appendChild(row);
  });
  $('genus').textContent = `genus ${genus()}`;
  refreshAlphaValues();
}

function refreshAlphaValues() {
  const rows = $('alphaList').children;
  state.alphas.forEach(([re, im], i) => {
    const [r, t] = rows[i].querySelectorAll('input');
    if (document.activeElement !== r) r.value = fmt(Math.hypot(re, im), 3);
    if (document.activeElement !== t) t.value = fmt(Math.atan2(im, re), 3);
  });
  if (document.activeElement !== $('theta0')) $('theta0').value = fmt(wrapPi(state.theta0), 3);
}

$('theta0').addEventListener('change', () => {
  state.theta0 = +$('theta0').value || 0;
  widget.set({ theta0: state.theta0 });
  changed(false);
});
$('addAlpha').addEventListener('click', () => {
  if (genus() >= MAX_GENUS) return;
  // place the new point where it is well separated from the others
  let best = null, bestD = -1;
  for (let k = 0; k < 24; k++) {
    const z = polar(0.5, (2 * PI * k) / 24 + 0.1);
    const d = Math.min(...state.alphas.map((a) => Math.hypot(a[0] - z[0], a[1] - z[1])), 9);
    if (d > bestD) { bestD = d; best = z; }
  }
  state.alphas.push(best);
  genusChanged();
  changed(false);
});
$('remAlpha').addEventListener('click', () => {
  if (!genus()) return;
  state.alphas.pop();
  genusChanged();
  changed(false);
});

// ------------------------------------------------------------------ sliders

function slider(parent, { label, min, max, step, get, set, digits = 3, play, recompute = true }) {
  const after = recompute ? changed : writeHashSoon;
  const row = document.createElement('div');
  row.className = 'row';
  row.innerHTML = `${label ? `<label>${label}</label>` : ''}<input type="range"><input type="number">`;
  const [rng, num] = row.querySelectorAll('input');
  Object.assign(rng, { min, max, step });
  num.step = step;
  const refresh = () => {
    const v = get();
    rng.value = v;
    if (document.activeElement !== num) num.value = fmt(v, digits);
  };
  rng.addEventListener('input', () => { set(+rng.value); refresh(); after(true); });
  num.addEventListener('change', () => {
    const v = +num.value;
    if (!Number.isFinite(v)) return;
    set(v);
    refresh();
    after(false);
  });
  if (play) {
    const b = document.createElement('button');
    b.type = 'button';
    b.className = 'play';
    b.textContent = '▶';
    b.title = 'animate';
    b.addEventListener('click', () => toggleAnim(play, b));
    row.appendChild(b);
  }
  parent.appendChild(row);
  refresh();
  return { refresh, row };
}

let tauSliders = [];
function buildTau() {
  const box = $('tauList');
  box.innerHTML = '';
  if (anim && anim.key === 'tau') stopAnim();
  tauSliders = state.tau.map((_, i) => slider(box, {
    label: `τ<sub>${i + 1}</sub>`, min: -2 * PI, max: 2 * PI, step: 0.01,
    get: () => state.tau[i], set: (v) => { state.tau[i] = v; },
    play: { key: 'tau', i },
  }));
  const g = genus();
  $('tauNote').textContent = g >= 3
    ? `Times along the ${g - 2} shape-changing isospectral flows (the other two flows are translations of the domain).`
    : 'For genus ≤ 2 every isospectral deformation is a translation of the domain: use the domain centre below.';
}

const domainSliders = [
  slider($('domainRows'), { label: 'centre x', min: -20, max: 20, step: 0.01, get: () => state.z0[0], set: (v) => { state.z0[0] = v; } }),
  slider($('domainRows'), { label: 'centre y', min: -20, max: 20, step: 0.01, get: () => state.z0[1], set: (v) => { state.z0[1] = v; } }),
  slider($('domainRows'), { label: 'width', min: 0.5, max: 40, step: 0.05, get: () => state.width, set: (v) => { state.width = Math.max(0.05, v); } }),
  slider($('domainRows'), { label: 'height', min: 0.5, max: 40, step: 0.05, get: () => state.height, set: (v) => { state.height = Math.max(0.05, v); } }),
];
const phiSlider = slider($('phiRow'), { label: 'grid angle', min: -PI, max: PI, step: 0.01, get: () => state.phi, set: (v) => { state.phi = v; } });
const appearanceSliders = [
  slider($('appearanceRows'), {
    label: 'line spacing', min: 0.02, max: 2, step: 0.01, get: () => state.gridStep,
    set: (v) => { state.gridStep = Math.max(0.005, v); viewer.setStyle({ gridStep: state.gridStep }); },
    recompute: false,
  }),
  slider($('appearanceRows'), {
    label: 'line width', min: 0.2, max: 2.5, step: 0.05, get: () => state.lineW,
    set: (v) => { state.lineW = v; viewer.setStyle({ lineWidth: v }); },
    recompute: false,
  }),
];
function writeHashSoon() { clearTimeout(changed.hashTimer); changed.hashTimer = setTimeout(writeHash, 300); }

$('curv').addEventListener('change', () => {
  state.curv = $('curv').checked;
  $('phiRow').classList.toggle('hide', state.curv);
  changed(false);
});

// ------------------------------------------------------------------ mean curvature (a homothety in R^3)

function recolour() {
  if (!last) return;
  const { vals, range } = colourValues(last);
  if (vals) viewer.setColourValues(vals);
  viewer.setStyle({ mode: { side: 0, u: 1, K: 2 }[state.colour], range });
}
const Hslider = slider($('Hrow'), {
  label: 'H', min: 0.05, max: 3, step: 0.01, get: () => state.H,
  set: (v) => {
    state.H = Math.max(0.01, v);
    viewer.setScale(1 / (2 * state.H));
    recolour();
    if (last) show(last);
  },
  recompute: false,
});

// ------------------------------------------------------------------ Whitham deformation

// curve: the Whitham curve through the anchor (created lazily); key: the alphas it last produced, so
// that any other edit of the spectral data re-anchors; base: the domain at s = 0.
const whitham = { curve: null, s: 0, key: null, base: null, z: [1, 0], applied: null };
const domainKey = () => JSON.stringify([state.z0, state.width, state.height]);

function whithamCheckAnchor() {
  if (whitham.key !== null && whitham.key !== JSON.stringify(state.alphas)) whithamRecentre();
}
function whithamRecentre() {
  family.key = null; // the family's s = 0 moves here: retrace
  whitham.curve = null;
  whitham.s = 0;
  whitham.key = null;
  whitham.z = [1, 0];
  whithamSlider.refresh();
  whithamNote();
}
function whithamNote() {
  const g = genus();
  $('whithamBox').querySelectorAll('input, button').forEach((el) => { el.disabled = g < 2; });
  if (g < 2) { $('whithamNote').textContent = 'Needs genus ≥ 2.'; return; }
  const F = family.data;
  if (!F || F.points.length < 2) {
    $('whithamNote').textContent = 'Moves the branch points keeping the conformal type of the period lattice.';
    return;
  }
  // W at the current point: on the family, interpolated in s
  const P = F.points, s = whitham.s;
  let i = 1;
  while (i < P.length - 1 && P[i].s < s) i++;
  const f = P[i].s > P[i - 1].s ? Math.min(1, Math.max(0, (s - P[i - 1].s) / (P[i].s - P[i - 1].s))) : 0;
  const lines = [`𝒲 = ${fmt(P[i - 1].W + f * (P[i].W - P[i - 1].W), 4)}`];
  for (const c of F.critical) lines.push(`𝒲 critical at s = ${fmt(P[c].s, 3)}, 𝒲 = ${fmt(P[c].W, 4)}`);
  $('whithamNote').textContent = lines.join('\n');
}

// the family (Whitham curve through the current data) for the λ-plane, traced in its own worker so it
// never delays the surface; the latest request wins. Moving along the curve with the slider doesn't
// change the family, so only other edits of the spectral data (or recentring) trigger a new trace.
const family = { worker: null, busy: false, pending: null, key: null, data: null };
function requestFamily() {
  if (!$('showFamily').checked || genus() < 2) {
    family.pending = null;
    if (family.data) { family.data = null; family.key = null; widget.set({ family: null }); }
    return;
  }
  const key = JSON.stringify(state.alphas);
  if (key === family.key || (family.data && key === whitham.key)) return;
  family.key = key;
  family.pending = { key, alphas: state.alphas.map((a) => a.slice()) };
  pumpFamily();
}
function pumpFamily() {
  if (family.busy || !family.pending) return;
  if (!family.worker) {
    family.worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
    family.worker.onmessage = (e) => {
      family.busy = false;
      const { result, error } = e.data;
      // show every result, even if a newer request is queued: during a drag that keeps the paths live
      if ($('showFamily').checked && genus() >= 2) {
        family.data = error ? null : result;
        widget.set({ family: family.data });
        whithamNote();
      }
      pumpFamily();
    };
  }
  const job = family.pending;
  family.pending = null;
  family.busy = true;
  family.worker.postMessage({ id: 0, kind: 'family', params: { alphas: job.alphas, smax: 1.5, maxSteps: 80 } });
}
function whithamMove(s) {
  if (genus() < 2) return s;
  if (!whitham.curve) {
    try {
      whitham.curve = new WhithamCurve(state.alphas);
    } catch (e) {
      $('whithamNote').textContent = String(e.message || e);
      return 0;
    }
    whitham.base = { z0: state.z0.slice(), width: state.width, height: state.height };
    whitham.applied = domainKey();
  }
  // a domain edited by hand since the last move becomes the new base (pulled back to s = 0)
  if (whitham.applied !== domainKey()) {
    const [zr, zi] = whitham.z, m = Math.hypot(zr, zi);
    whitham.base = {
      z0: [state.z0[0] * zr - state.z0[1] * zi, state.z0[0] * zi + state.z0[1] * zr],
      width: state.width * m, height: state.height * m,
    };
  }
  const r = whitham.curve.at(s);
  state.alphas = r.alphas;
  whitham.s = r.s;
  whitham.z = r.z;
  whitham.key = JSON.stringify(state.alphas);
  if ($('whithamDomain').checked) {
    // Gamma(s) = Gamma(0) / z
    const [zr, zi] = r.z, d = zr * zr + zi * zi, m = Math.sqrt(d);
    const b = whitham.base;
    state.z0 = [(b.z0[0] * zr + b.z0[1] * zi) / d, (b.z0[1] * zr - b.z0[0] * zi) / d];
    state.width = b.width / m;
    state.height = b.height / m;
    for (const sl of domainSliders) sl.refresh();
  }
  whitham.applied = domainKey();
  widget.set({ alphas: state.alphas });
  refreshAlphaValues();
  whithamNote();
  return r.s;
}
const whithamSlider = slider($('whithamRow'), {
  label: '', min: -1.5, max: 1.5, step: 0.002, get: () => whitham.s,
  set: (v) => { whithamMove(v); },
  play: { key: 'whitham' },
});
{
  // the recentre button sits beside the play button, same size
  const b = document.createElement('button');
  b.type = 'button';
  b.className = 'play';
  b.id = 'whithamReset';
  b.textContent = '↻';
  b.title = 'recentre: make the current spectral curve the centre of the slider';
  b.addEventListener('click', () => { whithamRecentre(); requestFamily(); });
  whithamSlider.row.appendChild(b);
}
$('showFamily').addEventListener('change', () => { family.key = null; requestFamily(); whithamNote(); });
$('whithamDomain').addEventListener('change', () => { if (whitham.curve) { whithamMove(whitham.s); changed(false); } });

// ------------------------------------------------------------------ closing up

$('closeBtn').addEventListener('click', () => {
  const w = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  const info = $('closeInfo');
  info.textContent = 'searching for periods…';
  $('closeBtn').disabled = true;
  const done = () => { w.terminate(); $('closeBtn').disabled = false; };
  w.onerror = (e) => { info.textContent = `error: ${e.message}`; done(); };
  w.onmessage = (e) => {
    done();
    if (e.data.error) { info.textContent = e.data.error; return; }
    const parts = [];
    let changedAny = false;
    for (const d of ['s', 't']) {
      const r = e.data.result[d];
      const name = d === 's' ? 'width' : 'height';
      if (!r) { parts.push(`${d}: no period`); continue; }
      if (r.q) {
        state[name] = r.q * r.T;
        changedAny = true;
        parts.push(`${d}: closes after ${r.q === 1 ? 'one period' : `${r.q} periods`} (${fmt(r.q * r.T, 5)})`);
      } else {
        parts.push(`${d}: period ${fmt(r.T, 4)}, turns ${fmt(r.angle, 4)}, pitch ${fmt(r.pitch, 3)}`);
      }
    }
    info.textContent = parts.join(' · ');
    if (changedAny) {
      for (const s of domainSliders) s.refresh();
      saveHomeDomain();
      needFrame = true;
      changed(false);
    }
  };
  w.postMessage({ id: 0, kind: 'close', params: { ...params(true), hmax: Math.min(state.hmax, 0.01) } });
});

// ------------------------------------------------------------------ appearance

$('colour').addEventListener('change', () => {
  state.colour = $('colour').value;
  recolour();
  writeHashSoon();
});
$('grid').addEventListener('change', () => { state.grid = $('grid').checked; viewer.setStyle({ grid: state.grid }); });
$('front').addEventListener('input', () => viewer.setStyle({ front: $('front').value }));
$('back').addEventListener('input', () => viewer.setStyle({ back: $('back').value }));
$('bg').addEventListener('input', () => setBackground($('bg').value));
function setBackground(c) {
  state.bg = c;
  $('stage').style.background = c;
  viewer.setBackground(c);
}
$('adapt').addEventListener('change', () => { state.adapt = $('adapt').checked; changed(false); });
$('res').addEventListener('change', () => { state.res = +$('res').value; changed(false); });
$('acc').addEventListener('change', () => { state.hmax = +$('acc').value; changed(false); });

// ------------------------------------------------------------------ animation

let anim = null, animBtn = null, animLast = 0;
function toggleAnim(a, btn) {
  const same = anim && anim.key === a.key && anim.i === a.i;
  stopAnim();
  if (same) return;
  anim = a;
  animBtn = btn;
  btn.classList.add('on');
  btn.textContent = '■';
  animLast = performance.now();
  requestAnimationFrame(animStep);
}
function stopAnim() {
  if (animBtn) { animBtn.classList.remove('on'); animBtn.textContent = '▶'; }
  anim = null;
  animBtn = null;
  changed(false);
}
function animStep(now) {
  if (!anim) return;
  const dt = Math.min(0.1, (now - animLast) / 1000);
  animLast = now;
  if (anim.key === 'theta0') {
    state.theta0 = wrapPi(state.theta0 + 0.3 * dt);
    widget.set({ theta0: state.theta0 });
    refreshAlphaValues();
  } else if (anim.key === 'whitham') {
    if (!busy) {
      anim.dir = anim.dir || 1;
      const want = whitham.s + anim.dir * 0.15 * dt;
      const got = whithamMove(Math.max(-1.5, Math.min(1.5, want)));
      if (Math.abs(got - want) > 1e-9 || Math.abs(got) >= 1.5) anim.dir = -anim.dir;
      whithamSlider.refresh();
    }
  } else {
    let v = state.tau[anim.i] + 0.4 * dt;
    if (v > 2 * PI) v -= 4 * PI;
    state.tau[anim.i] = v;
    tauSliders[anim.i].refresh();
  }
  // only ask for a new surface when the worker has caught up
  if (!busy) changed(true);
  requestAnimationFrame(animStep);
}
$('playLam').addEventListener('click', () => toggleAnim({ key: 'theta0' }, $('playLam')));

// ------------------------------------------------------------------ presets and export

const presetSel = $('preset');
presetSel.innerHTML = '<option value="">choose a preset…</option>' +
  PRESETS.map((p, i) => `<option value="${i}">${p.name}</option>`).join('');
presetSel.addEventListener('change', () => {
  const p = PRESETS[+presetSel.value];
  if (!p) return;
  if (anim) stopAnim();
  Object.assign(state, structuredClone(VIEW_DEFAULTS), structuredClone(p.s));
  syncTau();
  refreshAll();
  saveHomeDomain();
  viewer.resetOrientation();
  needFrame = true;
  changed(false);
  presetSel.value = '';
});

$('frameBtn').addEventListener('click', () => viewer.frame());
$('viewReset').addEventListener('click', () => { viewer.resetOrientation(); viewer.frame(); });

// 'reset domain' returns to the domain of the last preset, link or close-up
let homeDomain = null;
const saveHomeDomain = () => { homeDomain = structuredClone({ z0: state.z0, width: state.width, height: state.height, curv: state.curv, phi: state.phi }); };
$('domainReset').addEventListener('click', () => {
  Object.assign(state, structuredClone(homeDomain));
  for (const sl of [...domainSliders, phiSlider]) sl.refresh();
  $('curv').checked = state.curv;
  $('phiRow').classList.toggle('hide', state.curv);
  if (whitham.curve) whitham.applied = null; // becomes the new base at the next Whitham move
  needFrame = true;
  changed(false);
});
$('pngBtn').addEventListener('click', () => download(viewer.snapshot(), 'cmc-surface.png'));
$('objBtn').addEventListener('click', () => {
  const url = URL.createObjectURL(new Blob([viewer.toOBJ()], { type: 'text/plain' }));
  download(url, 'cmc-surface.obj');
  setTimeout(() => URL.revokeObjectURL(url), 1000);
});
$('linkBtn').addEventListener('click', async () => {
  writeHash();
  try {
    await navigator.clipboard.writeText(location.href);
    $('linkBtn').textContent = 'copied';
  } catch {
    $('linkBtn').textContent = 'see address bar';
  }
  setTimeout(() => { $('linkBtn').textContent = 'copy link'; }, 1500);
});
function download(href, name) {
  const a = document.createElement('a');
  a.href = href;
  a.download = name;
  a.click();
}

// ------------------------------------------------------------------ start

function refreshAll() {
  buildAlphaList();
  buildTau();
  widget.set({ alphas: state.alphas, theta0: state.theta0, divisor: [] });
  for (const s of [...domainSliders, phiSlider, ...appearanceSliders, Hslider]) s.refresh();
  $('res').value = String(state.res);
  if (!$('res').value) { $('res').add(new Option(String(state.res), String(state.res))); $('res').value = String(state.res); }
  $('acc').value = String(state.hmax);
  if (!$('acc').value) { $('acc').add(new Option(`h ≤ ${state.hmax}`, String(state.hmax))); $('acc').value = String(state.hmax); }
  viewer.setScale(1 / (2 * state.H));
  $('closeInfo').textContent = '';
  whithamRecentre();
  $('curv').checked = state.curv;
  $('phiRow').classList.toggle('hide', state.curv);
  $('colour').value = state.colour;
  viewer.setStyle({ gridStep: state.gridStep, lineWidth: state.lineW, grid: state.grid });
}

// a link pasted into the address bar of an open page (our own writes use replaceState, which is silent)
window.addEventListener('hashchange', () => {
  if (!readHash()) return;
  if (anim) stopAnim();
  syncTau();
  refreshAll();
  saveHomeDomain();
  needFrame = true;
  changed(false);
});

readHash();
syncTau();
makeWorker();
refreshAll();
saveHomeDomain();
setBackground(state.bg);
changed(false);
