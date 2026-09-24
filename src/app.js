// Glue: state, panel controls, progressive recomputation in a worker, presets, URL state, export.

import { Viewer } from './viewer.js';
import { SpectralWidget, MAX_GENUS } from './spectral-widget.js';
import { shapeFlows, kappa0, hopfArg } from './cmc/cmc.js';

const $ = (id) => document.getElementById(id);
const polar = (r, t) => [r * Math.cos(t), r * Math.sin(t)];
const PI = Math.PI;

// ------------------------------------------------------------------ state

const VIEW_DEFAULTS = {
  tau: [], z0: [0, 0], width: 12, height: 12, curv: true, phi: 0,
  gridStep: 0.25, lineW: 0.6,
};
// Widths of the rotational examples are exact closing periods of the parallels, pi / sqrt(det X(lam0)).
const PRESETS = [
  { name: 'Round cylinder (g = 0)', s: { alphas: [], theta0: 0, width: PI, height: 10, gridStep: PI / 16 } },
  { name: 'Unduloid (g = 1)', s: { alphas: [polar(0.45, 0)], theta0: PI, width: 2.90682, height: 14 } },
  { name: 'Unduloid near a sphere chain', s: { alphas: [polar(0.08, 0)], theta0: PI, width: 1.645512, height: 14, gridStep: 0.1 } },
  { name: 'Nodoid (g = 1)', s: { alphas: [polar(0.45, 0)], theta0: 0, width: 10, height: 7.663434 } },
  { name: 'Twizzler (g = 1)', s: { alphas: [polar(0.45, 0)], theta0: PI / 2, width: 8, height: 8 } },
  { name: 'Bent tube (g = 2)', s: { alphas: [polar(0.49, 1), polar(0.49, -1)], theta0: PI / 2, width: 6, height: 6 } },
  { name: 'Nearly a cylinder (g = 2)', s: { alphas: [polar(0.85, 0), polar(0.85, PI / 2)], theta0: 0.7, width: 10, height: 10 } },
  { name: 'Threefold (g = 3)', s: { alphas: [polar(0.4, 0), polar(0.4, 2 * PI / 3), polar(0.4, -2 * PI / 3)], theta0: 0.3, tau: [0], width: 6, height: 6 } },
  { name: 'Genus 4', s: { alphas: [polar(0.54, 0.38), polar(0.4, PI / 2), polar(0.58, -2.6), polar(0.63, -1.25)], theta0: 1, tau: [0.8, -0.5], width: 6, height: 6 } },
];

const state = {
  alphas: [polar(0.49, 1), polar(0.49, -1)], theta0: PI / 2,
  ...VIEW_DEFAULTS, width: 6, height: 6,
  res: 240, hmax: 0.05, colour: 'side', grid: true,
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
  p.set('a', state.alphas.map(([re, im]) => `${fmt(re)},${fmt(im)}`).join(';'));
  p.set('l', fmt(state.theta0));
  if (state.tau.length) p.set('t', state.tau.map((x) => fmt(x)).join(','));
  p.set('c', state.z0.map((x) => fmt(x)).join(','));
  p.set('d', `${fmt(state.width)},${fmt(state.height)}`);
  if (!state.curv) p.set('p', fmt(state.phi));
  if (state.colour !== 'side') p.set('m', state.colour);
  p.set('g', fmt(state.gridStep));
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
  } catch (e) {
    console.warn('bad URL state', e);
  }
  return true;
}

// ------------------------------------------------------------------ computation scheduling

const viewer = new Viewer($('view'));
let worker = null, busy = false, busyFull = false, pending = null, reqId = 0, shownId = 0;
let fullTimer = null, needFrame = true, last = null;

function makeWorker() {
  worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  worker.onmessage = onResult;
  worker.onerror = (e) => showStatus(`worker error: ${e.message}`, true);
}

function gridDims(res) {
  const m = Math.max(state.width, state.height);
  return [Math.max(2, Math.round((res * state.width) / m)), Math.max(2, Math.round((res * state.height) / m))];
}

function params(full) {
  const [nx, ny] = gridDims(full ? state.res : 72);
  return {
    alphas: state.alphas, theta0: state.theta0, tau: state.tau, z0: state.z0,
    phi: state.curv ? null : state.phi, width: state.width, height: state.height,
    nx, ny, hmax: full ? state.hmax : Math.max(state.hmax, 0.1),
  };
}

function request(full) {
  pending = { id: ++reqId, full, params: params(full) };
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
  worker.postMessage({ id: job.id, params: job.params });
}

function onResult(e) {
  busy = false;
  const { id, result, error } = e.data;
  if (error) showStatus(error, true);
  else if (id > shownId) { shownId = id; show(result); }
  if (!pending) $('busy').classList.remove('on');
  pump();
}

/** Called after every state change. dragging: more changes are coming, so only preview. */
function changed(dragging = false) {
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
    vals[v] = state.colour === 'K' ? 0.25 * (1 - Math.exp(-4 * u)) : u;
  }
  // robust range: K has a long negative tail where u is small
  const sorted = vals.filter(Number.isFinite).sort();
  const q = (f) => sorted[Math.min(sorted.length - 1, Math.floor(f * sorted.length))];
  if (!sorted.length) return { vals, range: [0, 1] };
  return { vals, range: state.colour === 'K' ? [Math.min(q(0.1), -1e-3), 0.25] : [q(0.02), q(0.98)] };
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
  const size = Math.max(res.bbox[3] - res.bbox[0], res.bbox[4] - res.bbox[1], res.bbox[5] - res.bbox[2]);
  showStatus(
    `genus ${g} · κ₀ = ${fmt(kappa0(state.alphas), 3)} · arg Q = ${fmt(wrapPi(argQ), 3)}` +
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
  row.innerHTML = `<label>${label}</label><input type="range"><input type="number">`;
  const [rng, num] = row.querySelectorAll('input');
  Object.assign(rng, { min, max, step });
  num.step = step;
  const refresh = () => {
    const v = get();
    rng.value = v;
    if (document.activeElement !== num) num.value = fmt(v, digits);
  };
  rng.addEventListener('input', () => { set(+rng.value); num.value = fmt(+rng.value, digits); after(true); });
  num.addEventListener('change', () => {
    const v = +num.value;
    if (!Number.isFinite(v)) return;
    set(v);
    rng.value = v;
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

// ------------------------------------------------------------------ appearance

$('colour').addEventListener('change', () => {
  state.colour = $('colour').value;
  if (last) {
    const { vals, range } = colourValues(last);
    if (vals) viewer.setColourValues(vals);
    viewer.setStyle({ mode: { side: 0, u: 1, K: 2 }[state.colour], range });
  }
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
  viewer.resetOrientation();
  needFrame = true;
  changed(false);
  presetSel.value = '';
});

$('frameBtn').addEventListener('click', () => viewer.frame());
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
  for (const s of [...domainSliders, phiSlider, ...appearanceSliders]) s.refresh();
  $('curv').checked = state.curv;
  $('phiRow').classList.toggle('hide', state.curv);
  $('colour').value = state.colour;
  viewer.setStyle({ gridStep: state.gridStep, lineWidth: state.lineW, grid: state.grid });
}

readHash();
syncTau();
makeWorker();
refreshAll();
setBackground(state.bg);
changed(false);
