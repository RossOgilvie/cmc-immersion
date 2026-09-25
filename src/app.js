// Glue: state, panel controls, progressive recomputation in a worker, presets, URL state, export.

import { Viewer } from './viewer.js';
import { SpectralWidget, MAX_GENUS } from './spectral-widget.js';
import { shapeFlows, kappa0, hopfArg } from './cmc/cmc.js';
import { WhithamCurve, willmore } from './cmc/whitham.js';
import { commonRootsOnCircle } from './cmc/periods.js';

const $ = (id) => document.getElementById(id);
const polar = (r, t) => [r * Math.cos(t), r * Math.sin(t)];
const PI = Math.PI;

// ------------------------------------------------------------------ panel furniture

// borderless stroke glyphs for buttons
const GLYPH = {
  chev: '<svg width="10" height="10" viewBox="0 0 10 10" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round" stroke-linejoin="round"><path d="M3.5 2l3 3-3 3"/></svg>',
  play: '<svg width="14" height="14" viewBox="0 0 14 14"><path d="M4 2.5 L11.5 7 L4 11.5 Z" fill="currentColor"/></svg>',
  stop: '<svg width="14" height="14" viewBox="0 0 14 14"><rect x="3.5" y="3.5" width="7" height="7" fill="currentColor"/></svg>',
  reset: '<svg width="15" height="15" viewBox="0 0 15 15" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round" stroke-linejoin="round"><path d="M12 7.5a4.5 4.5 0 1 1-1.3-3.2"/><path d="M11.2 1.8v2.7H8.5"/></svg>',
  plus: '<svg width="12" height="12" viewBox="0 0 12 12" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round"><path d="M6 2v8M2 6h8"/></svg>',
  minus: '<svg width="12" height="12" viewBox="0 0 12 12" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round"><path d="M2 6h8"/></svg>',
  cross: '<svg width="11" height="11" viewBox="0 0 12 12" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round"><path d="M2.5 2.5l7 7M9.5 2.5l-7 7"/></svg>',
};
const glyphs = (root) => root.querySelectorAll('[data-glyph]').forEach((el) => { el.innerHTML = GLYPH[el.dataset.glyph]; });
glyphs(document);

// every section (and the help and branch-point list) collapses; the choice is remembered per browser
{
  let open = {};
  try { open = JSON.parse(localStorage.getItem('cmc-open')) || {}; } catch { /* no storage */ }
  document.querySelectorAll('.disc[aria-controls]').forEach((b) => {
    const target = $(b.getAttribute('aria-controls'));
    const setOpen = (v) => { b.setAttribute('aria-expanded', String(v)); target.hidden = !v; };
    const id = target.id;
    setOpen(id in open ? open[id] : b.getAttribute('aria-expanded') === 'true');
    b.addEventListener('click', () => {
      const v = b.getAttribute('aria-expanded') !== 'true';
      setOpen(v);
      open[id] = v;
      try { localStorage.setItem('cmc-open', JSON.stringify(open)); } catch { /* no storage */ }
      if (v) widget.draw(); // the canvas may have been laid out while hidden
    });
  });
}

/** Fixed decimals, a true minus sign, and no "−0.000". */
function fixed(x, d) {
  const t = Math.abs(x).toFixed(d);
  return (x < 0 && +t !== 0 ? '−' : '') + t;
}
/** Number typed into a field: accepts − and a trailing π (ignored: angle fields are already in π). */
function parseNum(str) {
  const t = String(str).replace(/−/g, '-').replace(/π/g, '').trim();
  if (t === '') return null;
  const v = Number(t);
  return Number.isFinite(v) ? v : null;
}

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
  alphas: [polar(0.49, 1), polar(0.49, -1)], theta0: 0,
  ...VIEW_DEFAULTS, width: 6, height: 6,
  res: 240, hmax: 0.05, colour: 'side', grid: true, H: 0.5, adapt: true,
  front: '#6f8fb0', back: '#e8c9a0', bg: '#f2eee3',
};

const genus = () => state.alphas.length;
function syncTau() {
  const n = shapeFlows(genus()).length;
  state.tau = Array.from({ length: n }, (_, i) => state.tau[i] || 0);
}

// ------------------------------------------------------------------ links
// The state goes into the URL only through "Copy link". A page opened from such a link (or with one
// pasted into its address bar) takes the state from the hash and then clears it from the address bar.

const fmt = (x, d = 4) => String(+x.toFixed(d));
function linkURL() {
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
  return `${location.origin}${location.pathname}${location.search}#${p.toString().replace(/%2C/g, ',').replace(/%3B/g, ';')}`;
}
const clearHash = () => history.replaceState(null, '', location.pathname + location.search);
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
  whithamSlider.refresh();
  requestFamily();
  updateCommonRoots();
  request(false);
  clearTimeout(fullTimer);
  fullTimer = setTimeout(() => request(true), dragging || anim ? 350 : 0);
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
    `Genus ${g} · H = ${fmt(state.H, 3)} · κ₀ = ${fmt(kappa0(state.alphas), 3)} · arg Q = ${fixed(wrapPi(argQ) / PI, 3)}π` +
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

function removeAlpha(i) {
  if (i < 0 || i >= genus()) return;
  state.alphas.splice(i, 1);
  genusChanged();
  changed(false);
}

function buildAlphaList() {
  const box = $('alphaList');
  box.innerHTML = '';
  state.alphas.forEach((_, i) => {
    const row = document.createElement('div');
    row.className = 'bprow';
    row.innerHTML = `<span class="name">α<sub>${i + 1}</sub></span>
      <span class="lab">r</span><input class="num" type="text" inputmode="decimal" autocomplete="off" aria-label="Modulus of α${i + 1}">
      <span class="lab">arg</span><input class="num" type="text" inputmode="decimal" autocomplete="off" aria-label="Argument of α${i + 1} in units of π">
      <span class="muted">π</span>
      <button class="gly" type="button" data-glyph="cross" aria-label="Remove α${i + 1}" title="Remove"></button>`;
    glyphs(row);
    const [r, t] = row.querySelectorAll('input');
    const update = () => {
      const [re, im] = state.alphas[i];
      const rv = parseNum(r.value), tv = parseNum(t.value);
      const rr = Math.min(0.97, Math.max(0.03, rv ?? Math.hypot(re, im)));
      state.alphas[i] = polar(rr, tv === null ? Math.atan2(im, re) : tv * PI);
      widget.set({ alphas: state.alphas });
      refreshAlphaValues();
      changed(false);
    };
    for (const inp of [r, t]) {
      inp.addEventListener('change', update);
      inp.addEventListener('keydown', (e) => { if (e.key === 'Enter') inp.blur(); });
    }
    row.querySelector('.gly').addEventListener('click', () => removeAlpha(i));
    box.appendChild(row);
  });
  $('genus').textContent = String(genus());
  $('remAlpha').disabled = !genus();
  $('addAlpha').disabled = genus() >= MAX_GENUS;
  refreshAlphaValues();
}

function refreshAlphaValues() {
  const rows = $('alphaList').children;
  state.alphas.forEach(([re, im], i) => {
    const [r, t] = rows[i].querySelectorAll('input');
    if (document.activeElement !== r) r.value = fixed(Math.hypot(re, im), 3);
    if (document.activeElement !== t) t.value = fixed(Math.atan2(im, re) / PI, 3);
  });
  if (document.activeElement !== $('theta0')) $('theta0').value = fixed(wrapPi(state.theta0) / PI, 3);
}

$('theta0').addEventListener('change', () => {
  const v = parseNum($('theta0').value);
  if (v !== null) state.theta0 = wrapPi(v * PI);
  widget.set({ theta0: state.theta0 });
  refreshAlphaValues();
  changed(false);
});
$('theta0').addEventListener('keydown', (e) => { if (e.key === 'Enter') $('theta0').blur(); });
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
$('remAlpha').addEventListener('click', () => removeAlpha(genus() - 1));

// ------------------------------------------------------------------ sliders

/**
 * A slider row: label, track, action buttons. The value floats above the thumb (display(v), default
 * fixed digits of v / unit plus suffix); clicking it opens a field for exact entry (enter(x) if given,
 * else set(x * unit)). snap() lists values that get a tick on the track and catch the thumb within a
 * few pixels. play: an animation toggle; reset: { title, onClick } for a second button.
 */
function slider(parent, { label, min: lo, max: hi, step, get, set, digits = 3, unit = 1, suffix = '', display, entry, enter,
  snap, play, reset, recompute = true, title }) {
  const after = recompute ? changed : () => {};
  const id = `sl${(slider.n = (slider.n || 0) + 1)}`;
  const row = document.createElement('div');
  row.className = 'sl';
  row.innerHTML = `<label for="${id}">${label}</label><div class="trk"><button type="button" class="tag" title="Click to type a value"></button>` +
    `<input type="range" class="rng" id="${id}"></div><div class="acts"></div>`;
  const rng = row.querySelector('input'), tag = row.querySelector('.tag'), trk = row.querySelector('.trk');
  const acts = row.querySelector('.acts');
  if (title) row.title = title;
  let min = lo, max = hi;
  Object.assign(rng, { min, max, step });
  const frac = (v) => Math.min(1, Math.max(0, (v - min) / (max - min)));
  const at = (f) => `calc(6.5px + (100% - 13px) * ${f.toFixed(5)})`;
  const text = (v) => (display ? display(v) : fixed(v / unit, digits) + suffix);
  let ticks = [];
  const refresh = () => {
    const v = get();
    rng.value = v;
    const f = frac(v);
    rng.style.setProperty('--p', `${(100 * f).toFixed(3)}%`);
    tag.textContent = text(v);
    tag.style.left = at(f);
    const marks = snap ? snap().filter((c) => c >= min && c <= max) : [];
    while (ticks.length > marks.length) ticks.pop().remove();
    while (ticks.length < marks.length) {
      const t = document.createElement('span');
      t.className = 'tick';
      trk.appendChild(t);
      ticks.push(t);
    }
    marks.forEach((c, k) => { ticks[k].style.left = at(frac(c)); });
  };
  rng.addEventListener('input', () => {
    let v = +rng.value;
    if (snap) {
      const px = (rng.clientWidth - 13) / (max - min);
      for (const c of snap()) if (Math.abs(v - c) * px < 3) v = c;
    }
    set(v);
    refresh();
    after(true);
  });
  tag.addEventListener('click', () => {
    const inp = document.createElement('input');
    inp.className = 'num tagin';
    inp.type = 'text';
    inp.inputMode = 'decimal';
    inp.autocomplete = 'off';
    inp.setAttribute('aria-label', `${rng.labels[0]?.textContent || 'value'}: exact value`);
    inp.value = entry ? entry() : fixed(get() / unit, digits);
    inp.style.left = tag.style.left;
    tag.hidden = true;
    trk.appendChild(inp);
    inp.focus();
    inp.select();
    let done = false;
    const finish = (commit) => {
      if (done) return;
      done = true;
      const v = parseNum(inp.value);
      inp.remove();
      tag.hidden = false;
      if (commit && v !== null) {
        if (enter) enter(v); else set(v * unit);
        after(false);
      }
      refresh();
    };
    inp.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') finish(true);
      else if (e.key === 'Escape') finish(false);
    });
    inp.addEventListener('blur', () => finish(true));
  });
  const button = (glyph, tip, onClick) => {
    const b = document.createElement('button');
    b.type = 'button';
    b.className = 'gly';
    b.innerHTML = GLYPH[glyph];
    b.title = tip;
    b.setAttribute('aria-label', tip);
    b.addEventListener('click', () => onClick(b));
    acts.appendChild(b);
    return b;
  };
  if (play) button('play', 'Animate', (b) => toggleAnim(play, b));
  if (reset) button('reset', reset.title, reset.onClick);
  parent.appendChild(row);
  refresh();
  const setRange = (a, b) => {
    if (a === min && b === max) return;
    [min, max] = [a, b];
    Object.assign(rng, { min, max });
    refresh();
  };
  return { refresh, row, setRange };
}

let tauSliders = [];
function buildTau() {
  const box = $('tauList');
  box.innerHTML = '';
  if (anim && anim.key === 'tau') stopAnim();
  tauSliders = state.tau.map((_, i) => slider(box, {
    label: `<i>τ</i><sub>${i + 1}</sub>`, min: -2 * PI, max: 2 * PI, step: 0.001,
    get: () => state.tau[i], set: (v) => { state.tau[i] = v; },
    play: { key: 'tau', i },
    reset: { title: `Set τ${i + 1} back to 0`, onClick: () => { state.tau[i] = 0; tauSliders[i].refresh(); changed(false); } },
  }));
  const g = genus();
  $('tauNote').textContent = g >= 3 ? ''
    : 'For genus ≤ 2 every isospectral deformation is a translation of the domain: use the domain centre below.';
}

const domainSliders = [
  slider($('domainRows'), { label: '<i>x</i><sub>0</sub>', title: 'Centre of the domain', min: -20, max: 20, step: 0.001, get: () => state.z0[0], set: (v) => { state.z0[0] = v; } }),
  slider($('domainRows'), { label: '<i>y</i><sub>0</sub>', title: 'Centre of the domain', min: -20, max: 20, step: 0.001, get: () => state.z0[1], set: (v) => { state.z0[1] = v; } }),
  slider($('domainRows'), { label: 'Width', min: 0.5, max: 40, step: 0.001, get: () => state.width, set: (v) => { state.width = Math.max(0.05, v); } }),
  slider($('domainRows'), { label: 'Height', min: 0.5, max: 40, step: 0.001, get: () => state.height, set: (v) => { state.height = Math.max(0.05, v); } }),
];
const phiSlider = slider($('phiRow'), { label: 'Angle', title: 'Angle of the grid in the z-plane', min: -PI, max: PI, step: 0.001, unit: PI, suffix: 'π', get: () => state.phi, set: (v) => { state.phi = v; } });
const appearanceSliders = [
  slider($('appearanceRows'), {
    label: 'Line spacing', min: 0.02, max: 2, step: 0.001, get: () => state.gridStep,
    set: (v) => { state.gridStep = Math.max(0.005, v); viewer.setStyle({ gridStep: state.gridStep }); },
    recompute: false,
  }),
  slider($('appearanceRows'), {
    label: 'Line width', min: 0.2, max: 2.5, step: 0.01, digits: 2, get: () => state.lineW,
    set: (v) => { state.lineW = v; viewer.setStyle({ lineWidth: v }); },
    recompute: false,
  }),
];

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
  label: '<i>H</i>', min: 0.05, max: 3, step: 0.01, digits: 2, get: () => state.H,
  set: (v) => {
    state.H = Math.max(0.01, v);
    viewer.setScale(1 / (2 * state.H));
    recolour();
    if (last) show(last);
  },
  recompute: false,
});

// ------------------------------------------------------------------ common roots of the differentials on S^1

// where every Theta_w vanishes: a Sym point there has no translational period in any direction
function updateCommonRoots() {
  let roots = [];
  if (genus() >= 1) {
    try { roots = commonRootsOnCircle(state.alphas, 1e-5).map((r) => r.lam); } catch { roots = []; }
  }
  widget.set({ commonRoots: roots });
}

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
  $('whithamNote').textContent = g < 2 ? 'Needs genus ≥ 2.' : '';
}

// W at the current spectral data, in the frame of the slider's centre (the curve's anchor)
function currentW() {
  if (genus() < 2) return null;
  const key = JSON.stringify([state.alphas, whitham.z]);
  if (currentW.key !== key) {
    currentW.key = key;
    try { currentW.val = willmore(state.alphas, whitham.z)[0]; } catch { currentW.val = null; }
  }
  return currentW.val;
}

// the family (Whitham curve through the current data), traced in its own worker so it never delays the
// surface; the latest request wins. It is traced in both directions until the flow meets an obstacle
// (a branch point reaching 0 or the unit circle, two colliding, or the continuation failing), and those
// ends are the ends of the slider. It is drawn in the λ-plane when the switch is on. Moving along the
// curve with the slider doesn't change the family, so only other edits of the spectral data (or
// recentring) trigger a new trace.
const family = { worker: null, busy: false, pending: null, key: null, data: null };
const WHITHAM_RANGE = [-1.5, 1.5]; // until the first family arrives
function showFamily() {
  widget.set({ family: $('showFamily').checked ? family.data : null });
}
function familyRange() {
  const P = family.data && family.data.points;
  return P && P.length > 1 ? [P[0].s, P[P.length - 1].s] : WHITHAM_RANGE;
}
function requestFamily() {
  if (genus() < 2) {
    family.pending = null;
    if (family.data) { family.data = null; family.key = null; showFamily(); }
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
      if (genus() >= 2) {
        family.data = error ? null : result;
        showFamily();
        whithamSlider.setRange(...familyRange());
        whithamSlider.refresh();
      }
      pumpFamily();
    };
  }
  const job = family.pending;
  family.pending = null;
  family.busy = true;
  family.worker.postMessage({ id: 0, kind: 'family', params: { alphas: job.alphas, smax: 30, maxSteps: 400 } });
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
// The critical points of W along the family (the ticks) are where B_a has a common root, a
// codimension-one condition, so the slider snaps onto them. Their positions come from the coarse
// family and are refined on the accurate curve (golden section on W), once per family.
function criticalPoints() {
  const F = family.data;
  if (!F || !whitham.curve || !F.critical.length) return [];
  if (whitham.critical && whitham.critical.family === F && whitham.critical.curve === whitham.curve) return whitham.critical.s;
  const W = (t) => whitham.curve.willmore(whitham.curve.at(t))[0];
  const list = F.critical.map((c) => {
    const P = F.points, sign = P[c].W > P[c - 1].W ? 1 : -1; // maximum or minimum
    let a = P[c - 1].s, b = P[c + 1].s;
    const gr = (Math.sqrt(5) - 1) / 2;
    let x1 = b - gr * (b - a), x2 = a + gr * (b - a), f1 = sign * W(x1), f2 = sign * W(x2);
    for (let it = 0; it < 40 && b - a > 1e-9; it++) {
      if (f1 > f2) { b = x2; x2 = x1; f2 = f1; x1 = b - gr * (b - a); f1 = sign * W(x1); }
      else { a = x1; x1 = x2; f1 = f2; x2 = a + gr * (b - a); f2 = sign * W(x2); }
    }
    return (a + b) / 2;
  });
  whitham.critical = { family: F, curve: whitham.curve, s: list };
  return list;
}
// the critical points for the slider's ticks: refined once there is a curve, else the coarse family's
function criticalS() {
  if (whitham.curve) return criticalPoints();
  const F = family.data;
  return F ? F.critical.map((c) => F.points[c].s) : [];
}

// exact entry of W: the nearest point (to the current one) of the curve where W takes that value,
// located on the family's samples (or the curve's, with the family off) and refined by secant steps
function whithamEnter(target) {
  if (genus() < 2) return;
  if (!whitham.curve) whithamMove(whitham.s);
  const C = whitham.curve;
  if (!C) return;
  const W = (t) => C.willmore(C.at(t))[0];
  const P = family.data ? family.data.points
    : Array.from({ length: 61 }, (_, k) => { const t = -1.5 + k * 0.05; return { s: t, W: W(t) }; }); // family not in yet
  let best = null;
  for (let i = 1; i < P.length; i++) {
    const a = P[i - 1].W - target, b = P[i].W - target;
    if (a * b > 0 || a === b) continue;
    const t = P[i - 1].s + (a / (a - b)) * (P[i].s - P[i - 1].s);
    if (best === null || Math.abs(t - whitham.s) < Math.abs(best - whitham.s)) best = t;
  }
  if (best === null) { // out of the family's range: go to the nearest value it takes
    best = P.reduce((m, p) => (Math.abs(p.W - target) < Math.abs(m.W - target) ? p : m)).s;
  } else {
    let t0 = best, t1 = best + 1e-3, f0 = W(t0) - target, f1 = W(t1) - target;
    for (let it = 0; it < 8 && Math.abs(f1) > 1e-10 && f1 !== f0; it++) {
      const t2 = t1 - (f1 * (t1 - t0)) / (f1 - f0);
      if (!Number.isFinite(t2) || Math.abs(t2 - best) > 0.2) break;
      [t0, f0, t1, f1] = [t1, f1, t2, W(t2) - target];
    }
    if (Math.abs(t1 - best) <= 0.2) best = t1;
  }
  const [a, b] = familyRange();
  whithamMove(Math.max(a, Math.min(b, best)));
}

const whithamSlider = slider($('whithamRow'), {
  label: '𝒲', min: -1.5, max: 1.5, step: 'any', get: () => whitham.s,
  set: (v) => { whithamMove(v); },
  display: () => { const W = currentW(); return W === null ? '–' : fixed(W, 4); },
  entry: () => { const W = currentW(); return W === null ? '' : fixed(W, 4); },
  enter: whithamEnter,
  snap: criticalS,
  play: { key: 'whitham' },
});
$('showFamily').addEventListener('change', showFamily);
$('whithamDomain').addEventListener('change', () => { if (whitham.curve) { whithamMove(whitham.s); changed(false); } });

// ------------------------------------------------------------------ closing up

$('closeBtn').addEventListener('click', () => {
  const w = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  const info = $('closeInfo');
  info.textContent = 'Searching for periods…';
  $('closeBtn').disabled = true;
  const done = () => { w.terminate(); $('closeBtn').disabled = false; };
  w.onerror = (e) => { info.textContent = `Error: ${e.message}`; done(); };
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
  btn.innerHTML = GLYPH.stop;
  animLast = performance.now();
  requestAnimationFrame(animStep);
}
function stopAnim() {
  if (animBtn) { animBtn.classList.remove('on'); animBtn.innerHTML = GLYPH.play; }
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
      const [a, b] = familyRange();
      const got = whithamMove(Math.max(a, Math.min(b, want)));
      if (Math.abs(got - want) > 1e-9 || got <= a || got >= b) anim.dir = -anim.dir;
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
presetSel.innerHTML = '<option value="">Choose a preset…</option>' +
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
  const url = linkURL();
  try {
    await navigator.clipboard.writeText(url);
    $('linkBtn').textContent = 'Copied';
    setTimeout(() => { $('linkBtn').textContent = 'Copy link'; }, 1500);
  } catch {
    window.prompt('Copy this link:', url); // no clipboard access (e.g. not https)
  }
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

// a link pasted into the address bar of an open page (clearing the hash with replaceState is silent)
window.addEventListener('hashchange', () => {
  const ok = readHash();
  clearHash();
  if (!ok) return;
  if (anim) stopAnim();
  syncTau();
  refreshAll();
  saveHomeDomain();
  needFrame = true;
  changed(false);
});

if (location.hash) { readHash(); clearHash(); }
syncTau();
makeWorker();
refreshAll();
saveHomeDomain();
setBackground(state.bg);
changed(false);
window.cmcStarted = true; // see index.html
