// Tests for src/cmc against cmc_reference.py (via test/fixtures.json) and the checks of §8.
// Run from the repository root:  node test/run.mjs

import { readFileSync } from 'node:fs';
import {
  aPoly, kappa0, killingField, flowKilling, isospectralBasis, divisor, FrameIntegrator,
  computeSurface, hopfArg, closingInfo,
} from '../src/cmc/cmc.js';
import { periodVector, periodLattice } from '../src/cmc/periods.js';
import { WhithamCurve } from '../src/cmc/whitham.js';

let failures = 0;
function check(name, value, tol) {
  const ok = value <= tol;
  if (!ok) failures++;
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}: ${value.toExponential(2)} (tol ${tol.toExponential(0)})`);
}
const maxDiff = (a, b) => a.reduce((m, x, i) => Math.max(m, Math.abs(x - b[i])), 0);
const setDist = (A, B) => A.reduce((m, [ar, ai]) =>
  Math.max(m, Math.min(...B.map(([br, bi]) => Math.hypot(ar - br, ai - bi)))), 0);

// ------------------------------------------------------------------ fixtures
console.log('fixtures (JS vs cmc_reference.py):');
const fixtures = JSON.parse(readFileSync(new URL('./fixtures.json', import.meta.url)));
for (const fx of fixtures) {
  const g = fx.alphas.length;
  const tag = `g=${g}`;
  check(`${tag} a(lam)`, maxDiff(Array.from(aPoly(fx.alphas)), fx.a), 1e-14);
  const Z = killingField(fx.alphas, fx.A);
  check(`${tag} zeta base`, maxDiff(Array.from(Z), fx.Zbase), 1e-10);
  if (fx.times) {
    const kap = kappa0(fx.alphas);
    isospectralBasis(g).forEach(([k, c], i) => flowKilling(Z, kap, k, c, fx.times[i], 0.002));
    check(`${tag} zeta after isospectral flows`, maxDiff(Array.from(Z), fx.Zinit), 1e-9);
  }
  check(`${tag} divisor`, setDist(divisor(Z).map((d) => d.mu), fx.divisor), 1e-9);
  const I = new FrameIntegrator(fx.alphas, fx.theta0);
  for (const p of fx.points) {
    const y = I.initialState(Z);
    I.march(y, 1, 0, p.z[0], 0.005);
    I.march(y, 0, 1, p.z[1], 0.005);
    const P = new Float64Array(3), N = new Float64Array(3);
    const u = I.evaluate(y, P, 0, N, 0);
    const err = Math.max(maxDiff(Array.from(P), p.f), maxDiff(Array.from(N), p.N), Math.abs(u - p.u));
    check(`${tag} f, N, u at z=${p.z}`, err, 1e-9);
  }
}

// ------------------------------------------------------------------ geometric invariants on the grid
console.log('finite-difference invariants on computeSurface (errors are O(h^2)):');
function geomCheck(label, params) {
  const S = computeSurface(params);
  const { nx, ny, pos, nrm, u } = S;
  const hs = params.width / (nx - 1), ht = params.height / (ny - 1);
  const P = (i, j, e) => pos[3 * (j * nx + i) + e];
  let conf = 0, orth = 0, metric = 0, Herr = 0, Qerr = 0, gauss = 0;
  // expected Hopf differential in the rotated coordinate w = s + i t: Q_w = e^{2 i phi} Q_z
  const argQ = hopfArg(params.alphas, params.theta0) + 2 * S.phi;
  for (let j = 2; j < ny - 2; j++) {
    for (let i = 2; i < nx - 2; i++) {
      const fs = [], ft = [], fss = [], ftt = [], fst = [], Nn = [];
      for (let e = 0; e < 3; e++) {
        fs[e] = (P(i + 1, j, e) - P(i - 1, j, e)) / (2 * hs);
        ft[e] = (P(i, j + 1, e) - P(i, j - 1, e)) / (2 * ht);
        fss[e] = (P(i + 1, j, e) - 2 * P(i, j, e) + P(i - 1, j, e)) / (hs * hs);
        ftt[e] = (P(i, j + 1, e) - 2 * P(i, j, e) + P(i, j - 1, e)) / (ht * ht);
        fst[e] = (P(i + 1, j + 1, e) - P(i + 1, j - 1, e) - P(i - 1, j + 1, e) + P(i - 1, j - 1, e)) / (4 * hs * ht);
        Nn[e] = nrm[3 * (j * nx + i) + e];
      }
      const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
      const E = dot(fs, fs), G = dot(ft, ft), F = dot(fs, ft);
      const uu = u[j * nx + i];
      conf = Math.max(conf, Math.abs(E - G) / E);
      orth = Math.max(orth, Math.abs(F) / E);
      metric = Math.max(metric, Math.abs(E / (4 * Math.exp(2 * uu)) - 1));
      const H = (dot(fss, Nn) + dot(ftt, Nn)) / (2 * E);
      Herr = Math.max(Herr, Math.abs(H - 0.5));
      const Qr = 0.25 * (dot(fss, Nn) - dot(ftt, Nn)), Qi = -0.5 * dot(fst, Nn);
      Qerr = Math.max(Qerr, Math.hypot(Qr - Math.cos(argQ), Qi - Math.sin(argQ)));
      // sinh-Gordon: Delta u + 2 sinh 2u = 0
      const lap = (u[j * nx + i + 1] - 2 * uu + u[j * nx + i - 1]) / (hs * hs)
        + (u[(j + 1) * nx + i] - 2 * uu + u[(j - 1) * nx + i]) / (ht * ht);
      gauss = Math.max(gauss, Math.abs(lap + 2 * Math.sinh(2 * uu)) / (1 + Math.abs(2 * Math.sinh(2 * uu))));
    }
  }
  const tol = 3e-3;
  check(`${label}: |E-G|/E`, conf, tol);
  check(`${label}: |F|/E`, orth, tol);
  check(`${label}: E/(4e^2u) - 1`, metric, tol);
  check(`${label}: H - 1/2`, Herr, tol);
  check(`${label}: Q error`, Qerr, tol);
  check(`${label}: sinh-Gordon residual`, gauss, 2e-2);
  check(`${label}: det zeta drift`, S.stats.detDefect, 1e-8);
}
const base = { width: 1.5, height: 1.5, nx: 61, ny: 61, hmax: 0.01 };
geomCheck('g=0', { ...base, alphas: [], theta0: 0.3 });
geomCheck('g=1', { ...base, alphas: [[0.4, 0.2]], theta0: 0.7, z0: [0.3, -0.4] });
geomCheck('g=2', { ...base, alphas: [[0.5, 0.1], [-0.3, 0.6]], theta0: 1.3, phi: 0.2 });
geomCheck('g=3', { ...base, alphas: [[0.3, 0], [0, 0.5], [-0.6, -0.2]], theta0: Math.PI / 2, tau: [0.4] });
geomCheck('g=4', { ...base, alphas: [[0.6, 0], [0, 0.4], [-0.5, -0.3], [0.2, -0.7]], theta0: 2, tau: [0.3, -0.5] });

// ------------------------------------------------------------------ Delaunay radii
console.log('Delaunay surfaces (g=1): neck and bulge radii of the parallels:');
function circumradius(a, b, c) {
  const d = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
  const A = d(b, c), B = d(a, c), C = d(a, b);
  const s = (A + B + C) / 2;
  return (A * B * C) / (4 * Math.sqrt(s * (s - A) * (s - B) * (s - C)));
}
for (const r of [0.3, 0.7]) {
  for (const sign of [-1, 1]) {
    const th = 0.9;
    const S = computeSurface({
      alphas: [[r * Math.cos(th), r * Math.sin(th)]], theta0: th + (sign < 0 ? Math.PI : 0),
      phi: 0, width: 1.4, height: 12, nx: 3, ny: 1201, hmax: 0.01,
    });
    let lo = Infinity, hi = 0;
    for (let j = 0; j < S.ny; j++) {
      const p = (i) => Array.from(S.pos.subarray(3 * (j * 3 + i), 3 * (j * 3 + i) + 3));
      const R = circumradius(p(0), p(1), p(2));
      lo = Math.min(lo, R); hi = Math.max(hi, R);
    }
    const exp = sign < 0 ? [2 * r / (1 + r), 2 / (1 + r)] : [2 * r / (1 - r), 2 / (1 - r)];
    check(`|alpha|=${r} ${sign < 0 ? 'unduloid' : 'nodoid'} radii [${lo.toFixed(4)}, ${hi.toFixed(4)}]`,
      Math.max(Math.abs(lo - exp[0]), Math.abs(hi - exp[1])), 1e-4);
  }
}
{
  // curvature-line grid: one family of parameter lines is the unit circles, the other straight lines
  const S = computeSurface({ alphas: [], theta0: 1.1, width: 2, height: 2, nx: 3, ny: 3, hmax: 0.01 });
  const p = (i, j) => Array.from(S.pos.subarray(3 * (3 * j + i), 3 * (3 * j + i) + 3));
  const Rs = circumradius(p(0, 1), p(1, 1), p(2, 1)), Rt = circumradius(p(1, 0), p(1, 1), p(1, 2));
  check('g=0 cylinder radius 1', Math.min(Math.abs(Rs - 1), Math.abs(Rt - 1)), 1e-9);
}

// ------------------------------------------------------------------ commutativity and divisor
console.log('commutativity around a closed rectangle, and base-point divisor:');
{
  const alphas = [[0.5, 0.1], [-0.3, 0.6], [0.1, -0.4]];
  const I = new FrameIntegrator(alphas, 0.8);
  const Z = killingField(alphas);
  flowKilling(Z, kappa0(alphas), 1, [0, 1], 0.7);
  const y0 = I.initialState(Z);
  const y = y0.slice();
  I.march(y, 1, 0, 1.3, 0.005); I.march(y, 0, 1, 0.9, 0.005);
  I.march(y, 1, 0, -1.3, 0.005); I.march(y, 0, 1, -0.9, 0.005);
  check('g=3 rectangle holonomy of (zeta, F, G)', maxDiff(Array.from(y), Array.from(y0)), 1e-9);
  const alphas2 = [[0.3, 0], [0, 0.5], [-0.2, -0.6]];
  check('base-point divisor = branch points', setDist(divisor(killingField(alphas2)).map((d) => d.mu), alphas2), 1e-12);
}

// ------------------------------------------------------------------ periods and closing (§9)
console.log('periods of Theta_w and closing conditions:');
{
  // g = 1: the y-period of the unduloid, against shooting (closingInfo)
  const T = Math.abs(1 / periodVector([[0.45, 0]], [0, 1])[0]);
  check('g=1 y-period from Theta_w vs shooting', Math.abs(T - closingInfo({ alphas: [[0.45, 0]], theta0: Math.PI }, 't').T), 1e-6);
  // g = 2: every lattice generator is a period of zeta (§9.4)
  const al = [[0.3, 0.4], [0.3, -0.4]];
  const { generators } = periodLattice(al);
  let worst = 0;
  for (const w of generators) {
    const I = new FrameIntegrator(al, 0), Z = killingField(al), y = I.initialState(Z);
    const L = Math.hypot(...w);
    I.march(y, w[0] / L, w[1] / L, L, 0.002);
    for (let e = 0; e < 32; e++) worst = Math.max(worst, Math.abs(y[e] - Z[e]));
  }
  check(`g=2 lattice generators (${generators.length}) are periods of zeta`, generators.length === 2 ? worst : Infinity, 1e-9);
  // Wente torus: closes around 3 w1 and 2 w2 - w1
  const W = [[0.1412634686, 0.1017768953], [0.1412634686, -0.1017768953]];
  // (generator signs are arbitrary: use T1 = |w1| on the real axis and T2/2 = |Im w2|)
  const [w1, w2] = periodLattice(W).generators;
  const T1 = Math.abs(w1[0]), T2h = Math.abs(w2[1]);
  let err = 0;
  for (const w of [[3 * T1, 0], [1.5 * T1, T2h], [0, 2 * T2h]]) {
    const I = new FrameIntegrator(W, 0), Z = killingField(W), y = I.initialState(Z);
    const L = Math.hypot(...w);
    I.march(y, w[0] / L, w[1] / L, L, 0.002);
    const P = new Float64Array(3), N = new Float64Array(3);
    I.evaluate(y, P, 0, N, 0);
    for (let e = 0; e < 32; e++) err = Math.max(err, Math.abs(y[e] - Z[e]));
    err = Math.max(err, Math.hypot(...P), Math.abs(Math.abs(y[32]) - 1));
  }
  check('Wente torus closes on the brick lattice <3 T1, 3 T1/2 + i T2/2>', err, 1e-7);
  // bubbleton: closes around the cylinder at the resonance point
  const r2 = (2 - Math.sqrt(3)) ** 2;
  const c = closingInfo({ alphas: [[r2, 0], [r2, 0]], theta0: 0, hmax: 0.002 }, 's');
  check('two-lobed bubbleton closes after one period pi', c && c.q === 1 ? Math.abs(c.T - Math.PI) : Infinity, 1e-5);
}

// ------------------------------------------------------------------ Whitham deformation
console.log('Whitham deformation keeps the period lattice up to Gamma(0)/z:');
for (const [name, al] of [
  ['Wente', [[0.1412634686, 0.1017768953], [0.1412634686, -0.1017768953]]],
  ['generic g=2', [[0.4, 0.25], [-0.3, 0.5]]],
]) {
  const C = new WhithamCurve(al);
  const L0 = periodLattice(al).generators;
  const cov = (L) => { // lattice invariants: |w1|, |w2|, |w1 x w2|, up to the choice of generators
    const [a, b] = L;
    return [Math.hypot(...a), Math.hypot(...b), Math.abs(a[0] * b[1] - a[1] * b[0])];
  };
  let worst = 0, ok = true;
  for (const s of [-0.1, 0.1, 0.25]) {
    const r = C.at(s);
    if (Math.abs(r.s - s) > 1e-3) { ok = false; continue; } // s is the arclength actually travelled
    // Gamma(s) z = Gamma(0): compare the invariants of Gamma(s) scaled by |z| (and |z|^2 for the area)
    const L = periodLattice(r.alphas).generators;
    if (L.length < 2 || L0.length < 2) { ok = false; continue; }
    const m = Math.hypot(...r.z);
    const [a, b, A] = cov(L), [a0, b0, A0] = cov(L0);
    worst = Math.max(worst, Math.abs(a * m - a0) / a0, Math.abs(b * m - b0) / b0, Math.abs(A * m * m - A0) / A0);
  }
  check(`${name}: lattice of the deformed curve = Gamma(0)/z`, ok ? worst : Infinity, 1e-6);
}
{
  // genus 3 (generally no lattice): the curve continues and keeps [ell] in CP^2 fixed
  const C = new WhithamCurve([[0.4, 0.25], [-0.3, 0.5], [0.1, -0.6]]);
  const r = C.at(0.2), q = C.at(-0.2);
  check('g=3: Whitham curve continues to s = ±0.2', Math.max(Math.abs(r.s - 0.2), Math.abs(q.s + 0.2)), 1e-3);
}

// ------------------------------------------------------------------ timing
console.log('timing:');
for (const [g, n] of [[2, 64], [2, 256], [4, 256]]) {
  const alphas = [[0.6, 0], [0, 0.4], [-0.5, -0.3], [0.2, -0.7]].slice(0, g);
  const S = computeSurface({ alphas, theta0: 1, tau: [0.3, 0.2], width: 20, height: 20, nx: n, ny: n, hmax: 0.05 });
  console.log(`  g=${g} ${n}x${n} over 20x20, hmax=0.05: ${S.stats.ms.toFixed(0)} ms`);
}

console.log(failures ? `\n${failures} FAILED` : '\nall passed');
process.exit(failures ? 1 : 0);
