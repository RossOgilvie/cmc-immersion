// Tests for src/cmc against cmc_reference.py (via test/fixtures.json) and the checks of §8.
// Run from the repository root:  node test/run.mjs

import { readFileSync } from 'node:fs';
import {
  aPoly, kappa0, killingField, flowKilling, isospectralBasis, divisor, FrameIntegrator,
  computeSurface, hopfArg,
} from '../src/cmc/cmc.js';

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

// ------------------------------------------------------------------ timing
console.log('timing:');
for (const [g, n] of [[2, 64], [2, 256], [4, 256]]) {
  const alphas = [[0.6, 0], [0, 0.4], [-0.5, -0.3], [0.2, -0.7]].slice(0, g);
  const S = computeSurface({ alphas, theta0: 1, tau: [0.3, 0.2], width: 20, height: 20, nx: n, ny: n, hmax: 0.05 });
  console.log(`  g=${g} ${n}x${n} over 20x20, hmax=0.05: ${S.stats.ms.toFixed(0)} ms`);
}

console.log(failures ? `\n${failures} FAILED` : '\nall passed');
process.exit(failures ? 1 : 0);
