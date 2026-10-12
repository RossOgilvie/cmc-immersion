// Tests for src/cmc against the Python reference values in test/fixtures.json and the checks of §8.
// Run from the repository root:  node test/run.mjs

import { readFileSync } from 'node:fs';
import {
  aPoly, kappa0, killingField, flowKilling, isospectralBasis, divisor, FrameIntegrator,
  computeSurface, hopfArg, closingInfo,
} from '../src/cmc/cmc.js';
import { periodVector, periodLattice, commonRootsOnCircle } from '../src/cmc/periods.js';
import { WhithamCurve } from '../src/cmc/whitham.js';
import { RootFlow, symData } from '../src/cmc/rootflow.js';
import { ellipticKE, nodalBoundary, f0, openEndTau } from '../src/cmc/genus2.js';
import {
  TauFlow, leafData, periodsOf, rationalPlanes, saturate, reduceFrame, classOrder, nearbyPeriods,
  PlaneFlow, latticeFrame, predictMoves,
} from '../src/cmc/leaf.js';

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
console.log('fixtures (JS vs the Python reference):');
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
  // the Willmore functional of the paper is real along the flow, and stationary where B_a has a
  // common root (the Wente torus has one, at lam0 = 1)
  const C = new WhithamCurve([[0.1412634686, 0.1017768953], [0.1412634686, -0.1017768953]]);
  let im = 0;
  for (const s of [-0.1, 0, 0.3]) { const W = C.willmore(C.at(s)); im = Math.max(im, Math.abs(W[1] / W[0])); }
  check('Willmore functional is real', im, 1e-10);
  const e = 1e-3, w = (s) => C.willmore(C.at(s))[0];
  check('dW/ds = 0 at the Wente torus (relative to the change over 0.2)',
    Math.abs((w(e) - w(-e)) / (2 * e)) / Math.abs((w(0.2) - w(0)) / 0.2), 1e-3);
}
{
  // common roots of B_a on S^1: the Wente torus has one at its Sym point lam0 = 1, generic data none
  const W = commonRootsOnCircle([[0.1412634686, 0.1017768953], [0.1412634686, -0.1017768953]], 1e-5);
  check('Wente: common root of the differentials at lam = 1', W.length === 1 ? Math.hypot(W[0].lam[0] - 1, W[0].lam[1]) : Infinity, 1e-9);
  check('generic genus 2: no common root', commonRootsOnCircle([[0.4, 0.25], [-0.3, 0.5]], 1e-5).length, 0);
}
{
  // genus 3 (generally no lattice): the curve continues and keeps [ell] in CP^2 fixed
  const C = new WhithamCurve([[0.4, 0.25], [-0.3, 0.5], [0.1, -0.6]]);
  const r = C.at(0.2), q = C.at(-0.2);
  check('g=3: Whitham curve continues to s = ±0.2', Math.max(Math.abs(r.s - 0.2), Math.abs(q.s + 0.2)), 1e-3);
}

// ------------------------------------------------------------------ root-preserving flow (CKKS), genus 2
{
  const PI = Math.PI, W = [[0.1412634686, 0.1017768953], [0.1412634686, -0.1017768953]];
  const phi = symData(W, 0).phi;
  check('Wente: Sym integrals phi = (pi/3, pi/3)', Math.hypot(phi[0] - PI / 3, phi[1] - PI / 3), 1e-7);
  const run = (f, target) => { let r; for (let i = 0; i < 500; i++) { r = f.moveToward(target, { budget: 50 }); if (r.reached || r.blocked) break; } return r; };
  const f = new RootFlow(W, 0);
  const r = run(f, [PI / 4, PI / 2]);
  const cr = commonRootsOnCircle(r.alphas, 1e-5);
  check('flow reaches phi = (pi/4, pi/2)', r.reached ? Math.hypot(r.phi[0] - PI / 4, r.phi[1] - PI / 2) : Infinity, 1e-9);
  check('flow keeps the common root at lam0 = 1', cr.length ? Math.min(...cr.map((c) => Math.hypot(c.lam[0] - 1, c.lam[1]))) : Infinity, 1e-8);
  run(f, [0.4 * PI, 0.2 * PI]);
  const back = run(f, [PI / 3, PI / 3]);
  check('flow around a loop returns to the Wente data (no holonomy)', Math.hypot(...back.alphas.flat().map((v, i) => v - W.flat()[i])), 1e-7);
  const end = run(new RootFlow(W, 0), [0.6 * PI, 0.6 * PI]);
  check('flow stops just short of the hypotenuse phi_1 + phi_2 = pi', end.blocked ? Math.abs(end.phi[0] + end.phi[1] - PI) / PI : Infinity, 5e-3);
}

// ------------------------------------------------------------------ root-preserving flow, genus 4
// Data on S^4 with real coefficients (the critical point of W on the Whitham curve through a symmetric
// start), lam0 = 1. The flow keeps the common root and the period plane; phi_2 = 0 by the symmetry.
{
  const PI = Math.PI;
  const A = [[0.241595657003, 0.131253423469], [0.241595657003, -0.131253423469], [-0.226370306867, 0.288849491705], [-0.226370306867, -0.288849491705]];
  const d0 = symData(A, 0);
  check('genus 4: phi_2 = 0 for real coefficients', Math.abs(d0.phi[1]), 1e-9);
  const P = d0.P.map((Pl, l) => (d0.phi[l] < 0 ? Pl.map((v) => -v) : Pl));
  const f = new RootFlow(A, 0, P);
  const target = [f.phi[0] + 0.004 * PI, f.phi[1] + 0.004 * PI];
  let r; for (let i = 0; i < 500; i++) { r = f.moveToward(target, { budget: 50, step: 0.005 }); if (r.reached || r.blocked) break; }
  check('genus 4: flow reaches a nearby target', r.reached ? Math.hypot(r.phi[0] - target[0], r.phi[1] - target[1]) : Infinity, 1e-9);
  const cr = commonRootsOnCircle(r.alphas, 1e-5);
  check('genus 4: flow keeps the common root at lam0 = 1', cr.length ? Math.min(...cr.map((c) => Math.hypot(c.lam[0] - 1, c.lam[1]))) : Infinity, 1e-8);
  const d = symData(r.alphas, 0, P, f.ref);
  const planeErr = Math.max(...[0, 1].flatMap((l) => { const c = d.Bm.map((row) => row[l]), n = Math.hypot(...c); return f.normals.map((N) => Math.abs(c.reduce((acc, v, j) => acc + v * N[j], 0) / n)); }));
  check('genus 4: flow keeps the period plane', planeErr, 1e-9);
}

// ------------------------------------------------------------------ leaves: the tau pad, LLL, locking periods
console.log('leaves of the Whitham foliation (leaf.js):');
{
  // tau at fixed psi = sum logit|alpha_j|: the periods of the frame and psi stay fixed, and a loop in tau
  // returns to the start
  const run = (f, t) => { let r; for (let i = 0; i < 2000; i++) { r = f.moveToward(t, { budget: 50 }); if (r.reached || r.blocked) break; } return r; };
  for (const al of [[[0.4, 0.25], [-0.3, 0.5]], [[0.4, 0.25], [-0.3, 0.5], [0.1, -0.6]]]) {
    const g = al.length, f = new TauFlow(al, [[1, 0], [0, 1]]);
    const r = run(f, [0.01, 1.01]);
    const d = leafData(r.alphas, f.ref);
    const perr = Math.max(...[0, 1].flatMap((k) => periodsOf(d.ell, r.w[k]).map((v, j) => Math.abs(v - f.P[k][j]))));
    check(`g=${g} tau pad reaches tau = 0.01 + 1.01 i`, r.reached ? Math.hypot(r.tau[0] - 0.01, r.tau[1] - 1.01) : Infinity, 1e-9);
    check(`g=${g} tau pad keeps the periods of the frame`, perr, 1e-9);
    const psi = r.alphas.reduce((s, [a, b]) => { const q = Math.hypot(a, b); return s + Math.log(q / (1 - q)); }, 0);
    check(`g=${g} tau pad keeps psi`, Math.abs(psi - f.psi0), 1e-9);
    run(f, [-0.01, 0.99]);
    const back = run(f, [0, 1]);
    check(`g=${g} tau pad around a loop returns to the start (no holonomy)`, Math.hypot(...back.alphas.flat().map((v, i) => v - al.flat()[i])), 1e-8);
  }
}
{
  // LLL finds a rational plane, with a basis of its lattice: span(m1, m2) with m1 = (1, 2, 0, -1), m2 = (0, 1, 3, 1)
  const m1 = [1, 2, 0, -1], m2 = [0, 1, 3, 1];
  const v1 = m1.map((v, j) => 0.7 * v + 0.2 * m2[j]), v2 = m1.map((v, j) => -0.3 * v + 1.1 * m2[j]);
  const P = rationalPlanes(v1.map((v, j) => [v, -v2[j]]));
  const best = P.find((p) => p.dist < 1e-10);
  // the found basis spans the same lattice: its coordinates in (m1, m2) form a unimodular matrix
  const coords = best ? best.m.map((m) => { const a = m[0], b = m[1] - 2 * a; return [a, b, m.every((v, j) => v === a * m1[j] + b * m2[j])]; }) : null;
  check('LLL recognises a rational plane with a basis of its lattice',
    coords && coords.every((c) => c[2]) ? Math.abs(Math.abs(coords[0][0] * coords[1][1] - coords[0][1] * coords[1][0]) - 1) : Infinity, 0);
  const S = saturate([2, 4, 0, -2], [1, 3, 3, 0]); // 2 m1 and m1 + m2: index 2 in the lattice
  const area = (a, b) => Math.sqrt(a.reduce((s, x) => s + x * x, 0) * b.reduce((s, x) => s + x * x, 0) - a.reduce((s, x, j) => s + x * b[j], 0) ** 2);
  check('saturate: the lattice of the plane (same area as m1, m2)', Math.abs(area(...S) - area(m1, m2)), 1e-9);
}
{
  // reduction into F_6, F_3, F_1 keeps the lattice and the classes mod 2
  let worst = 0;
  for (const kind of ['F6', 'F3', 'F1']) for (const t of [[3.7, 0.05], [-0.45, 0.3], [0.2, 0.01], [-7.1, 2]]) {
    const [b1, b2] = reduceFrame([{ w: [1, 0], m: [1, 0] }, { w: t, m: [0, 1] }], kind);
    const [x, y] = [(b2.w[0] * b1.w[0] + b2.w[1] * b1.w[1]) / (b1.w[0] ** 2 + b1.w[1] ** 2), (b2.w[1] * b1.w[0] - b2.w[0] * b1.w[1]) / (b1.w[0] ** 2 + b1.w[1] ** 2)];
    const inside = kind === 'F1' ? Math.abs(x) <= 0.5 + 1e-9 && x * x + y * y >= 1 - 1e-9
      : Math.abs(x) <= 1 + 1e-9 && Math.hypot(x - 0.5, y) >= 0.5 - 1e-9 && Math.hypot(x + 0.5, y) >= 0.5 - 1e-9 && (kind === 'F6' || x * x + y * y >= 1 - 1e-9);
    const det = b1.m[0] * b2.m[1] - b1.m[1] * b2.m[0];
    const parity = kind === 'F6' ? (b1.m[0] % 2 !== 0) && (b1.m[1] % 2 === 0) && (b2.m[0] % 2 === 0) && (b2.m[1] % 2 !== 0) : true;
    if (!inside || det !== 1 || !parity) worst = Infinity;
  }
  check('reduction into F_6 (by Gamma(2)), F_3 and F_1, with determinant 1', worst, 0);
  const labels = new Set([[0, 1.6], [0.78, 1.6], [0, 0.75], [0.3, 0.62], [0.75, 0.82], [0.8, 0.45]].map(classOrder));
  check('the six copies of F_1 in F_6 have six different class orders', 6 - labels.size, 0);
}
{
  // lock two integer periods into the plane (genus 3): the data become doubly periodic, and zeta closes up
  const al = [[0.4, 0.25], [-0.3, 0.5], [0.1, -0.6]];
  const go = (f) => { let s; for (let i = 0; i < 500; i++) { s = f.advance(100); if (s.done || s.blocked) break; } return s; };
  const d0 = leafData(al), pts = nearbyPeriods(d0.ell);
  const f1 = new PlaneFlow(al, [pts[0].m], d0.ell), s1 = go(f1);
  const d1 = f1.data(), second = nearbyPeriods(d1.ell).find((p) => p.dist > 1e-8);
  const pred = predictMoves(al, d0.ell, [pts[0].m])[0], moved = Math.hypot(...s1.alphas.flat().map((v, i) => v - al.flat()[i]));
  check('predicted move of a lock vs the actual move (relative)', Math.abs(pred - moved) / moved, 0.05);
  check('lock one period: it lies in the plane', s1.done ? nearbyPeriods(d1.ell).find((p) => p.m.join() === pts[0].m.join()).dist : Infinity, 1e-9);
  const f2 = new PlaneFlow(s1.alphas, [pts[0].m, second.m], d1.ell), s2 = go(f2);
  const d2 = f2.data(), L = latticeFrame(d2.ell, ...saturate(pts[0].m, second.m));
  const intErr = Math.max(...L.w.flatMap((w) => periodsOf(d2.ell, w).map((v) => Math.abs(v - Math.round(v)))));
  check('lock two periods: the lattice basis has integer periods', s2.done ? intErr : Infinity, 1e-9);
  const Z = killingField(s2.alphas), I = new FrameIntegrator(s2.alphas, 0);
  let zerr = 0;
  for (const w of L.w) {
    const y0 = I.initialState(Z), y = y0.slice(), r = Math.hypot(...w);
    I.march(y, w[0] / r, w[1] / r, r, 0.002);
    for (let k = 0; k < 8 * 5; k++) zerr = Math.max(zerr, Math.abs(y[k] - y0[k]));
  }
  check('lock two periods: zeta(z + w) = zeta(z) for both basis vectors (shooting)', zerr, 1e-8);
}

// ------------------------------------------------------------------ genus 2: the boundary of S^2 in F_3
console.log('genus 2: the boundary of S^2 (genus2.js) against the periods:');
{
  const { K, E } = ellipticKE([0.6, 0]);
  check('K(0.6), E(0.6) by the AGM', Math.max(Math.abs(K[0] - 1.7507538029157526), Math.abs(E[0] - 1.4180833944487245)), 1e-13);
  // nodal end: open the node at -e^{i phi}, t = cot(phi/2) = t0(k), slightly; simple root at -r
  const k = 0.5, { t0, tau } = nodalBoundary(k), r = (1 - k) / (1 + k), phi = 2 * Math.atan(1 / t0);
  const al = [[-(1 - 1e-4) * Math.cos(phi), -(1 - 1e-4) * Math.sin(phi)], [-r, 0]];
  const L = latticeFrame(leafData(al).ell, [1, 0], [0, 1]);
  const err = Math.min(Math.hypot(L.tau[0] - tau[0], L.tau[1] - tau[1]), Math.hypot(L.tau[0] + tau[0], L.tau[1] - tau[1]));
  check('nodal end: tau of the paper = tau of the lattice (node opened by 1e-4)', err, 1e-4);
  // open end: on the real locus of f0, rescaled roots s alpha, s / alpha with s = 1e-3
  const rho = 0.6, g = (th) => f0([rho * Math.cos(th), rho * Math.sin(th)])[1];
  let lo = 1e-6, hi = Math.PI / 2 - 1e-6;
  const glo = g(lo);
  for (let i = 0; i < 60; i++) { const m = (lo + hi) / 2; if (Math.sign(g(m)) === Math.sign(glo)) lo = m; else hi = m; }
  const a = [rho * Math.cos(lo), rho * Math.sin(lo)], s = 1e-3, a2 = a[0] ** 2 + a[1] ** 2;
  const d = leafData([[s * a[0], s * a[1]], [(s * a[0]) / a2, (-s * a[1]) / a2]]), M = latticeFrame(d.ell, [1, 0], [0, 1]);
  const t2 = openEndTau(a);
  const err2 = Math.min(Math.hypot(M.tau[0] - t2[0], M.tau[1] - t2[1]), Math.hypot(M.tau[0] + t2[0], M.tau[1] - t2[1]));
  check('open end: tau of the formula = tau of the lattice (s = 1e-3)', err2, 2e-3);
  const lin = (w) => [w[0] * d.p1[1][0] + w[1] * d.pi[1][0], w[0] * d.p1[1][1] + w[1] * d.pi[1][1]];
  const [u, v] = [lin(M.w[0]), lin(M.w[1])];
  check('open end: on the real locus of f0 the ratio of linear coefficients is real', Math.abs(u[0] * v[1] - u[1] * v[0]) / Math.hypot(...u) / Math.hypot(...v), 1e-3);
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
