// Whitham deformations (Hasse–Ogilvie–Schmidt, Lemma "Whitham deformation"): move the branch points
// alpha_j so that the period plane BPer[B_a] with its complex structure J stays fixed.
//
// With the period map of periods.js, w -> (Im oint_{C_j} Theta_w / 2π)_j over the cycles C_j around
// {0, alpha_j} (the A-periods, around {alpha_j, 1/conj(alpha_j)}, vanish identically), write
//   L_alpha(w)_j = Re(ell_j w),   ell_j = v1_j - i v2_j,   v1 = L(1), v2 = L(i).
// The plane with its complex structure is the complex line [ell] in CP^{g-1}, so the leaf is
//   ell(alpha(t)) = z(t) ell(alpha(0)),  z(t) in C^*,
// which is 2g-2 real equations on the 2g real coordinates of alpha. The remaining freedom besides the
// Whitham flow is the rotation alpha -> e^{i phi} alpha; it is fixed by keeping arg a(0), i.e. the sum
// of arg alpha_j, constant (the normalisation a(0) = 1 of the paper). The lattice of periods of zeta
// then moves as Gamma(t) = Gamma(0) / z(t).

import { basicPeriods, periodVector } from './periods.js';

const M = 96; // quadrature points per cycle (the segment rule has converged far beyond 1e-12 here)

/** ell_j over the cycles around {0, alpha_j}, with signs aligned to ref (cycle orientations are arbitrary). */
export function ellVector(alphas, ref = null) {
  const P = basicPeriods(alphas, M);
  const v1 = periodVector(alphas, [1, 0], P), v2 = periodVector(alphas, [0, 1], P);
  return alphas.map((_, j) => {
    let re = v1[2 * j], im = -v2[2 * j];
    if (ref && re * ref[j][0] + im * ref[j][1] < 0) { re = -re; im = -im; }
    return [re, im];
  });
}

const toX = (alphas) => alphas.flat();
const toAlphas = (x) => Array.from({ length: x.length / 2 }, (_, j) => [x[2 * j], x[2 * j + 1]]);
const wrap = (t) => t - 2 * Math.PI * Math.round(t / (2 * Math.PI));
const cdiv = ([ar, ai], [br, bi]) => { const d = br * br + bi * bi; return [(ar * br + ai * bi) / d, (ai * br - ar * bi) / d]; };

/**
 * The Whitham curve through given spectral data, parametrised by (pseudo-)arclength s in alpha-space.
 * at(s) continues the curve from the anchor to s (caching what it has computed) and returns
 * { alphas, z, s } where s may fall short of the request if the curve leaves the allowed region.
 */
export class WhithamCurve {
  constructor(alphas, { h = 0.01 } = {}) {
    this.g = alphas.length;
    if (this.g < 2) throw new Error('Whitham deformations need genus >= 2');
    this.h = h;
    this.x0 = toX(alphas);
    this.arg0 = alphas.map(([re, im]) => Math.atan2(im, re));
    this.ell0 = ellVector(alphas);
    // divide by the largest component to form the ratios
    let k = 0;
    this.ell0.forEach((l, j) => { if (Math.hypot(...l) > Math.hypot(...this.ell0[k])) k = j; });
    this.k = k;
    this.r0 = this.ell0.map((l) => cdiv(l, this.ell0[k]));
    const t0 = this.nullVector(this.jacobian(this.x0, this.ell0));
    if (t0[0] < 0) for (let i = 0; i < t0.length; i++) t0[i] = -t0[i];
    const start = { s: 0, x: this.x0, t: t0, ell: this.ell0 };
    // both branches share the tangent; the sign of the step picks the direction
    this.branch = { 1: [start], [-1]: [{ ...start }] };
  }

  /** The defining equations G(x) = 0 (2g-1 of them) and the ell vector at x. */
  equations(x, ref) {
    const alphas = toAlphas(x);
    const ell = ellVector(alphas, ref);
    const G = [];
    ell.forEach((l, j) => {
      if (j === this.k) return;
      const r = cdiv(l, ell[this.k]);
      G.push(r[0] - this.r0[j][0], r[1] - this.r0[j][1]);
    });
    let gauge = 0;
    alphas.forEach(([re, im], j) => { gauge += wrap(Math.atan2(im, re) - this.arg0[j]); });
    G.push(gauge);
    return { G, ell };
  }

  jacobian(x, ref) {
    const { G } = this.equations(x, ref);
    const eps = 1e-6;
    const J = G.map(() => new Array(x.length).fill(0));
    for (let i = 0; i < x.length; i++) {
      const xp = x.slice(); xp[i] += eps;
      const xm = x.slice(); xm[i] -= eps;
      const Gp = this.equations(xp, ref).G, Gm = this.equations(xm, ref).G;
      for (let r = 0; r < G.length; r++) J[r][i] = (Gp[r] - Gm[r]) / (2 * eps);
    }
    return J;
  }

  /** Unit vector orthogonal to the rows of J (J has one row fewer than columns). */
  nullVector(J) {
    const n = J[0].length;
    const rows = [];
    for (const r of J) { // Gram–Schmidt on the rows
      const v = r.slice();
      for (const q of rows) { const d = dot(v, q); for (let i = 0; i < n; i++) v[i] -= d * q[i]; }
      const nv = Math.hypot(...v);
      if (nv > 1e-12) rows.push(v.map((c) => c / nv));
    }
    let best = null, bn = -1;
    for (let k = 0; k < n; k++) {
      const v = new Array(n).fill(0); v[k] = 1;
      for (const q of rows) { const d = dot(v, q); for (let i = 0; i < n; i++) v[i] -= d * q[i]; }
      const nv = Math.hypot(...v);
      if (nv > bn) { bn = nv; best = v.map((c) => c / nv); }
    }
    return best;
  }

  /** One predictor–corrector step of length h from point p along its tangent. Null if it fails. */
  step(p, h) {
    const n = p.x.length;
    const xp = p.x.map((v, i) => v + h * p.t[i]);
    let x = xp.slice(), ell = p.ell;
    for (let it = 0; it < 12; it++) {
      const eq = this.equations(x, ell);
      ell = eq.ell;
      const res = Math.hypot(...eq.G);
      const J = this.jacobian(x, ell);
      // augmented system: G = 0 and t . (x - xp) = 0
      const A = [...J, p.t.slice()], b = [...eq.G.map((v) => -v), -dot(p.t, x.map((v, i) => v - xp[i]))];
      const dx = solve(A, b, n);
      if (!dx.every(Number.isFinite)) return null;
      for (let i = 0; i < n; i++) x[i] += dx[i];
      if (res < 1e-12 && Math.hypot(...dx) < 1e-12) break;
      if (it === 11 && res > 1e-8) return null;
    }
    const alphas = toAlphas(x);
    const bad = alphas.some(([re, im]) => { const r = Math.hypot(re, im); return r < 0.02 || r > 0.98; })
      || alphas.some((a, i) => alphas.some((b, j) => j > i && Math.hypot(a[0] - b[0], a[1] - b[1]) < 2e-3));
    if (bad) return null;
    const eq = this.equations(x, ell);
    if (Math.hypot(...eq.G) > 1e-8) return null;
    let t = this.nullVector(this.jacobian(x, eq.ell));
    if (dot(t, p.t) < 0) t = t.map((v) => -v);
    return { s: p.s + Math.sign(h) * Math.hypot(...x.map((v, i) => v - p.x[i])), x, t, ell: eq.ell };
  }

  at(s) {
    const dir = s >= 0 ? 1 : -1;
    const pts = this.branch[dir];
    // continue until we pass s
    while (dir * pts[pts.length - 1].s < dir * s && !pts.done) {
      const q = this.step(pts[pts.length - 1], dir * this.h);
      if (!q) { pts.done = true; break; }
      pts.push(q);
    }
    // the last cached point not beyond s, then a partial step onto s
    let i = pts.length - 1;
    while (i > 0 && dir * pts[i].s > dir * s) i--;
    let p = pts[i];
    if (dir * p.s < dir * s) {
      const q = this.step(p, s - p.s);
      if (q) p = q;
    }
    const z = cdiv(p.ell[this.k], this.ell0[this.k]);
    return { alphas: toAlphas(p.x), z, s: p.s };
  }
}

const dot = (a, b) => a.reduce((acc, v, i) => acc + v * b[i], 0);

function solve(A, b, n) {
  const Mx = A.map((r, i) => [...r, b[i]]);
  for (let c = 0; c < n; c++) {
    let piv = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(Mx[r][c]) > Math.abs(Mx[piv][c])) piv = r;
    [Mx[c], Mx[piv]] = [Mx[piv], Mx[c]];
    for (let r = 0; r < n; r++) {
      if (r === c) continue;
      const f = Mx[r][c] / Mx[c][c];
      for (let k = c; k <= n; k++) Mx[r][k] -= f * Mx[c][k];
    }
  }
  return Mx.map((r, i) => r[n] / r[i]);
}
