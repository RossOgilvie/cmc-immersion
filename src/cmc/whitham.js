// Whitham deformations (Hasse–Ogilvie–Schmidt, Lemma "Whitham deformation"): move the branch points
// alpha_j so that the period plane BPer[B_a] with its complex structure J stays fixed.
//
// With the period map of periods.js, w -> (Im oint_{C_j} Theta_w / 2π)_j over the cycles C_j around
// {0, alpha_j} (the A-periods, around {alpha_j, 1/conj(alpha_j)}, vanish identically), write
//   L_alpha(w)_j = Re(ell_j w),   ell_j = v1_j - i v2_j,   v1 = L(1), v2 = L(i).
// The plane with its complex structure is the complex line [ell] in CP^{g-1}, so the leaf is
//   ell(alpha(t)) = z(t) ell(alpha(0)),  z(t) in C^*.
// We solve this with z as an unknown (no division by a component of ell), plus the gauge
// sum_j arg alpha_j = const, which fixes the rotation alpha -> e^{i phi} alpha (the paper's
// normalisation a(0) = 1): 2g+1 real equations in the 2g+2 unknowns (alpha, z), a curve.
// The lattice of periods of zeta moves as Gamma(t) = Gamma(0) / z(t).
//
// Continuation: predictor along the tangent, chord Newton on the pseudo-arclength system with the
// Jacobian of the step's start point (forward differences in alpha; the z-columns are exact).
// The parameter s is pseudo-arclength in alpha-space (arclength to O(h^3) per step).
//
// The Willmore functional of the paper, W = 2i Res_{lam=0} q_{b2} Theta_{b1} for the parallel frame
// b_k = p_{w_k / z} (w_1 = 1, w_2 = i at the anchor), becomes with lam = kappa^2 and
// Theta_b = 2 b(kappa^2) / (kappa^2 sqrt(a(kappa^2))) dkappa:
//   W = 8i (b_0^{(1)} b_1^{(2)} - b_0^{(2)} b_1^{(1)}) / a_0,
// where b_0, b_1 are the coefficients of 1 and lam (the a_1 terms cancel). It is real, and its
// critical points along the curve are where B_a has a common root.

import { differentials } from './periods.js';
import { aPoly } from './cmc.js';

const M = 96; // minimum quadrature points per cycle (basicPeriods refines near other branch points)
// the imaginary part of the period of p over a B-cycle, divided by 2 pi
const imPeriod = (p, Pj) => p.reduce((acc, c, k) => acc + c[0] * Pj[k][1] + c[1] * Pj[k][0], 0) / (2 * Math.PI);

/** ell_j over the cycles around {0, alpha_j}, with signs aligned to ref (cycle orientations are arbitrary). */
export function ellVector(alphas, ref = null) {
  // only the B-cycles (around {0, alpha_j}); the A-periods vanish by the reality condition
  const D = differentials(alphas, M);
  const p1 = D.poly([1, 0]), p2 = D.poly([0, 1]);
  return alphas.map((_, j) => {
    let re = imPeriod(p1, D.B[j]), im = -imPeriod(p2, D.B[j]);
    if (ref && re * ref[j][0] + im * ref[j][1] < 0) { re = -re; im = -im; }
    return [re, im];
  });
}

const toAlphas = (x) => Array.from({ length: x.length / 2 }, (_, j) => [x[2 * j], x[2 * j + 1]]);
const wrap = (t) => t - 2 * Math.PI * Math.round(t / (2 * Math.PI));
const cmul = ([ar, ai], [br, bi]) => [ar * br - ai * bi, ar * bi + ai * br];
const cdiv = ([ar, ai], [br, bi]) => { const d = br * br + bi * bi; return [(ar * br + ai * bi) / d, (ai * br - ar * bi) / d]; };
const dot = (a, b) => a.reduce((acc, v, i) => acc + v * b[i], 0);

/**
 * The Whitham curve through given spectral data.
 * at(s) continues the curve from the anchor to s (caching) and returns { alphas, z, s }, with s short
 * of the request if the curve leaves the allowed region. trace() returns both branches for drawing.
 */
export class WhithamCurve {
  constructor(alphas, { h = 0.05 } = {}) {
    this.g = alphas.length;
    if (this.g < 2) throw new Error('Whitham deformations need genus >= 2');
    this.h = h;
    this.arg0 = alphas.map(([re, im]) => Math.atan2(im, re));
    this.ell0 = ellVector(alphas);
    const y0 = [...alphas.flat(), 1, 0];
    const J = this.jacobian(y0, this.ell0);
    let t = this.tangent(J);
    if (t[0] < 0) t = t.map((v) => -v);
    const start = { s: 0, y: y0, t, J, ell: this.ell0 };
    // both branches share the tangent; the sign of the step picks the direction
    this.branch = { 1: [start], [-1]: [{ ...start }] };
  }

  /** G(y) for y = (alpha as reals, zr, zi): the 2g real parts of ell - z ell0, then the gauge. */
  equations(y, ref) {
    const n = 2 * this.g;
    const alphas = toAlphas(y.slice(0, n));
    const ell = ellVector(alphas, ref);
    const z = [y[n], y[n + 1]];
    const G = [];
    ell.forEach((l, j) => { const zl = cmul(z, this.ell0[j]); G.push(l[0] - zl[0], l[1] - zl[1]); });
    let gauge = 0;
    alphas.forEach(([re, im], j) => { gauge += wrap(Math.atan2(im, re) - this.arg0[j]); });
    G.push(gauge);
    return { G, ell };
  }

  /** Jacobian of G at y: forward differences in alpha, exact in z. */
  jacobian(y, ref) {
    const n = 2 * this.g;
    const G0 = this.equations(y, ref).G;
    const J = G0.map(() => new Array(n + 2).fill(0));
    const eps = 1e-7;
    for (let i = 0; i < n; i++) {
      const yp = y.slice(); yp[i] += eps;
      const Gp = this.equations(yp, ref).G;
      for (let r = 0; r < G0.length; r++) J[r][i] = (Gp[r] - G0[r]) / eps;
    }
    this.ell0.forEach((l, j) => {
      J[2 * j][n] = -l[0]; J[2 * j + 1][n] = -l[1];         // d/dzr of -(z ell0)
      J[2 * j][n + 1] = l[1]; J[2 * j + 1][n + 1] = -l[0];  // d/dzi
    });
    return J;
  }

  /** Null vector of J (one row fewer than columns), scaled to unit length in alpha. */
  tangent(J) {
    const m = J[0].length, n = 2 * this.g;
    const rows = [];
    for (const r of J) { // Gram–Schmidt on the rows
      const v = r.slice();
      for (const q of rows) { const d = dot(v, q); for (let i = 0; i < m; i++) v[i] -= d * q[i]; }
      const nv = Math.hypot(...v);
      if (nv > 1e-12) rows.push(v.map((c) => c / nv));
    }
    let best = null, bn = -1;
    for (let k = 0; k < m; k++) {
      const v = new Array(m).fill(0); v[k] = 1;
      for (const q of rows) { const d = dot(v, q); for (let i = 0; i < m; i++) v[i] -= d * q[i]; }
      const nv = Math.hypot(...v);
      if (nv > bn) { bn = nv; best = v; }
    }
    const na = Math.hypot(...best.slice(0, n));
    return best.map((c) => c / na);
  }

  /** One predictor–corrector step of alpha-arclength h from point p, to residual tol. Null if it fails. */
  step(p, h, tol = 1e-11) {
    const n = 2 * this.g, m = n + 2;
    const yp = p.y.map((v, i) => v + h * p.t[i]);
    const y = yp.slice();
    let ell = p.ell, res = Infinity;
    // pseudo-arclength in alpha: the correction stays in the hyperplane t_alpha . (y - yp) = 0,
    // so s advances by exactly h (= alpha-arclength to O(h^3))
    const ta = p.t.map((v, i) => (i < n ? v : 0));
    const A = [...p.J, ta]; // chord Newton: the start point's Jacobian throughout
    for (let it = 0; it < 30; it++) {
      const eq = this.equations(y, ell);
      ell = eq.ell;
      res = Math.hypot(...eq.G);
      const b = [...eq.G.map((v) => -v), -dot(ta, y.map((v, i) => v - yp[i]))];
      const dy = solve(A, b, m);
      if (!dy.every(Number.isFinite)) return null;
      for (let i = 0; i < m; i++) y[i] += dy[i];
      // an iterate outside the disc can't converge to an admissible point (and its periods are costly)
      for (let j = 0; j < n; j += 2) { const r = Math.hypot(y[j], y[j + 1]); if (!(r > 1e-4 && r < 0.9999)) return null; }
      if (Math.hypot(...dy) < tol && res < tol) break;
    }
    if (!(res < Math.max(1e-9, 10 * tol))) return null;
    const alphas = toAlphas(y.slice(0, n));
    const bad = alphas.some(([re, im]) => { const r = Math.hypot(re, im); return r < 1e-3 || r > 0.999; })
      || alphas.some((a, i) => alphas.some((c, j) => j > i && Math.hypot(a[0] - c[0], a[1] - c[1]) < 2e-3));
    if (bad) return null;
    const J = this.jacobian(y, ell);
    let t = this.tangent(J);
    if (dot(t, p.t) < 0) t = t.map((v) => -v);
    return { s: p.s + h, y, t, J, ell };
  }

  /**
   * Extend branch dir (±1) until |s| >= smax or maxSteps; step length adapts on failure. Points are
   * computed to residual tol; a coarse branch (for drawing) is kept apart from the accurate one.
   */
  extend(dir, smax, maxSteps = Infinity, tol = 1e-11) {
    const key = tol < 1e-10 ? dir : `${dir}coarse`;
    if (!this.branch[key]) this.branch[key] = [{ ...this.branch[1][0] }];
    const pts = this.branch[key];
    let h = this.h, steps = 0;
    while (!pts.done && dir * pts[pts.length - 1].s < smax && steps < maxSteps) {
      const q = this.step(pts[pts.length - 1], dir * h, tol);
      if (q) { pts.push(q); steps++; h = Math.min(this.h, h * 1.5); }
      else if (h > this.h / 64) h /= 4;
      else pts.done = true;
    }
    return pts;
  }

  at(s) {
    const dir = s >= 0 ? 1 : -1;
    const pts = this.extend(dir, dir * s);
    let i = pts.length - 1;
    while (i > 0 && dir * pts[i].s > dir * s) i--;
    let p = pts[i];
    if (dir * p.s < dir * s) {
      const q = this.step(p, s - p.s);
      if (q) p = q;
    }
    return this.describe(p);
  }

  describe(p) {
    const n = 2 * this.g;
    return { alphas: toAlphas(p.y.slice(0, n)), z: [p.y[n], p.y[n + 1]], s: p.s };
  }

  /** The Willmore functional at a point of the curve, in the frame parallel to (1, i) at the anchor. */
  willmore({ alphas, z }) {
    return willmore(alphas, z);
  }

  /**
   * Both branches up to |s| <= smax (at most maxSteps each), as [{ s, alphas, z, W }] sorted by s,
   * plus the critical points of W (sign changes of its differences).
   */
  trace(smax = 1.5, maxSteps = 80, tol = 1e-8) {
    const neg = this.extend(-1, smax, maxSteps, tol), pos = this.extend(1, smax, maxSteps, tol);
    const pts = [...neg.slice(1).reverse(), ...pos].map((p) => {
      const d = this.describe(p);
      return { ...d, W: this.willmore(d)[0] };
    });
    const critical = [];
    for (let i = 1; i + 1 < pts.length; i++) {
      const a = pts[i].W - pts[i - 1].W, b = pts[i + 1].W - pts[i].W;
      if (a * b < 0) critical.push(i);
    }
    return { points: pts, critical, ends: { neg: !!neg.done, pos: !!pos.done } };
  }
}

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

/**
 * The Willmore functional W = 8i (b0^1 b1^2 - b0^2 b1^1) / a0 of the spectral data, with b^k the
 * differentials p_{w_k / z}: z = 1 is the frame parallel to (1, i) here, and a point of a Whitham curve
 * passes its own z to stay in the frame of the curve's anchor. Returns [re, im]; im vanishes up to
 * quadrature error.
 */
export function willmore(alphas, z = [1, 0]) {
  const D = differentials(alphas, M);
  const b1 = D.poly(cdiv([1, 0], z)), b2 = D.poly(cdiv([0, 1], z));
  const a = aPoly(alphas);
  const d = [cmul(b1[0], b2[1])[0] - cmul(b2[0], b1[1])[0], cmul(b1[0], b2[1])[1] - cmul(b2[0], b1[1])[1]];
  return cmul([0, 8], cdiv(d, [a[0], a[1]]));
}
