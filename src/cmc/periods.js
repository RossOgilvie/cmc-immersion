// Periods of the differentials Theta_w and the lattice of periods of zeta (§9 of
// cmc_spectral_to_immersion.md).
//
// Theta_w = p_w(lam) dlam / (lam nu),  nu^2 = lam a(lam),  deg p_w <= g+1, with
//   p_w(0) = -w a_0 / (2 kappa0),   [lam^{g+1}] p_w = conj(w) a_{2g} / (2 kappa0),
// and the g middle coefficients fixed by Re oint_gamma Theta_w = 0 on all cycles.
// w is a period of zeta iff oint_gamma Theta_w lies in 2 pi i Z for all cycles.

import { polyEval } from './poly.js';
import { aPoly, kappa0 } from './cmc.js';

/**
 * Cycles: loops around the pairs of branch points {0, alpha_j} and {alpha_j, 1/conj(alpha_j)}.
 * Each loop is collapsed onto the segment [e1, e2] and parametrised by lam = m - h cos(phi), which
 * cancels the square-root singularities: sqrt((lam-e1)(lam-e2)) = ±i h sin(phi) against
 * dlam = h sin(phi) dphi. The resulting integrand is smooth, even and 2π-periodic in phi, so the
 * trapezoid rule converges geometrically. (Periods are defined up to sign, which the conditions
 * below don't see.)
 *
 * The loops around {0, alpha_j} can't be collapsed for dlam/(lam nu) ~ lam^{-3/2}, so that one is
 * replaced modulo the exact form d(nu/lam) = (lam a' - a)/(2 lam nu) dlam:
 *   dlam/(lam nu) = (1/a0) (a'(lam) - (a(lam) - a0)/lam) dlam/nu - (2/a0) d(nu/lam).
 * Returns, for each cycle, the periods of lam^k dlam/(lam nu) for k = 0..g+1 as [re, im] pairs.
 */
export function basicPeriods(alphas, M = 256) {
  const g = alphas.length;
  const a = aPoly(alphas);
  // P(lam) = a'(lam) - (a(lam) - a0)/lam, coefficients by increasing power (degree <= 2g - 1)
  const P = [];
  for (let m = 0; m < 2 * g; m++) P.push([(m + 1) * a[2 * m + 2] - a[2 * m + 2], (m + 1) * a[2 * m + 3] - a[2 * m + 3]]);
  const a0r = a[0], a0i = a[1], a02 = a0r * a0r + a0i * a0i;
  // the factors of lam a(lam) = -lam prod (lam - alpha_j)(1 - conj(alpha_j) lam), one per branch point:
  // factor(lam) = c (lam - e)
  const roots = [{ e: [0, 0], c: [-1, 0] }]; // the overall -1 rides on the factor at 0
  for (const [ar, ai] of alphas) {
    const r2 = ar * ar + ai * ai;
    roots.push({ e: [ar, ai], c: [1, 0] });
    roots.push({ e: [ar / r2, ai / r2], c: [-ar, ai] }); // 1 - conj(a) lam = -conj(a) (lam - 1/conj(a))
  }
  const cycles = [];
  for (let j = 0; j < g; j++) cycles.push([0, 1 + 2 * j], [1 + 2 * j, 2 + 2 * j]);

  return cycles.map(([i1, i2]) => {
    const e1 = roots[i1].e, e2 = roots[i2].e;
    const mr = (e1[0] + e2[0]) / 2, mi = (e1[1] + e2[1]) / 2;
    const hr = (e2[0] - e1[0]) / 2, hi = (e2[1] - e1[1]) / 2;
    // constant from the two removed factors
    let kr = 1, ki = 0;
    for (const i of [i1, i2]) { const [cr, ci] = roots[i].c; const t = kr * cr - ki * ci; ki = kr * ci + ki * cr; kr = t; }
    const J = Array.from({ length: 2 * g + 1 }, () => [0, 0]); // oint lam^m dlam / nu
    let sPrev = null;
    for (let k = 0; k <= M; k++) {
      const phi = (Math.PI * k) / M;
      const wgt = (k === 0 || k === M ? 0.5 : 1) * (Math.PI / M);
      const c = Math.cos(phi);
      const lr = mr - hr * c, li = mi - hi * c;
      // rest(lam) = k * prod over the other factors
      let rr = kr, ri = ki;
      roots.forEach(({ e, c: cc }, i) => {
        if (i === i1 || i === i2) return;
        const fr = cc[0] * (lr - e[0]) - cc[1] * (li - e[1]), fi = cc[0] * (li - e[1]) + cc[1] * (lr - e[0]);
        const t = rr * fr - ri * fi; ri = rr * fi + ri * fr; rr = t;
      });
      // sqrt(rest), continued in phi
      const m = Math.sqrt(Math.hypot(rr, ri)), ph = Math.atan2(ri, rr) / 2;
      let sr = m * Math.cos(ph), si = m * Math.sin(ph);
      if (sPrev && (sr - sPrev[0]) ** 2 + (si - sPrev[1]) ** 2 > (sr + sPrev[0]) ** 2 + (si + sPrev[1]) ** 2) { sr = -sr; si = -si; }
      sPrev = [sr, si];
      // integrand of lam^m dlam/nu: -i lam^m / sqrt(rest); loop = 2 * segment
      const s2 = sr * sr + si * si;
      let br = -si / s2, bi = -sr / s2;                // -i / sqrt(rest)
      br *= 2 * wgt; bi *= 2 * wgt;
      for (let j = 0; j <= 2 * g; j++) {
        J[j][0] += br; J[j][1] += bi;
        const t = br * lr - bi * li; bi = br * li + bi * lr; br = t;
      }
    }
    // periods of lam^k dlam/(lam nu): k >= 1 is J_{k-1}; k = 0 via the exact form above
    const out = [[0, 0]];
    for (let k = 1; k <= g + 1; k++) out.push(J[k - 1]);
    let sr0 = 0, si0 = 0;
    for (let m = 0; m < 2 * g; m++) { sr0 += P[m][0] * J[m][0] - P[m][1] * J[m][1]; si0 += P[m][0] * J[m][1] + P[m][1] * J[m][0]; }
    if (g === 0) { sr0 = 0; si0 = 0; }
    out[0] = [(sr0 * a0r + si0 * a0i) / a02, (si0 * a0r - sr0 * a0i) / a02];
    return out;
  });
}

/** Coefficients of p_w (complex, powers 0..g+1) for w = [x, y], normalised by Re periods = 0. */
export function thetaPoly(alphas, w, periods = basicPeriods(alphas)) {
  const g = alphas.length;
  const a = aPoly(alphas);
  const kap = kappa0(alphas);
  const [wr, wi] = w;
  const a0r = a[0], a0i = a[1], aTr = a[4 * g], aTi = a[4 * g + 1];
  const p = Array.from({ length: g + 2 }, () => [0, 0]);
  // p0 = -w a0 / (2 kap), p_{g+1} = conj(w) a_{2g} / (2 kap)
  p[0] = [-(wr * a0r - wi * a0i) / (2 * kap), -(wr * a0i + wi * a0r) / (2 * kap)];
  p[g + 1] = [(wr * aTr + wi * aTi) / (2 * kap), (wr * aTi - wi * aTr) / (2 * kap)];
  if (g === 0) return p;
  // Re [ fixed + sum_k (x_k + i y_k) Pi_k ] = 0 on each of the 2g cycles
  const n = 2 * g;
  const Amat = [], rhs = [];
  for (const P of periods) {
    let fr = 0;
    for (const k of [0, g + 1]) fr += p[k][0] * P[k][0] - p[k][1] * P[k][1];
    const row = [];
    for (let k = 1; k <= g; k++) row.push(P[k][0], -P[k][1]);
    Amat.push(row); rhs.push(-fr);
  }
  const sol = solve(Amat, rhs, n);
  for (let k = 1; k <= g; k++) p[k] = [sol[2 * k - 2], sol[2 * k - 1]];
  return p;
}

/** Imaginary parts of the periods of Theta_w divided by 2 pi, one per cycle. */
export function periodVector(alphas, w, periods = basicPeriods(alphas)) {
  const p = thetaPoly(alphas, w, periods);
  return periods.map((P) => {
    let im = 0;
    for (let k = 0; k < p.length; k++) im += p[k][0] * P[k][1] + p[k][1] * P[k][0];
    return im / (2 * Math.PI);
  });
}

/**
 * The lattice Gamma of periods of zeta: w with periodVector(w) in Z^{2g}.
 * Searches integer vectors in the span of the period vectors of w = 1 and w = i.
 * Returns up to two generators (reduced), or [] if none found; for g <= 1 some directions are
 * continuous (flagged with continuous: [x, y] direction).
 */
export function periodLattice(alphas, { N = 8, tol = 1e-6 } = {}) {
  const periods = basicPeriods(alphas);
  const v1 = periodVector(alphas, [1, 0], periods), v2 = periodVector(alphas, [0, 1], periods);
  const n = v1.length;
  // least squares for x v1 + y v2 = m
  const g11 = dot(v1, v1), g12 = dot(v1, v2), g22 = dot(v2, v2);
  const det = g11 * g22 - g12 * g12;
  const found = [];
  if (Math.abs(det) < 1e-18 * (g11 * g22 + 1e-300)) return { v1, v2, generators: [], degenerate: true };
  const m = new Array(n).fill(-N);
  const step = () => { for (let i = 0; i < n; i++) { if (m[i] < N) { m[i]++; return true; } m[i] = -N; } return false; };
  do {
    if (m.every((x) => x === 0)) continue;
    const b1 = dot(v1, m), b2 = dot(v2, m);
    const x = (g22 * b1 - g12 * b2) / det, y = (g11 * b2 - g12 * b1) / det;
    let res = 0;
    for (let i = 0; i < n; i++) res = Math.max(res, Math.abs(x * v1[i] + y * v2[i] - m[i]));
    if (res < tol) found.push([x, y]);
  } while (step());
  // reduce to a basis: shortest vector, then shortest independent one
  found.sort((p, q) => Math.hypot(...p) - Math.hypot(...q));
  const gens = [];
  if (found.length) {
    gens.push(found[0]);
    const w1 = found[0];
    const second = found.find((w) => Math.abs(w1[0] * w[1] - w1[1] * w[0]) > 1e-6 * Math.hypot(...w1) * Math.hypot(...w));
    if (second) gens.push(gaussReduce(w1, second));
  }
  return { v1, v2, generators: gens };
}

function gaussReduce(w1, w2) {
  // make w2 as short as possible modulo w1
  const k = Math.round((w1[0] * w2[0] + w1[1] * w2[1]) / (w1[0] ** 2 + w1[1] ** 2));
  return [w2[0] - k * w1[0], w2[1] - k * w1[1]];
}
const dot = (a, b) => a.reduce((s, x, i) => s + x * b[i], 0);

function solve(A, b, n) {
  const M = A.map((r, i) => [...r, b[i]]);
  for (let c = 0; c < n; c++) {
    let piv = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r;
    [M[c], M[piv]] = [M[piv], M[c]];
    for (let r = 0; r < n; r++) {
      if (r === c) continue;
      const f = M[r][c] / M[c][c];
      for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k];
    }
  }
  return M.map((r, i) => r[n] / r[i]);
}
