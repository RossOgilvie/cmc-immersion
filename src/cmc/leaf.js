// The leaves of the Whitham foliation (Hasse–Ogilvie–Schmidt): the level sets of a -> BPer_a[B_a] in
// Gr(2, R^g), on which the period plane W of the differentials is fixed. Each leaf is 3-dimensional.
//
// A frame of the plane is a pair (w_1, w_2) of displacements, with b_k = p_{w_k}; along a leaf the frame
// moves so that its periods P_k = L(w_k) stay fixed (a parallel frame). Relative to a frame:
//
//   tau = w_2 / w_1 = b_2(0) / b_1(0),     W = 2i Res_0 q_{b_2} Theta_{b_1}   (the Willmore functional).
//
// The Whitham equations give, for a deformation with compatibility polynomial Q = A + r lam + conj(A) lam^2,
//   dW/dt = 4 Q'(0) = 4 r   and   dtau/dt = i a(0) A / (2 b_1(0)^2),
// so (tau, W) are coordinates on the leaf away from the set S^g where B_a has a common root: the tau-flow
// of whitham.js moves W at fixed tau (A = 0). W makes a poor slice across the tau-fibres, though (it varies
// little along a fibre, and is often critical on it), so TauFlow below moves tau at fixed
// psi = sum_j logit |alpha_j| instead, which in genus 2 is monotone on every fibre (see PLAN).
//
// Doubly periodic case: W is rational and Lambda = W ∩ Z^g is a lattice (always for g = 2, and the
// lattice of periods of zeta is then Gamma = L^{-1} Lambda). The cycles here (around {alpha_j, 1/conj alpha_j}
// and {0, alpha_j}) form an adapted basis of homology, so the paper's invariants apply: the triple
// (n_1, n_2, n_3), and tau of a lattice basis modulo a congruence subgroup, in F_1, F_3 or F_6.
// rationalPlanes() finds rational planes near W by LLL lattice reduction (used to recognise a rational
// plane); nearbyPeriods() lists integer vectors near W, and PlaneFlow moves the branch points so that W
// contains one or two of them ("locks them in"; two make W rational).

import { differentials } from './periods.js';
import { aPoly } from './cmc.js';

const M = 96; // minimum quadrature points per cycle, as in whitham.js
const imPeriod = (p, Pj) => p.reduce((acc, c, k) => acc + c[0] * Pj[k][1] + c[1] * Pj[k][0], 0) / (2 * Math.PI);
const cmul = ([ar, ai], [br, bi]) => [ar * br - ai * bi, ar * bi + ai * br];
const cdiv = ([ar, ai], [br, bi]) => { const d = br * br + bi * bi; return [(ar * br + ai * bi) / d, (ai * br - ar * bi) / d]; };
const dot = (a, b) => a.reduce((s, v, i) => s + v * b[i], 0);
const wrap = (t) => t - 2 * Math.PI * Math.round(t / (2 * Math.PI));

/**
 * The B-periods and the basic differentials at alphas: ell_j with L(w)_j = Re(ell_j w) (as ellVector of
 * whitham.js, signs aligned to ref), p1 = p_1 and pi = p_i (p_w = Re w p1 + Im w pi), and a(0).
 */
export function leafData(alphas, ref = null) {
  const D = differentials(alphas, M);
  const p1 = D.poly([1, 0]), pi = D.poly([0, 1]);
  const ell = alphas.map((_, j) => {
    let re = imPeriod(p1, D.B[j]), im = -imPeriod(pi, D.B[j]);
    if (ref && re * ref[j][0] + im * ref[j][1] < 0) { re = -re; im = -im; }
    return [re, im];
  });
  const a = aPoly(alphas);
  return { ell, p1, pi, a0: [a[0], a[1]] };
}

/** The period vector L(w) in R^g: the B-periods of Theta_w over 2 pi. */
export const periodsOf = (ell, [x, y]) => ell.map(([re, im]) => re * x - im * y);

/** The Willmore functional in the frame (w1, w2): W = 8i (b_0^1 b_1^2 - b_0^2 b_1^1) / a_0 (see whitham.js). */
export function willmoreFrame({ p1, pi, a0 }, w1, w2) {
  const b = (w) => [0, 1].map((k) => [w[0] * p1[k][0] + w[1] * pi[k][0], w[0] * p1[k][1] + w[1] * pi[k][1]]);
  const [b10, b11] = b(w1), [b20, b21] = b(w2);
  const u = cmul(b10, b21), v = cmul(b20, b11);
  return cmul([0, 8], cdiv([u[0] - v[0], u[1] - v[1]], a0))[0];
}

export const tauOf = ([w1, w2]) => cdiv(w2, w1);

/** The w with L(w) = m (least squares; exact when m lies in the period plane). */
export function frameVector(ell, m) {
  const v1 = ell.map((e) => e[0]), v2 = ell.map((e) => -e[1]);
  const g11 = dot(v1, v1), g12 = dot(v1, v2), g22 = dot(v2, v2), det = g11 * g22 - g12 * g12;
  const r1 = dot(v1, m), r2 = dot(v2, m);
  return [(g22 * r1 - g12 * r2) / det, (g11 * r2 - g12 * r1) / det];
}

// ------------------------------------------------------------------ small linear algebra

function solveSq(A, b) {
  const n = b.length, Mx = A.map((r, i) => [...r, b[i]]);
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(Mx[r][c]) > Math.abs(Mx[p][c])) p = r;
    [Mx[c], Mx[p]] = [Mx[p], Mx[c]];
    if (Math.abs(Mx[c][c]) < 1e-300) return null;
    for (let r = c + 1; r < n; r++) {
      const f = Mx[r][c] / Mx[c][c];
      for (let k = c; k <= n; k++) Mx[r][k] -= f * Mx[c][k];
    }
  }
  const x = new Array(n).fill(0);
  for (let r = n - 1; r >= 0; r--) {
    let s = Mx[r][n];
    for (let k = r + 1; k < n; k++) s -= Mx[r][k] * x[k];
    x[r] = s / Mx[r][r];
  }
  return x;
}

/** Orthonormal basis (Gram–Schmidt) of the span of the vectors; those dependent on earlier ones are dropped. */
function orthonormal(vs, eps = 1e-12) {
  const out = [];
  for (const v0 of vs) {
    let v = v0.slice();
    for (const u of out) { const d = dot(v, u); v = v.map((x, i) => x - d * u[i]); }
    const n = Math.hypot(...v);
    if (n > eps * Math.max(1, Math.hypot(...v0))) out.push(v.map((x) => x / n));
  }
  return out;
}

/** Orthonormal vectors completing the orthonormal set Q to a basis of R^g. */
function complement(Q, g) {
  const out = [...Q];
  for (let k = 0; k < g && out.length < g; k++) {
    let v = Array.from({ length: g }, (_, j) => (j === k ? 1 : 0));
    for (const u of out) { const d = dot(v, u); v = v.map((x, j) => x - d * u[j]); }
    const n = Math.hypot(...v);
    if (n > 0.3) out.push(v.map((x) => x / n));
  }
  return out.slice(Q.length);
}

/** Singular values [small, large] of the 2 x 2 matrix [[a, b], [c, d]]. */
function sv2(a, b, c, d) {
  const S1 = a * a + b * b + c * c + d * d, D = Math.abs(a * d - b * c);
  const big = Math.sqrt((S1 + Math.sqrt(Math.max(0, S1 * S1 - 4 * D * D))) / 2);
  return [big > 0 ? D / big : 0, big];
}

/**
 * Singular values of the rows `rows` of J restricted to the null space of the other rows (assumed
 * 2-dimensional): how well those two functions serve as coordinates on the solution set of the others.
 */
function restrictedSV(J, rows) {
  const n = J[0].length, others = J.filter((_, k) => !rows.includes(k));
  const scale = Math.max(...others.map((r) => Math.hypot(...r)));
  const B = orthonormal(others, 1e-6 * scale / Math.max(1, scale));
  const T = complement(B, n);
  if (T.length !== 2) return null;
  const [r1, r2] = rows.map((k) => T.map((t) => dot(J[k], t)));
  return sv2(r1[0], r1[1], r2[0], r2[1]);
}

// ------------------------------------------------------------------ the tau pad: moving tau at fixed psi

/** Which obstacle the branch points are near, if any: 'circle', 'zero' or 'collision' (else null). */
function endReason(al) {
  const rs = al.map(([re, im]) => Math.hypot(re, im));
  if (Math.max(...rs) > 0.95) return 'circle';
  if (Math.min(...rs) < 0.03) return 'zero';
  for (let i = 0; i < al.length; i++) for (let j = i + 1; j < al.length; j++) {
    if (Math.hypot(al[i][0] - al[j][0], al[i][1] - al[j][1]) < 0.02 * Math.max(rs[i], rs[j])) return 'collision';
  }
  return null;
}

// Unknowns per branch point (logit |alpha|, arg alpha), as in rootflow.js: alpha -> 0 and |alpha| -> 1
// are at infinity, and the gauge sum_j arg alpha_j = const is linear.
const R_MIN = 0.02, R_MAX = 0.98; // as for the tau-flow of whitham.js
const flat = (alphas) => alphas.flatMap(([re, im]) => { const r = Math.hypot(re, im); return [Math.log(r / (1 - r)), Math.atan2(im, re)]; });
const unflat = (x) => Array.from({ length: x.length / 2 }, (_, j) => {
  const r = 1 / (1 + Math.exp(-x[2 * j])), t = x[2 * j + 1];
  return [r * Math.cos(t), r * Math.sin(t)];
});

/**
 * Moves along the leaf through the given data so that tau reaches a target while
 * psi = sum_j logit |alpha_j| stays fixed. The tau-fibres (the tau-flow) are intervals that end where a
 * branch point reaches the unit circle or all reach 0, where psi -> +inf or -inf; in genus 2 psi was
 * monotone on every fibre tried, so psi = const meets each fibre once. (Holding W fixed instead, the
 * deformations with Q'(0) = 0, reaches the end of a fibre soon after leaving the start.) In genus >= 3
 * most fibres end at the circle at both ends, psi turns once on them, and the slice folds there.
 * w is the frame at alphas ([w1, w2], each [re, im]); it moves with the data (its periods stay fixed),
 * and tau is taken relative to it.
 *
 * Newton's method on the square system, in the unknowns (alpha, w_1, w_2), 2g + 4 reals:
 *   L(w_k) = P_k (2g),  gauge (1),  w_2 / w_1 = tau (2),  psi = psi_0 (1).
 * The unknowns for alpha are (logit |alpha_j|, arg alpha_j), so both psi and the gauge are linear.
 */
export class TauFlow {
  constructor(alphas, w) {
    this.g = alphas.length;
    if (this.g < 2) throw new Error('The τ pad needs genus ≥ 2');
    this.x = flat(alphas);
    this.w = w.map((v) => v.slice());
    const d = leafData(alphas);
    this.ref = d.ell; // cycle orientations, carried along
    this.P = this.w.map((wk) => periodsOf(d.ell, wk));
    this.psi0 = this.x.reduce((s, v, i) => s + (i % 2 ? 0 : v), 0);
    this.gauge0 = this.x.reduce((s, v, i) => s + (i % 2 ? v : 0), 0);
    this.tau = tauOf(this.w);
  }

  get alphas() { return unflat(this.x); }

  /** leafData at the unknowns x, or null outside the admissible region. */
  data(x) {
    const alphas = unflat(x), rs = alphas.map(([re, im]) => Math.hypot(re, im));
    if (!rs.every((r) => r > R_MIN && r < R_MAX)) return null;
    for (let i = 0; i < alphas.length; i++) for (let j = i + 1; j < alphas.length; j++) {
      if (Math.hypot(alphas[i][0] - alphas[j][0], alphas[i][1] - alphas[j][1]) < 2e-3 * Math.max(rs[i], rs[j])) return null;
    }
    return leafData(alphas, this.ref);
  }

  /** The residuals at (x, w) given d = data(x). */
  F(d, x, w, target) {
    const F = [];
    for (let k = 0; k < 2; k++) periodsOf(d.ell, w[k]).forEach((v, j) => F.push(v - this.P[k][j]));
    F.push(x.reduce((s, v, i) => s + (i % 2 ? v : 0), 0) - this.gauge0);
    const t = cdiv(w[1], w[0]);
    F.push(t[0] - target[0], t[1] - target[1]);
    F.push(x.reduce((s, v, i) => s + (i % 2 ? 0 : v), 0) - this.psi0);
    return F;
  }

  /** Jacobian of F at y = (x, w) (forward differences; the w-columns reuse d). */
  jacobian(y, d, F0, target) {
    const n = 2 * this.g, eps = 1e-7, J = F0.map(() => new Array(y.length).fill(0));
    for (let i = 0; i < y.length; i++) {
      const yp = y.slice(); yp[i] += eps;
      const x = yp.slice(0, n), w = [yp.slice(n, n + 2), yp.slice(n + 2)];
      const dp = i < n ? this.data(x) : d;
      if (!dp) return null;
      this.F(dp, x, w, target).forEach((v, k) => { J[k][i] = (v - F0[k]) / eps; });
    }
    return J;
  }

  /** Newton from the current point to tau = target. Returns { x, w, d } or null. */
  correct(target, tol = 1e-11) {
    const n = 2 * this.g;
    let y = [...this.x, ...this.w.flat()];
    const split = (v) => [v.slice(0, n), [v.slice(n, n + 2), v.slice(n + 2)]];
    let [x, w] = split(y), d = this.data(x);
    if (!d) return null;
    let F = this.F(d, x, w, target), res = Math.hypot(...F);
    for (let it = 0; it < 12 && res > tol; it++) {
      const J = this.jacobian(y, d, F, target);
      if (!J) return null;
      const dy = solveSq(J, F.map((v) => -v));
      if (!dy || !dy.every(Number.isFinite)) return null;
      let lam = 1, next = null;
      for (let k = 0; k < 8; k++, lam /= 2) { // damped: halve until the residual decreases
        const yt = y.map((v, i) => v + lam * dy[i]);
        const [xt, wt] = split(yt), dt = this.data(xt);
        if (!dt) continue;
        const Ft = this.F(dt, xt, wt, target), rt = Math.hypot(...Ft);
        if (rt < res) { next = { y: yt, x: xt, w: wt, d: dt, F: Ft, res: rt }; break; }
      }
      if (!next) break;
      ({ y, x, w, d, F, res } = next);
    }
    return res < 1e-8 ? { x, w, d } : null;
  }

  /**
   * Moves towards the target tau along the straight segment, in steps of at most `step` in the
   * hyperbolic metric of the upper half plane, for at most `budget` ms (as RootFlow.moveToward).
   * Returns { alphas, w, tau, reached, blocked, why }.
   */
  moveToward(target, { step = 0.04, budget = 12 } = {}) {
    const t0 = performance.now();
    if (!this.lastTarget || this.lastTarget[0] !== target[0] || this.lastTarget[1] !== target[1]) {
      this.lastTarget = target.slice();
      this.h = step;
    }
    let h = this.h, blocked = false;
    while (performance.now() - t0 < budget) {
      const dt = [target[0] - this.tau[0], target[1] - this.tau[1]], dist = Math.hypot(...dt);
      if (dist < 1e-10) break;
      const cap = step * this.tau[1]; // hyperbolic step
      h = Math.min(h, cap);
      const f = Math.min(1, h / dist);
      const want = [this.tau[0] + f * dt[0], this.tau[1] + f * dt[1]];
      const r = this.correct(want);
      if (r && Math.hypot(...r.x.map((v, i) => v - this.x[i])) < 0.5) {
        this.x = r.x;
        this.w = r.w;
        this.ref = r.d.ell;
        this.tau = tauOf(this.w);
        h = Math.min(cap, h * 2);
      } else if (h > cap / 1024) {
        h /= 4;
      } else {
        blocked = true;
        break;
      }
    }
    this.h = h;
    if (blocked) this.lastTarget = null;
    const reached = Math.hypot(target[0] - this.tau[0], target[1] - this.tau[1]) < 1e-9;
    return { alphas: this.alphas, w: this.w.map((v) => v.slice()), tau: this.tau.slice(), reached, blocked, why: blocked ? this.why() : null };
  }

  /**
   * Why the pad stopped: 'circle' | 'zero' | 'collision' (a branch point near the unit circle or 0, or
   * two meeting), 'fold' (tau restricted to the level set of psi in the leaf nearly drops rank: there
   * psi is critical along a tau-fibre, which happens in genus >= 3), or 'solver'.
   */
  why() {
    const e = endReason(this.alphas);
    if (e) return e;
    const s = this.tauSingularValues();
    return s && s[0] < 0.02 * s[1] ? 'fold' : 'solver';
  }

  /** Singular values of dtau on the tangent plane of {psi = psi_0} in the leaf, at the current point. */
  tauSingularValues() {
    const n = 2 * this.g, y = [...this.x, ...this.w.flat()], d = this.data(this.x);
    if (!d) return null;
    const F0 = this.F(d, this.x, this.w, this.tau), J = this.jacobian(y, d, F0, this.tau);
    return J ? restrictedSV(J, [n + 1, n + 2]) : null;
  }
}

// ------------------------------------------------------------------ rational planes near W (LLL)

/**
 * LLL reduction (delta = 0.99) of the lattice spanned by embed(m) for the integer vectors m in ms. The
 * reduction is done on the integer vectors, so the result is again a basis of the lattice they span.
 */
export function lll(ms, embed, delta = 0.99) {
  const k = ms.length, m = ms.map((r) => r.slice());
  let b = m.map(embed), Bs, mu, B;
  const gso = () => {
    Bs = []; mu = []; B = [];
    for (let i = 0; i < k; i++) {
      let v = b[i].slice();
      mu.push([]);
      for (let j = 0; j < i; j++) {
        mu[i][j] = B[j] > 0 ? dot(b[i], Bs[j]) / B[j] : 0;
        v = v.map((x, t) => x - mu[i][j] * Bs[j][t]);
      }
      Bs.push(v);
      B.push(dot(v, v));
    }
  };
  gso();
  let i = 1;
  for (let guard = 0; i < k && guard < 100000; guard++) {
    for (let j = i - 1; j >= 0; j--) {
      const q = Math.round(mu[i][j]);
      if (q) { m[i] = m[i].map((v, t) => v - q * m[j][t]); b[i] = embed(m[i]); gso(); }
    }
    if (B[i] >= (delta - mu[i][i - 1] ** 2) * B[i - 1]) i++;
    else {
      [m[i], m[i - 1]] = [m[i - 1], m[i]];
      [b[i], b[i - 1]] = [b[i - 1], b[i]];
      gso();
      i = Math.max(1, i - 1);
    }
  }
  return m;
}

/** sin of the largest principal angle between the planes with orthonormal bases Q and R. */
function planeAngle(Q, R) {
  const [s] = sv2(dot(Q[0], R[0]), dot(Q[0], R[1]), dot(Q[1], R[0]), dot(Q[1], R[1]));
  return Math.sqrt(Math.max(0, 1 - Math.min(1, s) ** 2));
}

/** The triple (n_1 >= n_2 >= n_3) of the lattice with basis m1, m2: the tally of nonzero (m1_j, m2_j) mod 2. */
export function invariantTriple(m1, m2) {
  const n = [0, 0, 0];
  m1.forEach((v, j) => {
    const p = Math.abs(v) % 2, q = Math.abs(m2[j]) % 2;
    if (p || q) n[p && q ? 2 : p ? 0 : 1]++;
  });
  return n.sort((a, b) => b - a);
}

/** The first two vectors of the LLL-reduced basis of Z^g mapped by m -> (m, C N m), for C = 10^(k/2), k = 1..scales. */
function lllSweep(ell, scales = 18) {
  const g = ell.length;
  const Q = orthonormal([ell.map((v) => v[0]), ell.map((v) => -v[1])]);
  if (Q.length < 2) return [];
  const N = complement(Q, g);
  const e = (i) => Array.from({ length: g }, (_, j) => (j === i ? 1 : 0));
  const out = [];
  for (let k = 1; k <= scales; k++) {
    const C = 10 ** (k / 2);
    const red = lll(Array.from({ length: g }, (_, i) => e(i)), (m) => [...m, ...N.map((nv) => C * dot(nv, m))]);
    out.push([red[0], red[1]]);
  }
  return out;
}

/**
 * The integer vectors that LLL finds near the period plane (both vectors of every scale of lllSweep, and
 * their negatives, without repeats), as nearbyPeriods: [{ m, dist, w }]. They are longer than the box of
 * nearbyPeriods reaches, and much closer to the plane.
 */
export function lllPeriods(ell, { maxNorm = 1e5 } = {}) {
  if (ell.length < 3) return [];
  const v1 = ell.map((e) => e[0]), v2 = ell.map((e) => -e[1]), seen = new Set(), out = [];
  for (const pair of lllSweep(ell)) for (const m0 of pair) for (const m of [m0, m0.map((v) => -v)]) {
    const key = m.join(',');
    if (seen.has(key) || Math.hypot(...m) > maxNorm) continue;
    seen.add(key);
    const P = project(v1, v2, m);
    out.push({ m, dist: Math.hypot(...m.map((v, j) => v - P[j])), w: frameVector(ell, m) });
  }
  return out;
}

/**
 * Rational planes near the period plane W = span(L(1), L(i)) ⊂ R^g, from the B-periods ell.
 *
 * With N the (g-2) x g matrix of orthonormal normals of W and a weight C, an integer vector m is mapped
 * to (m, C N m) in R^{2g-2}: its length is small when m is short and nearly in W. LLL reduction of the
 * image of Z^g finds such vectors, and its first two vectors m_1, m_2 are part of a basis of Z^g, so they
 * are a basis of Lambda* = span(m_1, m_2) ∩ Z^g. Large C asks for planes closer to W, at a larger height.
 * Sweeping C gives a sequence of approximations, like the convergents of a continued fraction.
 *
 * Returns [{ m: [m1, m2], height, dist, n }] sorted by height, keeping only those closer to W than every
 * plane of smaller height. height = area of the fundamental domain of Lambda*; dist = sin of the largest
 * angle between W and the plane; n = the invariant triple. For g = 2, W = R^2 and Lambda = Z^2.
 */
export function rationalPlanes(ell, { maxHeight = 1e5, scales = 18 } = {}) {
  const g = ell.length;
  const e = (i) => Array.from({ length: g }, (_, j) => (j === i ? 1 : 0));
  if (g === 2) return [{ m: [e(0), e(1)], height: 1, dist: 0, n: [1, 1, 0] }];
  if (g < 2) return [];
  const Q = orthonormal([ell.map((v) => v[0]), ell.map((v) => -v[1])]);
  if (Q.length < 2) return [];
  const found = [];
  for (const [m1, m2] of lllSweep(ell, scales)) {
    const R = orthonormal([m1, m2]);
    if (R.length < 2) continue;
    const height = Math.sqrt(Math.max(0, dot(m1, m1) * dot(m2, m2) - dot(m1, m2) ** 2));
    if (height > maxHeight) continue;
    // the same plane twice: keep the first (the lattices are equal, both being saturated)
    if (found.some((f) => planeAngle(f.R, R) < 1e-9)) continue;
    found.push({ m: [m1, m2], R, height, dist: planeAngle(Q, R), n: invariantTriple(m1, m2) });
  }
  found.sort((a, b) => a.height - b.height || a.dist - b.dist);
  const out = [];
  for (const f of found) if (!out.length || f.dist < out[out.length - 1].dist) out.push(f);
  return out.map(({ m, height, dist, n }) => ({ m, height, dist, n }));
}

// ------------------------------------------------------------------ lattice bases and the domains F_6, F_3, F_1

const oddCount = (m) => m.reduce((s, v) => s + (Math.abs(v) % 2), 0);

/**
 * Reduces the lattice basis (w_1, w_2) (with period vectors m_1, m_2) so that tau = w_2 / w_1 lies in
 * F_6 (by Gamma(2), which preserves the classes mod 2), F_3 (also tau -> -1/tau, swapping w_1 and w_2)
 * or F_1 (all of SL(2, Z)). The determinant is kept at 1, so W is unchanged.
 */
export function reduceFrame(basis, kind) {
  let [b1, b2] = basis.map(({ w, m }) => ({ w: w.slice(), m: m.slice() }));
  const comb = (a, p, b, q) => ({ w: [p * a.w[0] + q * b.w[0], p * a.w[1] + q * b.w[1]], m: a.m.map((v, j) => p * v + q * b.m[j]) });
  for (let it = 0; it < 200; it++) {
    const [x, y] = cdiv(b2.w, b1.w);
    if (kind === 'F1') {
      const k = Math.round(x);
      if (k) { b2 = comb(b2, 1, b1, -k); continue; }
      if (x * x + y * y < 1 - 1e-12) { [b1, b2] = [b2, comb(b1, -1, b1, 0)]; continue; } // tau -> -1/tau
      break;
    }
    const k = 2 * Math.round(x / 2);
    if (k && Math.abs(x) > 1 + 1e-12) { b2 = comb(b2, 1, b1, -k); continue; }
    if (Math.hypot(x + 0.5, y) < 0.5 - 1e-12) { b1 = comb(b1, 1, b2, 2); continue; }  // tau -> tau / (2 tau + 1)
    if (Math.hypot(x - 0.5, y) < 0.5 - 1e-12) { b1 = comb(b1, 1, b2, -2); continue; } // tau -> tau / (1 - 2 tau)
    if (kind === 'F3' && x * x + y * y < 1 - 1e-12) { [b1, b2] = [b2, comb(b1, -1, b1, 0)]; continue; }
    break;
  }
  return [b1, b2];
}

/**
 * A basis of the lattice Lambda = span(m1, m2) ∩ Z^g adapted to its invariants (the proof of the
 * theorem on CA(g)-orbits of complex structures): the three nonzero classes of Lambda / 2 Lambda are told
 * apart by how many of their B-periods are odd (n_2 + n_3, n_1 + n_3, n_1 + n_2 for v_1, v_2, v_3).
 * With all three counts distinct, the basis is (v_1, v_2) and tau is defined modulo Gamma(2) (domain F_6);
 * with two equal, the basis is the tied pair, which may be swapped (F_3); with all equal, F_1.
 * Returns { w: [w1, w2], m: [m1, m2] (their period vectors), n, kind, tau } with tau in the domain.
 */
export function latticeFrame(ell, m1, m2) {
  const cls = [m1, m2, m1.map((v, j) => v + m2[j])].map((m) => ({ m, s: oddCount(m) })).sort((a, b) => a.s - b.s);
  let kind, pair;
  if (cls[0].s === cls[2].s) { kind = 'F1'; pair = [cls[0], cls[1]]; }
  else if (cls[0].s === cls[1].s) { kind = 'F3'; pair = [cls[0], cls[1]]; }
  else if (cls[1].s === cls[2].s) { kind = 'F3'; pair = [cls[1], cls[2]]; }
  else { kind = 'F6'; pair = [cls[0], cls[1]]; }
  let basis = pair.map(({ m }) => ({ w: frameVector(ell, m), m }));
  if (cdiv(basis[1].w, basis[0].w)[1] < 0) basis[1] = { w: basis[1].w.map((v) => -v), m: basis[1].m.map((v) => -v) };
  basis = reduceFrame(basis, kind);
  const w = basis.map((b) => b.w);
  return { w, m: basis.map((b) => b.m), n: invariantTriple(m1, m2), kind, tau: tauOf(w) };
}

/**
 * The region of F_6 that tau lies in, as the order of the classes [w_1] = 1, [w_2] = 2, [w_1 + w_2] = 3
 * of Lambda / 2 Lambda by the length of their shortest vectors ('123' on F_1 itself). The six regions
 * are the six copies of F_1 in F_6.
 */
export function classOrder([x, y]) {
  const best = [Infinity, Infinity, Infinity];
  for (let p = -9; p <= 9; p++) for (let q = -9; q <= 9; q++) {
    const pp = Math.abs(p) % 2, qq = Math.abs(q) % 2;
    if (!pp && !qq) continue;
    const c = pp && !qq ? 0 : !pp && qq ? 1 : 2;
    best[c] = Math.min(best[c], Math.hypot(p + q * x, q * y));
  }
  return [0, 1, 2].sort((a, b) => best[a] - best[b]).map((c) => c + 1).join('');
}

// ------------------------------------------------------------------ locking integer period vectors into the plane

/** Orthogonal projection of m onto span(v1, v2). */
function project(v1, v2, m) {
  const g11 = dot(v1, v1), g12 = dot(v1, v2), g22 = dot(v2, v2), det = g11 * g22 - g12 * g12;
  const r1 = dot(v1, m), r2 = dot(v2, m);
  const x = (g22 * r1 - g12 * r2) / det, y = (g11 * r2 - g12 * r1) / det;
  return m.map((_, j) => x * v1[j] + y * v2[j]);
}

/**
 * Integer vectors near the period plane W = span(L(1), L(i)), for the picture of the plane: every m in
 * the box |m_j| <= N (N from the genus, so the box has at most ~10^5 points), with its distance to W and
 * the displacement w with L(w) the projection of m. Returns the `count` closest (both of ±m), excluding 0;
 * among the periods (distance 0) the shortest.
 */
export function nearbyPeriods(ell, { count = 80 } = {}) {
  const g = ell.length, N = { 1: 12, 2: 8, 3: 8, 4: 5, 5: 4, 6: 3 }[g] || 3;
  const v1 = ell.map((e) => e[0]), v2 = ell.map((e) => -e[1]);
  const g11 = dot(v1, v1), g12 = dot(v1, v2), g22 = dot(v2, v2), det = g11 * g22 - g12 * g12;
  // dist^2 = |m|^2 - r^T G^{-1} r, r = (v1 . m, v2 . m), accumulated along the odometer without allocation.
  // This cancels badly near the plane (an error ~ 1e-16 |m|^2), so it only screens; the kept distances
  // are recomputed by projection.
  const tol = 1e-9, noise = 1e-12;
  let keep = [], cut = Infinity;
  const m = new Array(g).fill(-N);
  const next = () => { for (let i = 0; i < g; i++) { if (m[i] < N) { m[i]++; return true; } m[i] = -N; } return false; };
  do {
    let mm = 0, r1 = 0, r2 = 0, zero = true;
    for (let j = 0; j < g; j++) { const v = m[j]; if (v) zero = false; mm += v * v; r1 += v1[j] * v; r2 += v2[j] * v; }
    if (zero) continue;
    const d2 = Math.max(0, mm - (g22 * r1 * r1 - 2 * g12 * r1 * r2 + g11 * r2 * r2) / det);
    // possible periods are always kept (ordered by length below); others only below the running cut
    if (d2 > noise * mm && d2 >= cut) continue;
    keep.push({ m: m.slice(), d2, mm });
    if (keep.length > 4 * count) {
      keep.sort((a, b) => a.d2 - b.d2);
      keep = keep.slice(0, count);
      cut = keep[count - 1].d2;
    }
  } while (next());
  const out = keep.map(({ m: mv, mm }) => {
    const P = project(v1, v2, mv), dist = Math.hypot(...mv.map((v, j) => v - P[j]));
    return { m: mv, dist, w: frameVector(ell, mv), period: dist <= tol * Math.sqrt(mm) };
  });
  // periods (distance 0 up to rounding) first, nearest the origin first; then by distance
  out.sort((a, b) => (a.period && b.period ? Math.hypot(...a.w) - Math.hypot(...b.w) : a.period ? -1 : b.period ? 1 : a.dist - b.dist));
  return out.slice(0, count).map(({ m: mv, dist, w }) => ({ m: mv, dist, w }));
}

/**
 * A basis of the saturated lattice span_R(m1, m2) ∩ Z^g (the lattice Lambda of the rational plane
 * spanned by m1, m2): LLL with a large weight on the normals of the plane, whose first two vectors are
 * then the lattice's shortest (see rationalPlanes).
 */
export function saturate(m1, m2) {
  const g = m1.length, N = complement(orthonormal([m1, m2]), g);
  const C = 1e6 * Math.max(1, Math.hypot(...m1), Math.hypot(...m2));
  const e = (i) => Array.from({ length: g }, (_, j) => (j === i ? 1 : 0));
  const red = lll(Array.from({ length: g }, (_, i) => e(i)), (m) => [...m, ...N.map((nv) => C * dot(nv, m))]);
  return [red[0], red[1]];
}

/**
 * A proxy for the chance that locking a point in succeeds: the size |Δalpha| (in the λ-plane) of one
 * minimum-norm Gauss–Newton step of PlaneFlow towards it from the current data, i.e. the move predicted
 * by linearising. The equations for all points share one linearisation of the period plane (2g
 * evaluations of the periods), and each point's step is then a small least-squares problem.
 * ms: the candidate vectors; locked: vectors already in the plane, which stay there; ell: the B-periods at
 * alphas in the orientation of both. Returns one value per candidate (Infinity where the step is singular,
 * e.g. a candidate dependent on the locked ones).
 */
export function predictMoves(alphas, ell, ms, locked = []) {
  const g = alphas.length, y0 = alphas.flat(), n = 2 * g, eps = 1e-7;
  const vecs = (e) => [e.map((x) => x[0]), e.map((x) => -x[1])];
  const N0 = complement(orthonormal(vecs(ell)), g);
  const planes = [vecs(ell)];
  for (let i = 0; i < n; i++) {
    const yp = y0.slice(); yp[i] += eps;
    planes.push(vecs(leafData(Array.from({ length: g }, (_, j) => [yp[2 * j], yp[2 * j + 1]]), ell).ell));
  }
  // rows of the residual N0 (m - proj_W m) / |m| and of their derivatives in alpha
  const rows = (m) => {
    const nm = Math.hypot(...m);
    const r = planes.map(([v1, v2]) => { const P = project(v1, v2, m); return N0.map((nv) => m.reduce((s, v, j) => s + nv[j] * (v - P[j]), 0) / nm); });
    return r[0].map((F, k) => ({ F, J: Array.from({ length: n }, (_, i) => (r[i + 1][k] - F) / eps) }));
  };
  const gauge = { F: 0, J: alphas.flatMap(([x, y]) => { const r2 = x * x + y * y; return [-y / r2, x / r2]; }) };
  const lockRows = locked.flatMap(rows);
  return ms.map((m) => {
    const R = [...lockRows, ...rows(m), gauge], k = R.length;
    const JJt = R.map((a) => R.map((b) => dot(a.J, b.J)));
    const u = solveSq(JJt, R.map((r) => -r.F));
    if (!u || !u.every(Number.isFinite)) return Infinity;
    const step = Array.from({ length: n }, (_, i) => R.reduce((s, r, a) => s + r.J[i] * u[a], 0));
    const size = Math.hypot(...step);
    return k && Number.isFinite(size) ? size : Infinity;
  });
}

/**
 * Moves the branch points so that the period plane contains the integer vectors ms (one or two; with two
 * independent ones the plane becomes span(m1, m2), a rational plane). For each m the residual is
 * N0 (m - proj_W m) / |m|, with N0 the normals of the plane at the start (g - 2 rows each; near the start
 * they detect the part of m off the plane), plus the gauge. By continuation: the residual is taken from
 * F(alpha_0) to 0 as (1 - t) F(alpha_0), t: 0 -> 1, each step by minimum-norm Gauss–Newton in the branch
 * points alpha (so "nearest" means nearest in the λ-plane). Rows for an m already in the plane start at 0
 * and stay there. ref: the B-periods ell at alphas whose orientation the ms refer to.
 */
export class PlaneFlow {
  constructor(alphas, ms, ref) {
    this.g = alphas.length;
    this.y = alphas.flat();
    this.arg0 = alphas.map(([re, im]) => Math.atan2(im, re));
    this.ms = ms.map((m) => m.slice());
    this.N0 = complement(orthonormal([ref.map((e) => e[0]), ref.map((e) => -e[1])]), this.g);
    this.ref = ref;
    const r = this.residual(this.y);
    if (!r) throw new Error('Out of range');
    this.F0 = r.F;
    this.t = 0;
    this.h = 0.25;
  }

  get alphas() { return Array.from({ length: this.g }, (_, j) => [this.y[2 * j], this.y[2 * j + 1]]); }

  residual(y, ref = this.ref) {
    const alphas = Array.from({ length: this.g }, (_, j) => [y[2 * j], y[2 * j + 1]]);
    const rs = alphas.map(([re, im]) => Math.hypot(re, im));
    if (!rs.every((r) => r > R_MIN && r < R_MAX)) return null;
    // two branch points meeting: the curve degenerates, and the quadrature of the periods refines without
    // bound (a single evaluation can take seconds), so stop well before
    for (let i = 0; i < alphas.length; i++) for (let j = i + 1; j < alphas.length; j++) {
      if (Math.hypot(alphas[i][0] - alphas[j][0], alphas[i][1] - alphas[j][1]) < 0.01 * Math.max(rs[i], rs[j])) return null;
    }
    const d = leafData(alphas, ref);
    const v1 = d.ell.map((e) => e[0]), v2 = d.ell.map((e) => -e[1]);
    const F = [];
    for (const m of this.ms) {
      const P = project(v1, v2, m), n = Math.hypot(...m);
      for (const nv of this.N0) F.push(m.reduce((s, v, j) => s + nv[j] * (v - P[j]), 0) / n);
    }
    F.push(alphas.reduce((s, [re, im], j) => s + wrap(Math.atan2(im, re) - this.arg0[j]), 0));
    return { F, d };
  }

  /** Minimum-norm Gauss–Newton from y to the residual (1 - t) F0. Returns { y, d } or null. */
  solve(y0, t, tol = 1e-12, deadline = Infinity) {
    const target = this.F0.map((v, k) => (k < this.F0.length - 1 ? (1 - t) * v : 0));
    let y = y0.slice(), r = this.residual(y);
    if (!r) return null;
    let G = r.F.map((v, k) => v - target[k]), res = Math.hypot(...G);
    for (let it = 0; it < 15 && res > tol; it++) {
      const eps = 1e-7, cols = [];
      for (let i = 0; i < y.length; i++) {
        if (performance.now() > deadline) return null;
        const yp = y.slice(); yp[i] += eps;
        const rp = this.residual(yp, r.d.ell);
        if (!rp) return null;
        cols.push(rp.F.map((v, k) => (v - r.F[k]) / eps));
      }
      // dy = -J^T (J J^T)^{-1} G, J[k][i] = cols[i][k]
      const m = G.length;
      const JJt = Array.from({ length: m }, (_, a) => Array.from({ length: m }, (_, b) => cols.reduce((s, c) => s + c[a] * c[b], 0)));
      const u = solveSq(JJt, G.map((v) => -v));
      if (!u || !u.every(Number.isFinite)) return null;
      const dy = cols.map((c) => dot(c, u));
      let lam = 1, next = null;
      for (let k = 0; k < 8; k++, lam /= 2) {
        const yt = y.map((v, i) => v + lam * dy[i]);
        const rt = this.residual(yt, r.d.ell);
        if (!rt) continue;
        const Gt = rt.F.map((v, kk) => v - target[kk]), nt = Math.hypot(...Gt);
        if (nt < res) { next = { y: yt, r: rt, G: Gt, res: nt }; break; }
      }
      if (!next) break;
      ({ y, r, G, res } = next);
    }
    return res < 1e-9 ? { y, d: r.d } : null;
  }

  /**
   * Advances t towards 1 for about `budget` ms. A single step that takes longer than `stepLimit` ms (the
   * periods have become expensive: branch points close together or close to a cycle) ends the flow.
   * Returns { alphas, t, done, blocked, why } with why as for TauFlow ('circle' | 'zero' | 'collision' |
   * 'cost' | 'solver').
   */
  advance(budget = 12, stepLimit = 250) {
    const t0 = performance.now();
    let blocked = false, why = null;
    while (this.t < 1 && performance.now() - t0 < budget) {
      const tn = Math.min(1, this.t + this.h), ts = performance.now();
      const r = this.solve(this.y, tn, 1e-12, ts + stepLimit);
      if (performance.now() - ts > stepLimit) { blocked = true; why = 'cost'; break; }
      if (r && Math.hypot(...r.y.map((v, i) => v - this.y[i])) < 0.3) {
        this.y = r.y;
        this.ref = r.d.ell;
        this.t = tn;
        this.h = Math.min(0.5, this.h * 1.5);
      } else if (this.h > 1e-4) {
        this.h /= 4;
      } else {
        blocked = true;
        break;
      }
    }
    if (blocked && !why) why = endReason(this.alphas) || 'solver';
    return { alphas: this.alphas, t: this.t, done: this.t >= 1, blocked, why };
  }

  /** The B-periods at the current point, in the orientation of the ms. */
  data() { return leafData(this.alphas, this.ref); }
}
