// The root-preserving Whitham flow on S^g (Carberry–Kilian–Klein–Schmidt): deformations of the
// spectral curve that keep the plane of periods of the differentials Theta_b, b in B_a, and keep a
// common root of all of B_a at the Sym point lam0 on the unit circle. On S^g the flow is 2-dimensional,
// and its coordinates are the Sym integrals
//
//   phi_l = -i q_l(y),   q_l(y) = 1/2 int_{sigma(y)}^{y} Theta(b_l)  along the real locus over |lam| = 1,
//
// for the basis (b_1, b_2) of B_a whose B-periods are 2 pi i delta_kl (so the periods, and with them
// the flow's normalisation, stay fixed). For g = 2 the map to (phi_1, phi_2) is a global diffeomorphism
// onto a triangle (CKKS, Theorem 1.3), and the tori with Sym point lam0 are exactly phi in (pi Q)^2.
//
// Here the flow is not integrated as an ODE: since phi are coordinates, moving to a target phi is a
// root-finding problem in the branch points alpha (pinned common root at lam0 + phi = target), solved by
// Gauss–Newton, following a straight path in phi in small steps so the solution stays on one branch.
// Implemented for even g (the real locus over |lam| = 1 is then one circle, and going once around it
// runs from sigma(y) to y); the definition for odd g is open.
//
// Genus >= 4: the period plane (the span of the B-period vectors of B_a in R^g) is a real constraint,
// 2(g - 2) equations, and there is no canonical basis: b_1, b_2 are the differentials whose B-periods are
// fixed vectors P_1, P_2 spanning the plane (by default orthonormal; phi changes by a linear map with P).
// Unknowns 2g, equations: common root (2 effective) + plane (2g - 4) + phi (2), so the flow is still
// 2-dimensional. Its image in the phi-plane is not known, so the pad marks where the flow stops.

import { aPoly } from './cmc.js';
import { basicPeriods, thetaPolyB } from './periods.js';

const cmul = ([a, b], [c, d]) => [a * c - b * d, a * d + b * c];
function evalC(p, z) { // p: array of [re, im] by increasing power
  let r = [0, 0];
  for (let k = p.length - 1; k >= 0; k--) { r = cmul(r, z); r[0] += p[k][0]; r[1] += p[k][1]; }
  return r;
}
function evalF(a, z) { // a: interleaved Float64Array
  let r = [0, 0];
  for (let k = a.length / 2 - 1; k >= 0; k--) { r = cmul(r, z); r[0] += a[2 * k]; r[1] += a[2 * k + 1]; }
  return r;
}

const M_PERIODS = 128; // minimum quadrature points per cycle (basicPeriods refines near other branch points)

/** Gauss–Legendre nodes and weights on [-1, 1]. */
function gaussLegendre(n) {
  const x = [], w = [];
  for (let i = 0; i < n; i++) {
    let z = Math.cos((Math.PI * (i + 0.75)) / (n + 0.5)), dp = 0;
    for (let it = 0; it < 100; it++) {
      let p0 = 1, p1 = z;
      for (let k = 2; k <= n; k++) { const p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
      dp = (n * (z * p1 - p0)) / (z * z - 1);
      const dz = p1 / dp;
      z -= dz;
      if (Math.abs(dz) < 1e-16) break;
    }
    x.push(z); w.push(2 / ((1 - z * z) * dp * dp));
  }
  const order = x.map((_, i) => i).sort((i, j) => x[i] - x[j]);
  return { x: order.map((i) => x[i]), w: order.map((i) => w[i]) };
}
const GL = gaussLegendre(12);

/** Panels covering [theta0, theta0 + 2 pi], graded towards arg alpha_j at the scale |1 - |alpha_j||. */
function circlePanels(alphas, theta0) {
  const T = 2 * Math.PI, pts = [0, T];
  for (const [re, im] of alphas) {
    const d = Math.max(1e-7, Math.abs(1 - Math.hypot(re, im)));
    let c = Math.atan2(im, re) - theta0;
    c -= T * Math.floor(c / T);
    pts.push(c);
    for (let o = d; o < Math.PI; o *= 2) pts.push(c - o, c + o);
  }
  const u = [...new Set(pts.filter((t) => t >= 0 && t <= T).map((t) => +t.toFixed(15)))].sort((p, q) => p - q);
  const panels = [];
  for (let i = 1; i < u.length; i++) {
    const n = Math.max(1, Math.ceil((u[i] - u[i - 1]) / 0.3));
    for (let k = 0; k < n; k++) {
      const t0 = u[i - 1] + ((u[i] - u[i - 1]) * k) / n, t1 = u[i - 1] + ((u[i] - u[i - 1]) * (k + 1)) / n;
      if (t1 - t0 > 1e-14) panels.push([theta0 + t0, theta0 + t1]);
    }
  }
  return panels;
}

/** An orthonormal basis of the span of the columns of Bm (g x 2, rows [v1_j, v2_j]). */
export function planeBasis(Bm) {
  const c1 = Bm.map((r) => r[0]), c2 = Bm.map((r) => r[1]);
  const n1 = Math.hypot(...c1), P1 = c1.map((v) => v / n1);
  const d = c2.reduce((s, v, j) => s + v * P1[j], 0);
  const q = c2.map((v, j) => v - d * P1[j]), n2 = Math.hypot(...q);
  return [P1, q.map((v) => v / n2)];
}

/** Orthonormal vectors completing the plane spanned by the orthonormal pair P to a basis of R^g. */
function normals(P) {
  const g = P[0].length, out = [...P];
  for (let k = 0; k < g && out.length < g; k++) {
    let v = Array.from({ length: g }, (_, j) => (j === k ? 1 : 0));
    for (const u of out) { const d = v.reduce((s, x, j) => s + x * u[j], 0); v = v.map((x, j) => x - d * u[j]); }
    const n = Math.hypot(...v);
    if (n > 0.3) out.push(v.map((x) => x / n));
  }
  return out.slice(2);
}

/**
 * The spectral data at alphas with Sym point lam0 = e^{i theta0}: the normalised differentials p_1, p_i,
 * their B-periods Bm (g x 2, divided by 2 pi), the basis b_l with B-periods P_l, and the Sym integrals
 * phi_l. For g = 2, P = delta_kl and phi = |phi_l| (the sign depends only on orientation choices, and CKKS
 * choose it positive). For g >= 4, P defaults to an orthonormal basis of the period plane and phi keeps
 * its sign (the caller fixes the orientation of P once). Also the values of p_1, p_i at lam0 (zero on S^g).
 * The orientation of each B-cycle is arbitrary and can flip as alpha moves; with ref (an earlier Bm),
 * each row of Bm is aligned with it, so that Bm (and with it the plane and the basis) moves continuously.
 */
export function symData(alphas, theta0, P = null, ref = null) {
  const g = alphas.length;
  // only the B-cycles (loops around {0, alpha_j}) are needed: the A-periods vanish by the reality condition
  const Per = basicPeriods(alphas, M_PERIODS, 'B');
  const p1 = thetaPolyB(alphas, [1, 0], Per), pi = thetaPolyB(alphas, [0, 1], Per);
  const imPeriod = (p, Pj) => p.reduce((acc, c, k) => acc + c[0] * Pj[k][1] + c[1] * Pj[k][0], 0) / (2 * Math.PI);
  const Bm = Per.map((Pj) => [imPeriod(p1, Pj), imPeriod(pi, Pj)]);
  if (ref) Bm.forEach((r, j) => { if (r[0] * ref[j][0] + r[1] * ref[j][1] < 0) { r[0] = -r[0]; r[1] = -r[1]; } });
  const lam0 = [Math.cos(theta0), Math.sin(theta0)];
  const roots = [evalC(p1, lam0), evalC(pi, lam0)];
  if (g < 2 || g % 2) return { g, p1, pi, roots, Bm, phi: null };
  if (!P) P = g === 2 ? [[1, 0], [0, 1]] : planeBasis(Bm);
  // b_l = x p_1 + y p_i with B-periods P_l (least squares: exact when the plane is the one P spans)
  const G11 = Bm.reduce((s, r) => s + r[0] * r[0], 0), G12 = Bm.reduce((s, r) => s + r[0] * r[1], 0);
  const G22 = Bm.reduce((s, r) => s + r[1] * r[1], 0), det = G11 * G22 - G12 * G12;
  const basis = P.map((Pl) => {
    const r0 = Bm.reduce((s, r, j) => s + r[0] * Pl[j], 0), r1 = Bm.reduce((s, r, j) => s + r[1] * Pl[j], 0);
    const x = (G22 * r0 - G12 * r1) / det, y = (G11 * r1 - G12 * r0) / det;
    return p1.map((c, k) => [x * c[0] + y * pi[k][0], x * c[1] + y * pi[k][1]]);
  });
  // q_l(y) = 1/2 int over theta0 .. theta0 + 2 pi of i b(lam) / nu dtheta, nu = sqrt(lam a) continued.
  // When a branch pair alpha, 1/conj(alpha) straddles the circle closely the integrand has a peak of
  // width ~ 1 - |alpha| at arg alpha, so the panels are graded geometrically towards each arg alpha_j.
  const a = aPoly(alphas);
  const acc = [0, 0];
  // the sheet at the start: on the circle nu = i lam^{(g+1)/2} prod |lam - alpha_j|, with lam^{(g+1)/2} taken
  // as e^{i theta0 (g+1)/2}. (The principal square root of lam0 a(lam0) jumps when that value crosses the
  // cut, which would flip the sign of phi.)
  const psi = (theta0 * (g + 1)) / 2;
  let prev = [-Math.sin(psi), Math.cos(psi)];
  for (const [t0, t1] of circlePanels(alphas, theta0)) {
    const c = (t0 + t1) / 2, r = (t1 - t0) / 2;
    for (let k = 0; k < GL.x.length; k++) {
      const t = c + r * GL.x[k], lam = [Math.cos(t), Math.sin(t)];
      const s = cmul(lam, evalF(a, lam));
      const m = Math.sqrt(Math.hypot(s[0], s[1])), ph = Math.atan2(s[1], s[0]) / 2;
      let nu = [m * Math.cos(ph), m * Math.sin(ph)];
      if (prev && (nu[0] - prev[0]) ** 2 + (nu[1] - prev[1]) ** 2 > (nu[0] + prev[0]) ** 2 + (nu[1] + prev[1]) ** 2) nu = [-nu[0], -nu[1]];
      prev = nu;
      const w = r * GL.w[k], n2 = nu[0] * nu[0] + nu[1] * nu[1];
      for (let l = 0; l < 2; l++) {
        const bb = evalC(basis[l], lam);
        acc[l] += w * (bb[0] * nu[0] + bb[1] * nu[1]) / n2; // Im(i b / nu) = Re(b / nu)
      }
    }
  }
  // latticeScale: the length scale of the periods in w = x + iy (p_w = x p_1 + y p_i): for g = 2 the
  // square root of the area of the lattice Bm^{-1} Z^2; in general det(Bm^T Bm)^(-1/4)
  const phi = g === 2 ? [Math.abs(acc[0] / 2), Math.abs(acc[1] / 2)] : [acc[0] / 2, acc[1] / 2];
  return { g, p1, pi, roots, Bm, P, basis, phi, latticeScale: Math.abs(det) ** -0.25 };
}

/** True when lam0 is (numerically) a common root of B_a, i.e. the data lie on S^g with Sym point lam0. */
export function onS(alphas, theta0, tol = 1e-6) {
  const d = symData(alphas, theta0);
  const scale = Math.max(...d.pi.map(([re, im]) => Math.hypot(re, im)), ...d.p1.map(([re, im]) => Math.hypot(re, im)));
  return Math.hypot(...d.roots[0], ...d.roots[1]) / scale < tol;
}

// Unknowns: per branch point (logit |alpha|, arg alpha), so that alpha -> 0 and |alpha| -> 1 are both at
// infinity and finite-difference steps are relative.
const R_LO = 1e-6, R_HI = 0.999;
// genus 2: phi fills the triangle phi_1, phi_2 > 0, phi_1 + phi_2 < pi (CKKS, Theorem 1.3)
const edgeDistance = ([p1, p2]) => Math.min(p1, p2, (Math.PI - p1 - p2) / Math.SQRT2);
const EDGE_MARGIN = 0.005;
const flat = (alphas) => alphas.flatMap(([re, im]) => { const r = Math.hypot(re, im); return [Math.log(r / (1 - r)), Math.atan2(im, re)]; });
const unflat = (x) => Array.from({ length: x.length / 2 }, (_, j) => {
  const r = 1 / (1 + Math.exp(-x[2 * j])), t = x[2 * j + 1];
  return [r * Math.cos(t), r * Math.sin(t)];
});

/** Solves the small dense system A x = b (A square) by Gaussian elimination with partial pivoting. */
function solveSq(A, b) {
  const n = b.length, M = A.map((r, i) => [...r, b[i]]);
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    [M[c], M[p]] = [M[p], M[c]];
    if (Math.abs(M[c][c]) < 1e-300) return null;
    for (let r = c + 1; r < n; r++) {
      const f = M[r][c] / M[c][c];
      for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k];
    }
  }
  const x = new Array(n).fill(0);
  for (let r = n - 1; r >= 0; r--) {
    let s = M[r][n];
    for (let k = r + 1; k < n; k++) s -= M[r][k] * x[k];
    x[r] = s / M[r][r];
  }
  return x;
}

/**
 * Follows the root-preserving flow from the current data (on S^2 with Sym point lam0 = e^{i theta0})
 * towards target values of the Sym integrals. The common root stays pinned at lam0, which also fixes
 * the rotation freedom.
 */
export class RootFlow {
  /** P: the period vectors of b_1, b_2 (default as in symData). */
  constructor(alphas, theta0, P = null) {
    const g = alphas.length;
    if (g < 2 || g % 2) throw new Error('The root-preserving flow is implemented for even genus');
    this.theta0 = theta0;
    this.x = flat(alphas);
    const d = symData(alphas, theta0, P);
    this.P = d.P;
    this.phi = d.phi;
    this.ref = d.Bm; // cycle orientations, carried along (see symData)
    // the period plane stays fixed: its normals, from the starting data
    this.normals = g === 2 ? [] : normals(planeBasis(d.Bm));
    this.edge = g === 2 ? edgeDistance : () => Infinity; // the image of phi is known for g = 2 only
  }

  get alphas() { return unflat(this.x); }

  /** Residuals: the common root (4 reals, scaled), the period plane (2(g - 2)) and phi - target (2). */
  residual(x, target) {
    const alphas = unflat(x);
    const rs = alphas.map(([re, im]) => Math.hypot(re, im));
    if (!rs.every((r) => r > R_LO && r < R_HI)) return null;
    for (let i = 0; i < alphas.length; i++) for (let j = i + 1; j < alphas.length; j++) {
      if (Math.hypot(alphas[i][0] - alphas[j][0], alphas[i][1] - alphas[j][1]) < 1e-3 * Math.max(rs[i], rs[j])) return null;
    }
    const d = symData(alphas, this.theta0, this.P, this.ref);
    // the root equations relative to the size of the differentials (which scale with alpha)
    const s = 1 / Math.max(...d.p1.map(([re, im]) => Math.hypot(re, im)));
    const F = [d.roots[0][0] * s, d.roots[0][1] * s, d.roots[1][0] * s, d.roots[1][1] * s];
    // both period vectors (normalised) orthogonal to the normals of the starting plane
    for (let l = 0; l < 2; l++) {
      const c = d.Bm.map((r) => r[l]), n = Math.hypot(...c);
      for (const N of this.normals) F.push(c.reduce((acc, v, j) => acc + v * N[j], 0) / n);
    }
    F.push(d.phi[0] - target[0], d.phi[1] - target[1]);
    return { F, phi: d.phi, Bm: d.Bm };
  }

  /** Gauss–Newton from x to the point with Sym integrals target. Returns the new x and phi, or null. */
  correct(x0, target, tol = 1e-10) {
    let x = x0.slice();
    let r = this.residual(x, target);
    if (!r) return null;
    for (let it = 0; it < 12; it++) {
      const n0 = Math.hypot(...r.F);
      if (n0 < tol) return { x, phi: r.phi, Bm: r.Bm };
      const eps = 1e-7, J = [];
      for (let i = 0; i < x.length; i++) {
        const xp = x.slice(); xp[i] += eps;
        const rp = this.residual(xp, target);
        if (!rp) return null;
        J.push(rp.F.map((v, k) => (v - r.F[k]) / eps)); // column i
      }
      // normal equations (J^T J) dx = -J^T F
      const n = x.length, A = [], b = [];
      for (let i = 0; i < n; i++) {
        A.push(Array.from({ length: n }, (_, j) => J[i].reduce((s, v, k) => s + v * J[j][k], 0)));
        b.push(-J[i].reduce((s, v, k) => s + v * r.F[k], 0));
      }
      const dx = solveSq(A, b);
      if (!dx || !dx.every(Number.isFinite)) return null;
      // damped step: halve until the residual decreases
      let lam = 1, next = null;
      for (let k = 0; k < 8; k++, lam /= 2) {
        const xt = x.map((v, i) => v + lam * dx[i]);
        const rt = this.residual(xt, target);
        if (rt && Math.hypot(...rt.F) < n0) { next = { x: xt, r: rt }; break; }
      }
      if (!next) return Math.hypot(...r.F) < 1e-7 ? { x, phi: r.phi, Bm: r.Bm } : null;
      x = next.x; r = next.r;
    }
    return Math.hypot(...r.F) < 1e-7 ? { x, phi: r.phi, Bm: r.Bm } : null;
  }

  /**
   * Moves towards target along the straight segment in phi, in steps of at most `step`, for at most
   * `budget` ms. Returns { alphas, phi, reached, blocked, why }: blocked when the flow can't continue,
   * why = 'circle' | 'zero' | 'collision' (a branch point near the unit circle or 0, or two meeting:
   * an end of S^g) or 'fold' (the solver fails in the interior: phi folds over, and stops being a chart).
   */
  moveToward(target, { step = 0.03, budget = 12 } = {}) {
    const t0 = performance.now();
    // the step size carries over between calls (one call per animation frame) while the target is the same
    if (!this.lastTarget || this.lastTarget[0] !== target[0] || this.lastTarget[1] !== target[1]) {
      this.lastTarget = target.slice();
      this.h = step;
    }
    let h = this.h, blocked = false;
    while (performance.now() - t0 < budget) {
      const dphi = [target[0] - this.phi[0], target[1] - this.phi[1]];
      const dist = Math.hypot(dphi[0], dphi[1]);
      if (dist < 1e-9) break;
      // Near the edge of the triangle the chart phi degenerates (a small change of phi needs a large
      // change of alpha), so steps are kept below 0.3 x the distance to the edge, and the flow stops
      // EDGE_MARGIN short of it, from where it can always step back inside.
      const e0 = this.edge(this.phi), cap = Math.max(1e-7, Math.min(step, 0.3 * e0));
      h = Math.min(h, cap);
      const f = Math.min(1, h / dist);
      const want = [this.phi[0] + f * dphi[0], this.phi[1] + f * dphi[1]];
      const ew = this.edge(want);
      if (ew < EDGE_MARGIN && ew < e0) { blocked = true; break; }
      const res = this.correct(this.x, want);
      // reject jumps: a converged point far from the start is on another branch
      if (res && Math.hypot(...res.x.map((v, i) => v - this.x[i])) < 0.5) {
        this.x = res.x;
        this.phi = res.phi;
        this.ref = res.Bm;
        h = Math.min(step, h * 2);
      } else if (h > cap / 1024) {
        h /= 4;
      } else {
        blocked = true; // the solver can't continue (an end of S^2 inside the margin, or a failure)
        break;
      }
    }
    this.h = h;
    if (blocked) this.lastTarget = null; // a new request from here starts afresh
    const reached = Math.hypot(target[0] - this.phi[0], target[1] - this.phi[1]) < 1e-9;
    return { alphas: this.alphas, phi: this.phi.slice(), reached, blocked, why: blocked ? this.why() : null };
  }

  /** Why the flow stopped at the current point (see moveToward). */
  why() {
    if (this.alphas.length === 2 && this.edge(this.phi) < 2 * EDGE_MARGIN) return 'edge';
    const al = this.alphas, rs = al.map(([re, im]) => Math.hypot(re, im));
    if (Math.max(...rs) > 0.95) return 'circle';
    if (Math.min(...rs) < 0.02) return 'zero';
    for (let i = 0; i < al.length; i++) for (let j = i + 1; j < al.length; j++) {
      if (Math.hypot(al[i][0] - al[j][0], al[i][1] - al[j][1]) < 0.02 * Math.max(rs[i], rs[j])) return 'collision';
    }
    return 'fold';
  }

}
