// CMC (H = 1/2) immersions from spectral data via polynomial Killing fields.
// Section numbers refer to cmc_spectral_to_immersion.md; cmc_reference.py is the ground truth.
//
// Storage. A 2x2 complex matrix occupies 8 consecutive doubles
//   [m11.re, m11.im, m12.re, m12.im, m21.re, m21.im, m22.re, m22.im].
// A Laurent polynomial sum_{p=lo}^{hi} X_p lam^p is a Float64Array of (hi-lo+1) such slots,
// slot i holding X_{lo+i}. The Killing field zeta has lo = -1, hi = g, so n = g+2 slots.

import { polyMul, polyEval, polySharp, polyRoots } from './poly.js';

// ---------------------------------------------------------------- 2x2 complex matrix kernels

/** C[co] += s * A[ao] B[bo], s = +1 or -1. */
function mulAcc(A, ao, B, bo, C, co, s) {
  const a0r = A[ao], a0i = A[ao + 1], a1r = A[ao + 2], a1i = A[ao + 3];
  const a2r = A[ao + 4], a2i = A[ao + 5], a3r = A[ao + 6], a3i = A[ao + 7];
  const b0r = B[bo], b0i = B[bo + 1], b1r = B[bo + 2], b1i = B[bo + 3];
  const b2r = B[bo + 4], b2i = B[bo + 5], b3r = B[bo + 6], b3i = B[bo + 7];
  C[co]     += s * (a0r * b0r - a0i * b0i + a1r * b2r - a1i * b2i);
  C[co + 1] += s * (a0r * b0i + a0i * b0r + a1r * b2i + a1i * b2r);
  C[co + 2] += s * (a0r * b1r - a0i * b1i + a1r * b3r - a1i * b3i);
  C[co + 3] += s * (a0r * b1i + a0i * b1r + a1r * b3i + a1i * b3r);
  C[co + 4] += s * (a2r * b0r - a2i * b0i + a3r * b2r - a3i * b2i);
  C[co + 5] += s * (a2r * b0i + a2i * b0r + a3r * b2i + a3i * b2r);
  C[co + 6] += s * (a2r * b1r - a2i * b1i + a3r * b3r - a3i * b3i);
  C[co + 7] += s * (a2r * b1i + a2i * b1r + a3r * b3i + a3i * b3r);
}

// transpose permutation of the four entries: (r,s) -> (s,r)
const TR = [0, 2, 1, 3];

// ---------------------------------------------------------------- spectral data (§2)

/** a(lam) = -prod (lam - alpha)(1 - conj(alpha) lam), coefficients by increasing power. */
export function aPoly(alphas) {
  let a = Float64Array.of(-1, 0);
  for (const [ar, ai] of alphas) {
    a = polyMul(a, Float64Array.of(-ar, -ai, 1 + ar * ar + ai * ai, 0, -ar, ai));
  }
  return a;
}

/** kappa_0 = 2 sqrt|a_0| = 2 prod |alpha_i|^(1/2). */
export function kappa0(alphas) {
  let p = 1;
  for (const [ar, ai] of alphas) p *= Math.hypot(ar, ai);
  return 2 * Math.sqrt(p);
}

// ---------------------------------------------------------------- initial Killing field (§4)

/**
 * Polynomial Killing field with diagonal part A (complex coefficients [re,im] of powers 0..g-1,
 * satisfying A_k = -conj(A_{g-1-k})), via Fejér–Riesz. A = null gives the canonical base point.
 * Returns a Float64Array of g+2 matrix slots.
 */
export function killingField(alphas, A = null) {
  const g = alphas.length;
  const a = aPoly(alphas);
  const Ac = new Float64Array(2 * g);
  if (A) A.forEach(([re, im], k) => { Ac[2 * k] = re; Ac[2 * k + 1] = im; });
  let b;
  if (g === 0) {
    b = Float64Array.of(1, 0);
  } else {
    // q = lam A^2 - a
    const q = a.map((x) => -x);
    const A2 = polyMul(Ac, Ac);
    for (let k = 0; k < A2.length; k++) q[k + 2] += A2[k];
    const rin = polyRoots(q).filter(([re, im]) => Math.hypot(re, im) < 1);
    if (rin.length !== g) throw new Error('A violates |A| <= prod|lam - alpha| on the unit circle');
    b = Float64Array.of(1, 0);
    for (const [re, im] of rin) b = polyMul(b, Float64Array.of(-re, -im, 1, 0));
    // |beta|^2 = q_{2g} / (b_g conj(b_0)), real and positive
    const [qr, qi] = [q[4 * g], q[4 * g + 1]];
    const [dr, di] = [b[2 * g] * b[0] + b[2 * g + 1] * b[1], b[2 * g + 1] * b[0] - b[2 * g] * b[1]];
    const s = (qr * dr + qi * di) / (dr * dr + di * di);
    // gauge: b(0) > 0
    const ph = Math.atan2(b[1], b[0]);
    const sr = Math.sqrt(s) * Math.cos(ph), si = -Math.sqrt(s) * Math.sin(ph);
    for (let k = 0; k <= g; k++) {
      const re = b[2 * k], im = b[2 * k + 1];
      b[2 * k] = re * sr - im * si;
      b[2 * k + 1] = re * si + im * sr;
    }
  }
  const n = g + 2;
  const Z = new Float64Array(8 * n);
  for (let k = 0; k < g; k++) {
    const o = 8 * (k + 1);
    Z[o] = Ac[2 * k]; Z[o + 1] = Ac[2 * k + 1];
    Z[o + 6] = -Ac[2 * k]; Z[o + 7] = -Ac[2 * k + 1];
  }
  for (let k = 0; k <= g; k++) { Z[8 * k + 2] = b[2 * k]; Z[8 * k + 3] = b[2 * k + 1]; }
  const bs = polySharp(b, g);
  for (let k = 0; k <= g; k++) { Z[8 * (k + 1) + 4] = -bs[2 * k]; Z[8 * (k + 1) + 5] = -bs[2 * k + 1]; }
  return Z;
}

/**
 * Eigenline divisor readout (§4): points (mu_j, nu_j) where the eigenline of zeta is [0:1].
 * Returns [{mu:[re,im], nu:[re,im]}].
 */
export function divisor(Z) {
  const g = Z.length / 8 - 2;
  if (g === 0) return [];
  const b = new Float64Array(2 * (g + 1));
  for (let k = 0; k <= g; k++) { b[2 * k] = Z[8 * k + 2]; b[2 * k + 1] = Z[8 * k + 3]; }
  const A = new Float64Array(2 * g);
  for (let k = 0; k < g; k++) { A[2 * k] = Z[8 * (k + 1)]; A[2 * k + 1] = Z[8 * (k + 1) + 1]; }
  return polyRoots(b).map(([mr, mi]) => {
    const [Ar, Ai] = polyEval(A, mr, mi);
    return { mu: [mr, mi], nu: [-(mr * Ar - mi * Ai), -(mr * Ai + mi * Ar)] };
  });
}

// ---------------------------------------------------------------- flows (§§3, 5)

/**
 * G_{k,c} = c W_k + R(c W_k), written into G (2k+3 slots, powers -1-k .. 1+k).
 * Z has n slots; only its first k+2 are read.
 */
export function generator(Z, kap, k, cr, ci, G) {
  G.fill(0);
  const m = 2 * k + 2;
  const w = new Float64Array(8);
  for (let i = 0; i <= k + 1; i++) {
    for (let e = 0; e < 8; e++) w[e] = Z[8 * i + e] / kap;
    if (i === k + 1) { // pi_0 projection
      w[0] *= 0.5; w[1] *= 0.5; w[2] = 0; w[3] = 0; w[6] *= 0.5; w[7] *= 0.5;
    }
    const o = 8 * i, oR = 8 * (m - i);
    for (let e = 0; e < 4; e++) {
      const wr = w[2 * e], wi = w[2 * e + 1];
      G[o + 2 * e] += cr * wr - ci * wi;
      G[o + 2 * e + 1] += cr * wi + ci * wr;
      // R(cW)_{rs} = -conj(c W_{sr}) at the reflected power
      const t = TR[e], vr = w[2 * t], vi = w[2 * t + 1];
      G[oR + 2 * e] -= cr * vr - ci * vi;
      G[oR + 2 * e + 1] += cr * vi + ci * vr;
    }
  }
}

/** out (n slots) = [Z, G] truncated to powers -1..g; G has 2k+3 slots starting at power -1-k. */
export function commutator(Z, n, G, k, out) {
  out.fill(0, 0, 8 * n);
  const nG = 2 * k + 3;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < nG; j++) {
      const t = i + j - (k + 1);
      if (t < 0 || t >= n) continue;
      mulAcc(Z, 8 * i, G, 8 * j, out, 8 * t, 1);
      mulAcc(G, 8 * j, Z, 8 * i, out, 8 * t, -1);
    }
  }
}

/** Re-impose the degree structure, trace-freeness and reality zeta_k = -zeta_{g-1-k}^dagger. */
export function projectKilling(Z) {
  const n = Z.length / 8, g = n - 2;
  // trace-free
  for (let i = 0; i < n; i++) {
    const o = 8 * i;
    const hr = 0.5 * (Z[o] - Z[o + 6]), hi = 0.5 * (Z[o + 1] - Z[o + 7]);
    Z[o] = hr; Z[o + 1] = hi; Z[o + 6] = -hr; Z[o + 7] = -hi;
  }
  // nilpotent ends: zeta_{-1} = B E12, zeta_g = C E21
  Z[0] = Z[1] = Z[4] = Z[5] = Z[6] = Z[7] = 0;
  const e = 8 * (n - 1);
  Z[e] = Z[e + 1] = Z[e + 2] = Z[e + 3] = Z[e + 6] = Z[e + 7] = 0;
  // reality: slot i pairs with slot g+1-i
  for (let i = 0; 2 * i <= g + 1; i++) {
    const j = g + 1 - i, oi = 8 * i, oj = 8 * j;
    for (let r = 0; r < 4; r++) {
      const t = TR[r];
      // new_i[r] = (Z_i[r] - conj(Z_j[t])) / 2
      w4[2 * r] = 0.5 * (Z[oi + 2 * r] - Z[oj + 2 * t]);
      w4[2 * r + 1] = 0.5 * (Z[oi + 2 * r + 1] + Z[oj + 2 * t + 1]);
    }
    for (let r = 0; r < 8; r++) Z[oi + r] = w4[r];
    if (i !== j) {
      // new_j = -new_i^dagger
      for (let r = 0; r < 4; r++) {
        const t = TR[r];
        Z[oj + 2 * r] = -w4[2 * t];
        Z[oj + 2 * r + 1] = w4[2 * t + 1];
      }
    }
  }
}
const w4 = new Float64Array(8);

/** -det zeta(lam) - a(lam)/lam, the spectral-curve invariant defect at one lam. */
export function detDefect(Z, a, lr, li) {
  const n = Z.length / 8;
  // evaluate zeta(lam)
  const d2 = lr * lr + li * li;
  let pr = lr / d2, pi = -li / d2; // lam^p, starting at p = -1
  const M = [0, 0, 0, 0, 0, 0, 0, 0];
  for (let i = 0; i < n; i++) {
    for (let e = 0; e < 4; e++) {
      const zr = Z[8 * i + 2 * e], zi = Z[8 * i + 2 * e + 1];
      M[2 * e] += zr * pr - zi * pi;
      M[2 * e + 1] += zr * pi + zi * pr;
    }
    const t = pr * lr - pi * li; pi = pr * li + pi * lr; pr = t;
  }
  // -det = -(m11 m22 - m12 m21)
  const dr = -(M[0] * M[6] - M[1] * M[7] - (M[2] * M[4] - M[3] * M[5]));
  const di = -(M[0] * M[7] + M[1] * M[6] - (M[2] * M[5] + M[3] * M[4]));
  const [ar, ai] = polyEval(a, lr, li);
  const qr = (ar * lr + ai * li) / d2, qi = (ai * lr - ar * li) / d2;
  return Math.hypot(dr - qr, di - qi);
}

/** Isospectral flow (§5): integrate d zeta/dt = [zeta, G_{k,c}] for time t (in place). */
export function flowKilling(Z, kap, k, c, t, hmax = 0.01) {
  const n = Z.length / 8;
  const G = new Float64Array(8 * (2 * k + 3));
  const k1 = new Float64Array(8 * n), k2 = new Float64Array(8 * n), k3 = new Float64Array(8 * n),
    k4 = new Float64Array(8 * n), y = new Float64Array(8 * n);
  const rhs = (Y, out) => { generator(Y, kap, k, c[0], c[1], G); commutator(Y, n, G, k, out); };
  const N = Math.max(1, Math.ceil(Math.abs(t) / hmax));
  const h = t / N;
  for (let s = 0; s < N; s++) {
    rhs(Z, k1);
    for (let e = 0; e < 8 * n; e++) y[e] = Z[e] + 0.5 * h * k1[e];
    rhs(y, k2);
    for (let e = 0; e < 8 * n; e++) y[e] = Z[e] + 0.5 * h * k2[e];
    rhs(y, k3);
    for (let e = 0; e < 8 * n; e++) y[e] = Z[e] + h * k3[e];
    rhs(y, k4);
    for (let e = 0; e < 8 * n; e++) Z[e] += (h / 6) * (k1[e] + 2 * k2[e] + 2 * k3[e] + k4[e]);
  }
  projectKilling(Z);
  return Z;
}

/** Independent real isospectral flows (k, c) (§5). The first two (k = 0) are the x, y translations. */
export function isospectralBasis(g) {
  const B = [];
  for (let k = 0; k < Math.floor(g / 2); k++) B.push([k, [1, 0]], [k, [0, 1]]);
  if (g % 2 === 1) B.push([(g - 1) / 2, [0, 1]]);
  return B;
}

/** Shape-changing flows: the basis without the two k = 0 translations (for g >= 2). */
export function shapeFlows(g) {
  const B = isospectralBasis(g);
  return g >= 2 ? B.slice(2) : [];
}

// ---------------------------------------------------------------- frame integration (§§6–7)

/**
 * Integrator for the state (zeta, F, G = dF/dlam) at the Sym point lam0 = e^{i theta0},
 * stepping in the z-direction c = e^{i phi}.
 */
export class FrameIntegrator {
  constructor(alphas, theta0) {
    this.g = alphas.length;
    this.n = this.g + 2;
    this.a = aPoly(alphas);
    this.kap = kappa0(alphas);
    this.l0r = Math.cos(theta0); this.l0i = Math.sin(theta0);
    const len = 8 * this.n + 16;
    this.len = len;
    this.Gen = new Float64Array(24);
    this.X = new Float64Array(8);
    this.Xd = new Float64Array(8);
    this.k1 = new Float64Array(len); this.k2 = new Float64Array(len);
    this.k3 = new Float64Array(len); this.k4 = new Float64Array(len);
    this.tmp = new Float64Array(len);
  }

  /** Initial state from a Killing field: F = I, G = 0. */
  initialState(Z) {
    const y = new Float64Array(this.len);
    y.set(Z);
    const f = 8 * this.n;
    y[f] = 1; y[f + 6] = 1;
    return y;
  }

  rhs(y, dy, cr, ci) {
    const n = this.n, Gen = this.Gen, X = this.X, Xd = this.Xd;
    generator(y, this.kap, 0, cr, ci, Gen);
    commutator(y, n, Gen, 0, dy);
    // X(lam0) = G_{-1}/lam0 + G_0 + G_1 lam0,  X'(lam0) = -G_{-1}/lam0^2 + G_1  (|lam0| = 1)
    const lr = this.l0r, li = this.l0i;
    const ir = lr, ii = -li;                          // 1/lam0
    const i2r = ir * ir - ii * ii, i2i = 2 * ir * ii; // 1/lam0^2
    for (let e = 0; e < 4; e++) {
      const ar = Gen[2 * e], ai = Gen[2 * e + 1];
      const br = Gen[8 + 2 * e], bi = Gen[8 + 2 * e + 1];
      const cr2 = Gen[16 + 2 * e], ci2 = Gen[16 + 2 * e + 1];
      X[2 * e] = ar * ir - ai * ii + br + cr2 * lr - ci2 * li;
      X[2 * e + 1] = ar * ii + ai * ir + bi + cr2 * li + ci2 * lr;
      Xd[2 * e] = -(ar * i2r - ai * i2i) + cr2;
      Xd[2 * e + 1] = -(ar * i2i + ai * i2r) + ci2;
    }
    const f = 8 * n, gg = f + 8;
    for (let e = 0; e < 16; e++) dy[f + e] = 0;
    mulAcc(y, f, X, 0, dy, f, 1);
    mulAcc(y, gg, X, 0, dy, gg, 1);
    mulAcc(y, f, Xd, 0, dy, gg, 1);
  }

  /** Advance y in place by distance t along direction c = (cr, ci), |c| = 1, with RK4 steps <= hmax. */
  march(y, cr, ci, t, hmax) {
    if (t === 0) return;
    const N = Math.max(1, Math.ceil(Math.abs(t) / hmax));
    const h = t / N, len = this.len;
    const { k1, k2, k3, k4, tmp } = this;
    for (let s = 0; s < N; s++) {
      this.rhs(y, k1, cr, ci);
      for (let e = 0; e < len; e++) tmp[e] = y[e] + 0.5 * h * k1[e];
      this.rhs(tmp, k2, cr, ci);
      for (let e = 0; e < len; e++) tmp[e] = y[e] + 0.5 * h * k2[e];
      this.rhs(tmp, k3, cr, ci);
      for (let e = 0; e < len; e++) tmp[e] = y[e] + h * k3[e];
      this.rhs(tmp, k4, cr, ci);
      for (let e = 0; e < len; e++) y[e] += (h / 6) * (k1[e] + 2 * k2[e] + 2 * k3[e] + k4[e]);
    }
  }

  /** Re-impose the Killing-field structure and F in SU(2). */
  project(y) {
    const zeta = y.subarray(0, 8 * this.n);
    projectKilling(zeta);
    const f = 8 * this.n;
    // F = [[p, q], [-conj q, conj p]]
    let pr = 0.5 * (y[f] + y[f + 6]), pi = 0.5 * (y[f + 1] - y[f + 7]);
    let qr = 0.5 * (y[f + 2] - y[f + 4]), qi = 0.5 * (y[f + 3] + y[f + 5]);
    const s = 1 / Math.sqrt(pr * pr + pi * pi + qr * qr + qi * qi);
    pr *= s; pi *= s; qr *= s; qi *= s;
    y[f] = pr; y[f + 1] = pi; y[f + 2] = qr; y[f + 3] = qi;
    y[f + 4] = -qr; y[f + 5] = qi; y[f + 6] = pr; y[f + 7] = -pi;
  }

  /**
   * Geometry at state y: writes f (3), N (3) into out at offsets and returns u.
   * f = 4 i lam0 G F^{-1}, N = -F e3 F^{-1}, e^u = 2|B_{-1}|/kappa0.
   */
  evaluate(y, P, po, Nn, no) {
    const f = 8 * this.n, gg = f + 8;
    // F^{-1} = [[d, -b], [-c, a]] (det F = 1)
    const Fi = this.tmp;
    Fi[0] = y[f + 6]; Fi[1] = y[f + 7]; Fi[2] = -y[f + 2]; Fi[3] = -y[f + 3];
    Fi[4] = -y[f + 4]; Fi[5] = -y[f + 5]; Fi[6] = y[f]; Fi[7] = y[f + 1];
    const M = this.k1;
    for (let e = 0; e < 8; e++) M[e] = 0;
    mulAcc(y, gg, Fi, 0, M, 0, 1); // G F^{-1}
    // X = 4 i lam0 M; coordinates x1 = Im(x12 + x21)/2, x2 = Re(x12 - x21)/2, x3 = Im(x11 - x22)/2
    const sr = -4 * this.l0i, si = 4 * this.l0r; // 4 i lam0
    const xr = (e) => sr * M[2 * e] - si * M[2 * e + 1];
    const xi = (e) => sr * M[2 * e + 1] + si * M[2 * e];
    P[po] = 0.5 * (xi(1) + xi(2));
    P[po + 1] = 0.5 * (xr(1) - xr(2));
    P[po + 2] = 0.5 * (xi(0) - xi(3));
    // N = -F e3 F^{-1}; F e3 = [[i F11, -i F12], [i F21, -i F22]]
    const Fe = this.k2;
    Fe[0] = -y[f + 1]; Fe[1] = y[f];
    Fe[2] = y[f + 3]; Fe[3] = -y[f + 2];
    Fe[4] = -y[f + 5]; Fe[5] = y[f + 4];
    Fe[6] = y[f + 7]; Fe[7] = -y[f + 6];
    for (let e = 0; e < 8; e++) M[e] = 0;
    mulAcc(Fe, 0, Fi, 0, M, 0, -1);
    Nn[no] = 0.5 * (M[3] + M[5]);
    Nn[no + 1] = 0.5 * (M[2] - M[4]);
    Nn[no + 2] = 0.5 * (M[1] - M[7]);
    return Math.log((2 * Math.hypot(y[2], y[3])) / this.kap);
  }
}

// ---------------------------------------------------------------- the whole pipeline

/** Hopf differential Q = -e^{i arg a0}/lam0 (§6), returned as its argument. */
export function hopfArg(alphas, theta0) {
  const a = aPoly(alphas);
  return Math.PI + Math.atan2(a[1], a[0]) - theta0;
}

/**
 * Compute the immersion on an nx x ny grid.
 *
 * params:
 *   alphas   [[re, im], ...] branch points in the punctured unit disc
 *   theta0   Sym point lam0 = e^{i theta0}
 *   tau      shape times, one per shapeFlows(g)
 *   z0       [x, y] domain centre (the frame is transported from z = 0, so moving z0 slides the
 *            window over a fixed surface)
 *   phi      grid rotation: z = z0 + e^{i phi} (s + i t); null means curvature-line aligned
 *   width, height   extent in s and t
 *   nx, ny   grid resolution
 *   hmax     RK4 step bound
 *
 * Returns { nx, ny, pos, nrm (Float32Array 3N), u, st (Float32Array N / 2N), bbox, divisor, stats }.
 * Vertex (i, j) has index j*nx + i.
 */
export function computeSurface(params) {
  const t0 = (typeof performance !== 'undefined' ? performance : Date).now();
  const { alphas, theta0, tau = [], z0 = [0, 0], width, height, nx, ny, hmax = 0.02 } = params;
  const g = alphas.length;
  const kap = kappa0(alphas);
  const phi = params.phi ?? -0.5 * hopfArg(alphas, theta0);

  // zeta^init: base point moved along the shape flows
  const Z = killingField(alphas);
  shapeFlows(g).forEach(([k, c], i) => {
    if (tau[i]) flowKilling(Z, kap, k, c, tau[i], 0.002); // cheap: done once, tiny system
  });

  const I = new FrameIntegrator(alphas, theta0);
  const y0 = I.initialState(Z);
  const r0 = Math.hypot(z0[0], z0[1]);
  if (r0 > 0) I.march(y0, z0[0] / r0, z0[1] / r0, r0, hmax);
  I.project(y0);

  const N = nx * ny;
  const pos = new Float32Array(3 * N), nrm = new Float32Array(3 * N);
  const u = new Float32Array(N), st = new Float32Array(2 * N);
  const s = Array.from({ length: nx }, (_, i) => -width / 2 + (nx > 1 ? (i * width) / (nx - 1) : 0));
  const t = Array.from({ length: ny }, (_, j) => -height / 2 + (ny > 1 ? (j * height) / (ny - 1) : 0));
  const csr = Math.cos(phi), csi = Math.sin(phi); // s-direction
  const ctr = -csi, cti = csr;                   // t-direction = i e^{i phi}
  const P = new Float64Array(3), Nv = new Float64Array(3);
  let maxDet = 0;
  const lam = [Math.cos(theta0 + 0.3), Math.sin(theta0 + 0.3)];

  const sweep = (ystart, coords, cr, ci, visit) => {
    for (const dir of [1, -1]) {
      const y = ystart.slice();
      let at = 0;
      const idx = coords.map((_, i) => i).filter((i) => (dir > 0 ? coords[i] >= 0 : coords[i] < 0));
      if (dir < 0) idx.reverse();
      for (const i of idx) {
        I.march(y, cr, ci, coords[i] - at, hmax);
        at = coords[i];
        I.project(y);
        visit(i, y);
      }
    }
  };

  sweep(y0, s, csr, csi, (i, ys) => {
    sweep(ys, t, ctr, cti, (j, yt) => {
      const v = j * nx + i;
      u[v] = I.evaluate(yt, P, 0, Nv, 0);
      pos[3 * v] = P[0]; pos[3 * v + 1] = P[1]; pos[3 * v + 2] = P[2];
      nrm[3 * v] = Nv[0]; nrm[3 * v + 1] = Nv[1]; nrm[3 * v + 2] = Nv[2];
      st[2 * v] = s[i]; st[2 * v + 1] = t[j];
    });
    maxDet = Math.max(maxDet, detDefect(ys.subarray(0, 8 * (g + 2)), I.a, lam[0], lam[1]));
  });

  const bbox = [Infinity, Infinity, Infinity, -Infinity, -Infinity, -Infinity];
  for (let v = 0; v < N; v++) {
    for (let e = 0; e < 3; e++) {
      const x = pos[3 * v + e];
      if (Number.isFinite(x)) {
        if (x < bbox[e]) bbox[e] = x;
        if (x > bbox[e + 3]) bbox[e + 3] = x;
      }
    }
  }
  const t1 = (typeof performance !== 'undefined' ? performance : Date).now();
  return {
    nx, ny, pos, nrm, u, st, bbox, phi,
    divisor: divisor(y0.subarray(0, 8 * (g + 2))),
    stats: { ms: t1 - t0, detDefect: maxDet },
  };
}
