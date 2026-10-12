// Genus 2: the boundary of S^2 in F_3, from Hasse–Ogilvie–Schmidt, "... constrained Willmore tori II"
// (the genus 2 paper). S^2 (where B_a has a common root on the unit circle) meets the two ends of the
// τ-fibres in two families of curves:
//
//  - at the nodal end (a double root of a on S^1, normalisation of genus 1 with modulus k): the double
//    point at the root t_0(k) of df, t_0 = sqrt(E / (E - (1 - k^2) K)) (the theorem on the low genus
//    boundary of S^2), with
//        Im τ = (2K/π) sqrt((1 + k^2 t_0^2) / (1 + t_0^2)),
//        Re τ = (2/π) ∫_0^{t_0} ((E - (1-k^2)K) s^2 + E) / ((1 + s^2)^{3/2} (1 + k^2 s^2)^{1/2}) ds;
//    two curves, τ and -conj(τ), from ±1 (k -> 0) to ±1/2 + i∞ (k -> 1).
//  - at the open end (all branch points to 0, rescaled roots α, 1/α): the real locus of
//    f_0 = i (K(α) - E(α)) / E(α'), α' = sqrt(1 - α^2) (the lemma on the real locus), with
//    τ = τ_0 / (1 + τ_0), τ_0 = i K(α) / K(α'). (τ_0 / (1 + τ_0) is the ratio for the basis dual to the
//    cycles B_1^+, B_2^+; the differentials written in the paper have their periods on B_1^+ and
//    B_1^+ - B_2^+.) One curve through the unit circle, from i∞ to i∞, symmetric under τ -> -conj(τ).
//
// Checked against the viewer's own periods (test/run.mjs): opening the node slightly (|alpha| = 1 - 1e-4)
// gives the first τ to 5 digits; scaling the roots s alpha, s/alpha with s = 1e-3 gives the second to 1e-3
// (converging as s -> 0), and there the ratio of the linear coefficients of the lattice basis is real.

const mul = ([a, b], [c, d]) => [a * c - b * d, a * d + b * c];
const div = ([a, b], [c, d]) => { const q = c * c + d * d; return [(a * c + b * d) / q, (b * c - a * d) / q]; };
const add = ([a, b], [c, d]) => [a + c, b + d];
const sub = ([a, b], [c, d]) => [a - c, b - d];
const abs = ([a, b]) => Math.hypot(a, b);
const csqrt = ([a, b]) => { const r = Math.sqrt(Math.hypot(a, b)), t = Math.atan2(b, a) / 2; return [r * Math.cos(t), r * Math.sin(t)]; };

/**
 * Complete elliptic integrals K(k), E(k) of the first and second kind for complex modulus k (|k| < 1 or
 * off the cuts (-inf, -1], [1, inf)), by the arithmetic–geometric mean with the right choice of square
 * roots: K = π / (2 M(1, k')), E = K (1 - Σ_{n>=0} 2^{n-1} c_n^2), c_0 = k, c_{n+1} = (a_n - b_n) / 2.
 */
export function ellipticKE(k) {
  let a = [1, 0], b = csqrt(sub([1, 0], mul(k, k))), sum = mul([0.5, 0], mul(k, k)), w = 1;
  for (let n = 0; n < 60; n++) {
    const c = mul([0.5, 0], sub(a, b));
    const an = mul([0.5, 0], add(a, b));
    let bn = csqrt(mul(a, b));
    if (abs(sub(an, bn)) > abs(add(an, bn))) bn = [-bn[0], -bn[1]];
    sum = add(sum, mul([w, 0], mul(c, c)));
    w *= 2;
    a = an; b = bn;
    if (abs(c) < 1e-17) break;
  }
  const K = div([Math.PI / 2, 0], a);
  return { K, E: mul(K, sub([1, 0], sum)) };
}

/** τ reduced into F_3 = { |Re τ| <= 1, |τ| >= 1 } by τ -> τ + 2 and τ -> -1/τ. */
export function toF3(t) {
  for (let i = 0; i < 200; i++) {
    const [x, y] = t;
    if (Math.abs(x) > 1) { t = [x - 2 * Math.round(x / 2), y]; continue; }
    if (x * x + y * y < 1 - 1e-13) { t = div([-1, 0], t); continue; }
    break;
  }
  return t;
}

/** The nodal-end boundary of S^2 at modulus k in (0, 1): { t0, tau } (the curve with Re τ > 0). */
export function nodalBoundary(k) {
  const { K: [K], E: [E] } = ellipticKE([k, 0]);
  const D = E - (1 - k * k) * K, t0 = Math.sqrt(E / D);
  const f = (s) => (D * s * s + E) / ((1 + s * s) ** 1.5 * Math.sqrt(1 + k * k * s * s));
  // Simpson in u = atan(s), which maps [0, t0] to [0, atan t0] and tames the tail as k -> 0
  const U = Math.atan(t0), N = 400;
  let re = 0;
  for (let i = 0; i <= N; i++) {
    const u = (U * i) / N, s = Math.tan(u), wt = i === 0 || i === N ? 1 : i % 2 ? 4 : 2;
    re += (wt * f(s) * (1 + s * s) * U) / (3 * N);
  }
  re *= 2 / Math.PI;
  const im = ((2 * K) / Math.PI) * Math.sqrt((1 + k * k * t0 * t0) / (1 + t0 * t0));
  return { t0, tau: [re, im] };
}

/** f_0(α) = i (K(α) - E(α)) / E(α'), whose real locus is the open-end boundary of S^2. */
export function f0(al) {
  const { K, E } = ellipticKE(al), Ep = ellipticKE(csqrt(sub([1, 0], mul(al, al)))).E;
  return mul([0, 1], div(sub(K, E), Ep));
}

/** τ at the open end with rescaled roots α, 1/α (reduced into F_3). */
export function openEndTau(al) {
  const t0 = mul([0, 1], div(ellipticKE(al).K, ellipticKE(csqrt(sub([1, 0], mul(al, al)))).K));
  return toF3(div(t0, add([1, 0], t0)));
}

/**
 * Polylines in F_3: { nodal: [curve, mirror], open: [curve, mirror] }, each an array of τ. Computed once.
 * Each open-end piece runs from the unit circle (where the two pieces meet, as identified points of F_3)
 * up to i∞.
 * The open-end locus is found on circles |α| = ρ, ρ in (0, 1], by bisection in arg α ∈ (0, π/2) (in the
 * first quadrant it is a simple curve from 0 to the unit circle; the part beyond the circle is its image
 * under α -> 1/conj(α), which maps to the mirror image in F_3).
 */
let CURVES = null;
export function s2BoundaryCurves() {
  if (CURVES) return CURVES;
  const mirror = (pts) => pts.map(([x, y]) => [-x, y]);
  const nodal = [];
  for (let j = 0; j <= 120; j++) {
    const k = j < 60 ? 0.002 + (0.5 - 0.002) * (j / 60) ** 2 : 1 - 0.5 * 0.9 ** (j - 60);
    nodal.push(toF3(nodalBoundary(k).tau));
  }
  const open = [];
  for (let j = 0; j <= 90; j++) {
    const rho = 0.92 ** j; // from the unit circle towards 0
    const g = (th) => f0([rho * Math.cos(th), rho * Math.sin(th)])[1];
    let lo = 1e-6, hi = Math.PI / 2 - 1e-6;
    const glo = g(lo);
    if (!(glo * g(hi) < 0)) continue;
    for (let i = 0; i < 60; i++) { const m = (lo + hi) / 2; if (Math.sign(g(m)) === Math.sign(glo)) lo = m; else hi = m; }
    const th = (lo + hi) / 2;
    open.push(openEndTau([rho * Math.cos(th), rho * Math.sin(th)]));
  }
  // on the unit circle τ and -conj(τ) are the same point of F_3: take the representative continuous with
  // the neighbouring point (else the polyline would cut across the top of the circle)
  for (const pts of [nodal, open]) for (let i = 0; i < pts.length; i++) {
    const [x, y] = pts[i], q = pts[i ? i - 1 : i + 1];
    if (q && Math.abs(x * x + y * y - 1) < 1e-6 && Math.hypot(-x - q[0], y - q[1]) < Math.hypot(x - q[0], y - q[1])) pts[i] = [-x, y];
  }
  CURVES = { nodal: [nodal, mirror(nodal)], open: [open, mirror(open)] };
  return CURVES;
}
