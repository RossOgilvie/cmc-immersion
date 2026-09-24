// Complex polynomials with coefficients stored interleaved in a Float64Array:
// p = [re0, im0, re1, im1, ...] represents sum_k (re_k + i im_k) lam^k.
// Complex scalars in the public API are [re, im] pairs.

export function polyDeg(p) {
  return p.length / 2 - 1;
}

export function polyMul(p, q) {
  const n = p.length / 2, m = q.length / 2;
  const r = new Float64Array(2 * (n + m - 1));
  for (let i = 0; i < n; i++) {
    const pr = p[2 * i], pi = p[2 * i + 1];
    for (let j = 0; j < m; j++) {
      const qr = q[2 * j], qi = q[2 * j + 1];
      r[2 * (i + j)] += pr * qr - pi * qi;
      r[2 * (i + j) + 1] += pr * qi + pi * qr;
    }
  }
  return r;
}

export function polyEval(p, zr, zi) {
  let re = 0, im = 0;
  for (let k = p.length / 2 - 1; k >= 0; k--) {
    const t = re * zr - im * zi + p[2 * k];
    im = re * zi + im * zr + p[2 * k + 1];
    re = t;
  }
  return [re, im];
}

/** p^sharp(lam) = lam^deg * conj(p(1/conj lam)): reverse and conjugate the coefficients. */
export function polySharp(p, deg) {
  const r = new Float64Array(2 * (deg + 1));
  const n = p.length / 2;
  for (let k = 0; k <= deg; k++) {
    const j = deg - k;
    if (j < n) { r[2 * k] = p[2 * j]; r[2 * k + 1] = -p[2 * j + 1]; }
  }
  return r;
}

/** All complex roots by Aberth–Ehrlich iteration. Returns an array of [re, im]. */
export function polyRoots(p) {
  let n = p.length / 2 - 1;
  while (n > 0 && p[2 * n] === 0 && p[2 * n + 1] === 0) n--;
  if (n <= 0) return [];
  // monic copy
  const lr = p[2 * n], li = p[2 * n + 1], ld = lr * lr + li * li;
  const c = new Float64Array(2 * (n + 1));
  for (let k = 0; k <= n; k++) {
    const a = p[2 * k], b = p[2 * k + 1];
    c[2 * k] = (a * lr + b * li) / ld;
    c[2 * k + 1] = (b * lr - a * li) / ld;
  }
  // derivative
  const d = new Float64Array(2 * n);
  for (let k = 1; k <= n; k++) { d[2 * (k - 1)] = k * c[2 * k]; d[2 * (k - 1) + 1] = k * c[2 * k + 1]; }
  // Cauchy bound for the initial circle
  let R = 0;
  for (let k = 0; k < n; k++) R = Math.max(R, Math.hypot(c[2 * k], c[2 * k + 1]));
  R = Math.min(1 + R, 1e6);
  const zr = new Float64Array(n), zi = new Float64Array(n);
  for (let k = 0; k < n; k++) {
    const t = (2 * Math.PI * k) / n + 0.4;
    zr[k] = 0.5 * R * Math.cos(t);
    zi[k] = 0.5 * R * Math.sin(t);
  }
  for (let iter = 0; iter < 500; iter++) {
    let maxStep = 0;
    for (let k = 0; k < n; k++) {
      const [pr, pi] = polyEval(c, zr[k], zi[k]);
      const [dr, di] = polyEval(d, zr[k], zi[k]);
      // ratio = p / p'
      const dd = dr * dr + di * di;
      if (dd === 0) continue;
      const rr = (pr * dr + pi * di) / dd, ri = (pi * dr - pr * di) / dd;
      // sum 1/(z_k - z_j)
      let sr = 0, si = 0;
      for (let j = 0; j < n; j++) {
        if (j === k) continue;
        const ar = zr[k] - zr[j], ai = zi[k] - zi[j], a2 = ar * ar + ai * ai;
        sr += ar / a2; si -= ai / a2;
      }
      // w = ratio / (1 - ratio * s)
      const denr = 1 - (rr * sr - ri * si), deni = -(rr * si + ri * sr);
      const den2 = denr * denr + deni * deni;
      const wr = (rr * denr + ri * deni) / den2, wi = (ri * denr - rr * deni) / den2;
      zr[k] -= wr; zi[k] -= wi;
      maxStep = Math.max(maxStep, Math.hypot(wr, wi) / (1 + Math.hypot(zr[k], zi[k])));
    }
    if (maxStep < 1e-15) break;
  }
  const out = [];
  for (let k = 0; k < n; k++) out.push([zr[k], zi[k]]);
  return out;
}
