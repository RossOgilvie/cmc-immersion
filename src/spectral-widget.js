// Interactive editor for the spectral data in the lambda-plane.
//
//  - branch points alpha_i: filled dots, draggable within the punctured unit disc
//  - Sym point lam0 = e^{i theta0}: a rust dot, draggable along the unit circle; it snaps onto
//  - common roots on S^1 of the differentials Theta_w: small rust diamonds
//  - reflections 1/conj(alpha_i): faint rings (read-only)
//  - divisor points mu_j: hollow grey rings (read-only), clamped to the edge with an arrow if far out
//
// Shift+drag rotates all alpha_i and lam0 together (a rotation of the z-plane: same surface, turned).
// In genus 1 the diameter through lam0 is drawn: there the surface closes into a Delaunay cylinder
// (unduloid for alpha on the far side, nodoid on the near side), and alpha snaps onto it.
// Double-click on empty space adds a branch point; double-click on a branch point removes it.
// If given a Whitham family, draws the path of each branch point along it, coloured by the Willmore
// functional W (blue low, rust high), with a tick across the path at the critical points of W.
// Labels sit on the side of each point away from the origin and from its own trail.

const R_VIEW = 1.25;       // half-width of the view in lambda units
const R_MIN = 0.03, R_MAX = 0.97;
export const MAX_GENUS = 8;

export class SpectralWidget {
  /** onChange(state, dragging) with state = { alphas, theta0 } */
  constructor(canvas, onChange) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.onChange = onChange;
    this.alphas = [];
    this.theta0 = 0;
    this.divisor = [];
    this.commonRoots = [];
    this.drag = null;
    this.hover = null;
    this.colours = {};
    this._bind();
    new ResizeObserver(() => this.draw()).observe(canvas);
    // canvas text doesn't trigger web-font loading, and draws in the fallback until the font is in
    if (document.fonts) document.fonts.load("16px 'EB Garamond'", 'αλ01').then(() => this.draw(), () => {});
  }

  /** family: { points: [{ s, alphas, W }], critical: [indices], ends: { neg, pos } } or null */
  set(opts) {
    const { alphas, theta0, divisor } = opts;
    if (alphas) this.alphas = alphas.map((a) => a.slice());
    if (theta0 !== undefined) this.theta0 = theta0;
    if (divisor) this.divisor = divisor;
    if ('family' in opts) this.family = opts.family;
    if ('commonRoots' in opts) this.commonRoots = opts.commonRoots || [];
    this.draw();
  }

  _emit(dragging) {
    this.onChange({ alphas: this.alphas.map((a) => a.slice()), theta0: this.theta0 }, dragging);
  }

  // ---------------------------------------------------------------- coordinates

  _geom() {
    const w = this.canvas.clientWidth, h = this.canvas.clientHeight;
    const s = Math.min(w, h) / (2 * R_VIEW);
    return { w, h, s, cx: w / 2, cy: h / 2 };
  }
  _toPx([re, im]) {
    const { s, cx, cy } = this._geom();
    return [cx + re * s, cy - im * s];
  }
  _fromPx(x, y) {
    const { s, cx, cy } = this._geom();
    return [(x - cx) / s, -(y - cy) / s];
  }
  _eventPos(e) {
    const r = this.canvas.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top];
  }

  _pick(x, y) {
    const tol = 13;
    let best = null, bd = tol;
    const lp = this._toPx([Math.cos(this.theta0), Math.sin(this.theta0)]);
    const dl = Math.hypot(lp[0] - x, lp[1] - y);
    if (dl < bd) { best = { kind: 'lam' }; bd = dl; }
    this.alphas.forEach((a, i) => {
      const p = this._toPx(a);
      const d = Math.hypot(p[0] - x, p[1] - y);
      if (d < bd) { best = { kind: 'alpha', i }; bd = d; }
    });
    return best;
  }

  // ---------------------------------------------------------------- interaction

  _bind() {
    const c = this.canvas;
    c.addEventListener('pointerdown', (e) => {
      this.noSnap = e.altKey;
      const [x, y] = this._eventPos(e);
      const hit = this._pick(x, y);
      if (e.shiftKey) {
        const z = this._fromPx(x, y);
        this.drag = { kind: 'rotate', ang: Math.atan2(z[1], z[0]) };
      } else if (hit) {
        this.drag = hit;
      } else {
        return;
      }
      c.setPointerCapture(e.pointerId);
      e.preventDefault();
    });
    c.addEventListener('pointermove', (e) => {
      const [x, y] = this._eventPos(e);
      if (!this.drag) {
        const h = this._pick(x, y);
        const key = h ? `${h.kind}${h.i ?? ''}` : null;
        if (key !== this.hover) { this.hover = key; this.draw(); }
        c.style.cursor = h ? 'grab' : e.shiftKey ? 'alias' : 'default';
        return;
      }
      const z = this._fromPx(x, y);
      const d = this.drag;
      if (d.kind === 'alpha') {
        let r = Math.hypot(z[0], z[1]);
        const t = Math.atan2(z[1], z[0]);
        r = Math.min(R_MAX, Math.max(R_MIN, r));
        this.alphas[d.i] = this._snap([r * Math.cos(t), r * Math.sin(t)]);
      } else if (d.kind === 'lam') {
        this.theta0 = Math.atan2(z[1], z[0]);
        // snap onto a common root of the differentials (hold alt to avoid)
        for (const [re, im] of this.commonRoots) {
          const p = this._toPx([re, im]), q = this._toPx([Math.cos(this.theta0), Math.sin(this.theta0)]);
          if (!this.noSnap && Math.hypot(p[0] - q[0], p[1] - q[1]) < 8) this.theta0 = Math.atan2(im, re);
        }
      } else if (d.kind === 'rotate') {
        const a = Math.atan2(z[1], z[0]);
        const da = a - d.ang;
        d.ang = a;
        const cr = Math.cos(da), ci = Math.sin(da);
        this.alphas = this.alphas.map(([re, im]) => [re * cr - im * ci, re * ci + im * cr]);
        this.theta0 += da;
      }
      c.style.cursor = 'grabbing';
      this.draw();
      this._emit(true);
    });
    const up = () => {
      if (!this.drag) return;
      this.drag = null;
      this.canvas.style.cursor = 'default';
      this._emit(false);
    };
    c.addEventListener('pointerup', up);
    c.addEventListener('pointercancel', up);
    c.addEventListener('dblclick', (e) => {
      const [x, y] = this._eventPos(e);
      const hit = this._pick(x, y);
      if (hit && hit.kind === 'alpha') {
        this.alphas.splice(hit.i, 1);
      } else if (!hit) {
        const z = this._fromPx(x, y);
        const r = Math.hypot(z[0], z[1]);
        if (r > 1 || this.alphas.length >= MAX_GENUS) return;
        const rr = Math.min(R_MAX, Math.max(R_MIN, r));
        const t = Math.atan2(z[1], z[0]);
        this.alphas.push([rr * Math.cos(t), rr * Math.sin(t)]);
      } else {
        return;
      }
      this.draw();
      this._emit(false);
    });
  }

  /** Genus 1: snap alpha onto the Delaunay diameter when within a few pixels (hold alt to avoid). */
  _snap(z) {
    if (this.alphas.length !== 1 || this.noSnap) return z;
    const ur = Math.cos(this.theta0), ui = Math.sin(this.theta0);
    const along = z[0] * ur + z[1] * ui;
    const off = -z[0] * ui + z[1] * ur;
    if (Math.abs(off) * this._geom().s > 7) return z;
    const r = Math.sign(along || 1) * Math.min(R_MAX, Math.max(R_MIN, Math.abs(along)));
    return [r * ur, r * ui];
  }

  // ---------------------------------------------------------------- drawing

  _drawFamily(ctx, ink) {
    this._trails = null;
    const F = this.family;
    if (!F || !F.points || F.points.length < 2) return;
    const P = F.points, g = P[0].alphas.length;
    if (g !== this.alphas.length) return; // stale (genus changed)
    const Ws = P.map((p) => p.W);
    const lo = Math.min(...Ws), hi = Math.max(...Ws), span = hi - lo || 1;
    const c0 = rgb(this._col('--blue', '#3D5A80')), c1 = rgb(this._col('--accent', '#9E3F24'));
    const ramp = (t) => `rgb(${c0.map((v, i) => Math.round(v + t * (c1[i] - v))).join(',')})`; // blue low, rust high
    this._trails = [];
    ctx.save();
    ctx.lineCap = 'round';
    for (let j = 0; j < g; j++) {
      const px = P.map((p) => this._toPx(p.alphas[j]));
      this._trails.push(px);
      ctx.lineWidth = 2.4;
      ctx.globalAlpha = 0.85;
      for (let i = 1; i < P.length; i++) {
        ctx.strokeStyle = ramp(((P[i - 1].W + P[i].W) / 2 - lo) / span);
        ctx.beginPath(); ctx.moveTo(...px[i - 1]); ctx.lineTo(...px[i]); ctx.stroke();
      }
      // critical points of W: a short tick across the path
      ctx.globalAlpha = 1;
      ctx.lineWidth = 1.6;
      ctx.strokeStyle = ink;
      for (const c of F.critical || []) {
        const a = px[Math.max(0, c - 1)], b = px[Math.min(P.length - 1, c + 1)];
        let dx = b[0] - a[0], dy = b[1] - a[1];
        const d = Math.hypot(dx, dy) || 1; dx /= d; dy /= d;
        const [x, y] = px[c];
        ctx.beginPath(); ctx.moveTo(x - dy * 5, y + dx * 5); ctx.lineTo(x + dy * 5, y - dx * 5); ctx.stroke();
      }
    }
    ctx.restore();
  }

  draw() {
    const c = this.canvas, ctx = this.ctx;
    const dpr = Math.min(window.devicePixelRatio || 1, 3);
    const W = Math.round(c.clientWidth * dpr), H = Math.round(c.clientHeight * dpr);
    if (!W || !H) return;
    if (c.width !== W || c.height !== H) { c.width = W; c.height = H; }
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const { w, h, s, cx, cy } = this._geom();
    const ink = this._col('--ink', '#27231E'), muted = this._col('--muted', '#6A6257');
    const rule = this._col('--rule', '#CFC5B0'), accent = this._col('--accent', '#9E3F24');
    const paper = this._col('--paper', '#F2EEE3'), disc = this._col('--disc', 'rgba(234,228,214,0.55)');
    ctx.clearRect(0, 0, w, h);

    // axes and circles
    ctx.lineWidth = 1;
    ctx.strokeStyle = rule;
    ctx.beginPath();
    ctx.moveTo(0, cy); ctx.lineTo(w, cy); ctx.moveTo(cx, 0); ctx.lineTo(cx, h);
    ctx.stroke();
    ctx.fillStyle = disc;
    ctx.beginPath(); ctx.arc(cx, cy, s, 0, 2 * Math.PI); ctx.fill();
    ctx.setLineDash([2, 4]);
    ctx.beginPath(); ctx.arc(cx, cy, 0.5 * s, 0, 2 * Math.PI); ctx.stroke();
    ctx.setLineDash([]);
    ctx.strokeStyle = ink;
    ctx.lineWidth = 1.4;
    ctx.beginPath(); ctx.arc(cx, cy, s, 0, 2 * Math.PI); ctx.stroke();

    // genus 1: the Delaunay diameter
    if (this.alphas.length === 1) {
      const ur = Math.cos(this.theta0), ui = Math.sin(this.theta0);
      ctx.save();
      ctx.strokeStyle = muted;
      ctx.globalAlpha = 0.7;
      ctx.setLineDash([5, 4]);
      ctx.lineWidth = 1.2;
      ctx.beginPath();
      ctx.moveTo(cx - ur * s, cy + ui * s); ctx.lineTo(cx + ur * s, cy - ui * s);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.globalAlpha = 1;
      ctx.fillStyle = muted;
      ctx.font = `italic 13px ${FONT}`;
      ctx.textAlign = 'center';
      ctx.textBaseline = 'alphabetic';
      const lab = (t, txt) => {
        ctx.save();
        ctx.translate(cx + t * ur * s, cy - t * ui * s);
        let a = -this.theta0;
        if (Math.cos(a) < 0) a += Math.PI;
        ctx.rotate(a);
        ctx.fillText(txt, 0, -5);
        ctx.restore();
      };
      lab(0.62, 'Nodoids');
      lab(-0.62, 'Unduloids');
      ctx.restore();
    }

    // branch point at 0: a small grey cross
    ctx.strokeStyle = muted;
    ctx.lineWidth = 1.1;
    ctx.beginPath();
    ctx.moveTo(cx - 3, cy - 3); ctx.lineTo(cx + 3, cy + 3);
    ctx.moveTo(cx - 3, cy + 3); ctx.lineTo(cx + 3, cy - 3);
    ctx.stroke();

    // reflections 1/conj(alpha)
    ctx.strokeStyle = ink;
    ctx.lineWidth = 1;
    ctx.globalAlpha = 0.3;
    for (const [re, im] of this.alphas) {
      const r2 = re * re + im * im;
      const p = this._toPx([re / r2, im / r2]);
      if (p[0] < -5 || p[0] > w + 5 || p[1] < -5 || p[1] > h + 5) continue;
      ctx.beginPath(); ctx.arc(p[0], p[1], 3.5, 0, 2 * Math.PI); ctx.stroke();
    }
    ctx.globalAlpha = 1;

    // divisor
    ctx.strokeStyle = muted;
    ctx.fillStyle = muted;
    ctx.lineWidth = 1.3;
    for (const { mu } of this.divisor) {
      let p = this._toPx(mu);
      const m = 8;
      const out = p[0] < m || p[0] > w - m || p[1] < m || p[1] > h - m;
      if (out) {
        // clamp to the frame along the ray from the centre, and draw an arrowhead
        const dx = p[0] - cx, dy = p[1] - cy;
        const k = Math.min((w / 2 - m) / Math.abs(dx || 1e-9), (h / 2 - m) / Math.abs(dy || 1e-9));
        p = [cx + dx * k, cy + dy * k];
        const a = Math.atan2(dy, dx);
        ctx.beginPath();
        ctx.moveTo(p[0] + 6 * Math.cos(a), p[1] + 6 * Math.sin(a));
        ctx.lineTo(p[0] + 6 * Math.cos(a + 2.4), p[1] + 6 * Math.sin(a + 2.4));
        ctx.lineTo(p[0] + 6 * Math.cos(a - 2.4), p[1] + 6 * Math.sin(a - 2.4));
        ctx.closePath(); ctx.fill();
      } else {
        ctx.beginPath(); ctx.arc(p[0], p[1], 6, 0, 2 * Math.PI); ctx.stroke();
      }
    }

    this._drawFamily(ctx, ink);

    // branch points
    const lp = this._toPx([Math.cos(this.theta0), Math.sin(this.theta0)]);
    const pts = this.alphas.map((a) => this._toPx(a));
    pts.forEach((p, i) => {
      const hot = this.hover === `alpha${i}` || (this.drag && this.drag.i === i);
      ctx.fillStyle = ink;
      ctx.beginPath(); ctx.arc(p[0], p[1], hot ? 6 : 5, 0, 2 * Math.PI); ctx.fill();
    });
    pts.forEach((p, i) => {
      const others = [...pts.filter((_, j) => j !== i), lp];
      this._label(ctx, 'α', String(i + 1), p, this._labelDir(p, i, others), ink, 9);
    });

    // Sym point, labelled outside the circle
    const hot = this.hover === 'lam' || (this.drag && this.drag.kind === 'lam');
    ctx.fillStyle = accent;
    ctx.beginPath(); ctx.arc(lp[0], lp[1], hot ? 6.5 : 5, 0, 2 * Math.PI); ctx.fill();
    this._label(ctx, 'λ', '0', lp, [Math.cos(this.theta0), -Math.sin(this.theta0)], accent, 9);

    // common roots of the differentials on S^1: small diamonds (outlined in paper, so one sitting
    // under the Sym point still shows)
    for (const [re, im] of this.commonRoots) {
      const p = this._toPx([re, im]), r = 3.8;
      ctx.beginPath();
      ctx.moveTo(p[0], p[1] - r); ctx.lineTo(p[0] + r, p[1]);
      ctx.lineTo(p[0], p[1] + r); ctx.lineTo(p[0] - r, p[1]);
      ctx.closePath();
      ctx.fillStyle = accent; ctx.fill();
      ctx.lineWidth = 1.2; ctx.strokeStyle = paper; ctx.stroke();
    }
  }

  _col(name, fallback) {
    return getComputedStyle(this.canvas).getPropertyValue(name).trim() || fallback;
  }

  /**
   * Unit vector (pixel coordinates) for the label of the branch point at p: away from the origin, and
   * clear of its own trail and of nearby points.
   */
  _labelDir(p, i, others) {
    const { cx, cy } = this._geom();
    const avoid = []; // [ux, uy, weight]
    const trail = this._trails && this._trails[i];
    if (trail) {
      // the trail leaves the point in up to two directions: look ~14px along it on either side
      let k = 0, bd = Infinity;
      trail.forEach((q, m) => { const d = Math.hypot(q[0] - p[0], q[1] - p[1]); if (d < bd) { bd = d; k = m; } });
      for (const dir of [-1, 1]) {
        for (let m = k + dir; m >= 0 && m < trail.length; m += dir) {
          const dx = trail[m][0] - p[0], dy = trail[m][1] - p[1], d = Math.hypot(dx, dy);
          if (d >= 14) { avoid.push([dx / d, dy / d, 1.6]); break; }
        }
      }
    }
    for (const q of others) {
      const dx = q[0] - p[0], dy = q[1] - p[1], d = Math.hypot(dx, dy);
      if (d > 0 && d < 40) avoid.push([dx / d, dy / d, 0.9]);
    }
    const rx = p[0] - cx, ry = p[1] - cy, rl = Math.hypot(rx, ry);
    const radial = rl > 1 ? Math.atan2(ry, rx) : -Math.PI / 4;
    let best = [1, 0], bs = -Infinity;
    for (let k = 0; k < 24; k++) {
      const t = (k * Math.PI) / 12, u = [Math.cos(t), Math.sin(t)];
      let sc = Math.cos(t - radial);
      for (const [ax, ay, wt] of avoid) sc -= wt * Math.max(0, u[0] * ax + u[1] * ay);
      if (sc > bs) { bs = sc; best = u; }
    }
    return best;
  }

  /** Draws base with a subscript, centred gap px beyond p in the direction u. */
  _label(ctx, base, subText, p, u, colour, gap) {
    ctx.font = `16px ${FONT}`;
    const wb = ctx.measureText(base).width;
    ctx.font = `11px ${FONT}`;
    const ws = ctx.measureText(subText).width;
    const wd = wb + ws, ht = 13;
    const d = gap + Math.abs(u[0]) * wd / 2 + Math.abs(u[1]) * ht / 2;
    const x = p[0] + u[0] * d - wd / 2, y = p[1] + u[1] * d;
    ctx.fillStyle = colour;
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    ctx.font = `16px ${FONT}`;
    ctx.fillText(base, x, y);
    ctx.font = `11px ${FONT}`;
    ctx.fillText(subText, x + wb, y + 4);
  }
}

const FONT = "'EB Garamond', Georgia, serif";

/** '#rrggbb' -> [r, g, b] */
function rgb(hex) {
  const m = /^#?([0-9a-f]{6})$/i.exec(hex);
  const n = m ? parseInt(m[1], 16) : 0;
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
}
