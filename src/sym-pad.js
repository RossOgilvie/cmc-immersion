// The pad for the root-preserving flow: a square window on the plane of Sym integrals (phi_1, phi_2).
//
//  - the current point: a rust dot, draggable; dragging sets a target that the flow chases
//  - where the flow couldn't follow a target, a grey ring marks the target
//
// Two modes:
//  - 'triangle' (genus 2): the window is [0, pi]^2. Faint lines at rational multiples of pi with small
//    denominators; where two cross, the data close up into a torus, and a dragged target snaps onto such
//    a crossing (hold alt to avoid). The diagonal phi_1 = phi_2 (dashed) is the Wente family; the flow
//    fills the triangle phi_1 + phi_2 < pi (CKKS), so the far half of the square is shaded.
//  - 'free' (genus >= 4): the image of phi is not known. The window is centred on the anchor and zooms
//    with the wheel (and widens by itself when the point nears its edge); a faint trail shows where the
//    flow has been and small crosses where it stopped, which maps out the reachable region.

const Q_MAX = 6; // denominators of the grid lines and snap points
const PI = Math.PI;

export class SymPad {
  /** onTarget(phi, dragging) when the user drags */
  constructor(canvas, onTarget) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.onTarget = onTarget;
    this.phi = null;       // current point, or null (disabled)
    this.target = null;    // point being chased
    this.blocked = false;  // the flow stopped short of the target
    this.drag = false;
    this.mode = 'triangle';
    this.view = [0, PI, 0, PI]; // [x0, x1, y0, y1] in phi
    this.trail = [];       // free mode: points the flow has passed
    this.marks = [];       // free mode: points where it stopped
    this._bind();
    new ResizeObserver(() => this.draw()).observe(canvas);
    if (document.fonts) document.fonts.load("14px 'EB Garamond'", 'φπ01').then(() => this.draw(), () => {});
  }

  set(opts) {
    Object.assign(this, opts);
    if (this.mode === 'triangle') this.view = [0, PI, 0, PI];
    else if (this.phi) this._fit(this.phi);
    this.draw();
  }

  /** Free mode: a window of half-width R centred on phi. */
  centre(phi, R = 0.02 * PI) {
    this.view = [phi[0] - R, phi[0] + R, phi[1] - R, phi[1] + R];
  }

  // widen the window (about its centre) until p is well inside
  _fit(p) {
    for (let k = 0; k < 20; k++) {
      const [x0, x1, y0, y1] = this.view, mx = 0.08 * (x1 - x0), my = 0.08 * (y1 - y0);
      if (p[0] > x0 + mx && p[0] < x1 - mx && p[1] > y0 + my && p[1] < y1 - my) return;
      const cx = (x0 + x1) / 2, cy = (y0 + y1) / 2, R = x1 - x0;
      this.view = [cx - R, cx + R, cy - R, cy + R];
    }
  }

  _geom() {
    const w = this.canvas.clientWidth, h = this.canvas.clientHeight;
    const pad = this.mode === 'triangle' ? { l: 24, r: 10, t: 10, b: 22 } : { l: 24, r: 10, t: 10, b: 30 };
    const s = Math.min(w - pad.l - pad.r, h - pad.t - pad.b);
    return { w, h, s, x0: pad.l, y0: pad.t + s };
  }
  _toPx([p1, p2]) {
    const { s, x0, y0 } = this._geom(), [a, b, c, d] = this.view;
    return [x0 + ((p1 - a) / (b - a)) * s, y0 - ((p2 - c) / (d - c)) * s];
  }
  _fromPx(x, y) {
    const { s, x0, y0 } = this._geom(), [a, b, c, d] = this.view;
    const cl = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
    return [cl(a + ((x - x0) / s) * (b - a), a, b), cl(c + ((y0 - y) / s) * (d - c), c, d)];
  }

  _snap(phi) {
    const q = this._toPx(phi);
    let best = null, bd = 7;
    for (const u of rationals()) for (const v of rationals()) {
      const p = this._toPx([u * PI, v * PI]);
      const d = Math.hypot(p[0] - q[0], p[1] - q[1]);
      if (d < bd) { bd = d; best = [u * PI, v * PI]; }
    }
    return best || phi;
  }

  _bind() {
    const c = this.canvas;
    const pos = (e) => { const r = c.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; };
    const emit = (e, dragging) => {
      let phi = this._fromPx(...pos(e));
      if (this.mode === 'triangle' && !e.altKey) phi = this._snap(phi);
      this.target = phi;
      this.blocked = false;
      this.draw();
      this.onTarget(phi, dragging);
    };
    c.addEventListener('pointerdown', (e) => {
      if (!this.phi) return;
      this.drag = true;
      c.setPointerCapture(e.pointerId);
      e.preventDefault();
      emit(e, true);
    });
    c.addEventListener('pointermove', (e) => {
      if (!this.phi) { c.style.cursor = 'default'; return; }
      c.style.cursor = this.drag ? 'grabbing' : 'crosshair';
      if (this.drag) emit(e, true);
    });
    const up = (e) => {
      if (!this.drag) return;
      this.drag = false;
      emit(e, false);
    };
    c.addEventListener('pointerup', up);
    c.addEventListener('pointercancel', up);
    // free mode: the wheel zooms about the current point
    c.addEventListener('wheel', (e) => {
      if (this.mode !== 'free' || !this.phi) return;
      e.preventDefault();
      const f = Math.exp((e.deltaMode === 1 ? e.deltaY * 16 : e.deltaY) * 0.0015);
      const [a, b, cc, d] = this.view, [px, py] = this.phi;
      this.view = [px + (a - px) * f, px + (b - px) * f, py + (cc - py) * f, py + (d - py) * f];
      this.draw();
    }, { passive: false });
  }

  draw() {
    const c = this.canvas, ctx = this.ctx;
    const dpr = Math.min(window.devicePixelRatio || 1, 3);
    const W = Math.round(c.clientWidth * dpr), H = Math.round(c.clientHeight * dpr);
    if (!W || !H) return;
    if (c.width !== W || c.height !== H) { c.width = W; c.height = H; }
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const css = getComputedStyle(c);
    const col = (n, f) => css.getPropertyValue(n).trim() || f;
    const C = {
      ink: col('--ink', '#27231E'), muted: col('--muted', '#6A6257'), rule: col('--rule', '#CFC5B0'),
      accent: col('--accent', '#9E3F24'), shade: col('--shade', '#EAE4D6'),
    };
    const { w, h, s, x0, y0 } = this._geom();
    ctx.clearRect(0, 0, w, h);
    const on = !!this.phi;
    ctx.globalAlpha = on ? 1 : 0.45;
    const FONT = "'EB Garamond', Georgia, serif";
    if (this.mode === 'triangle') this._drawTriangle(ctx, C, FONT, on);
    else this._drawFree(ctx, C, FONT, on);
    ctx.globalAlpha = 1;
    if (!on) return;

    // target the flow couldn't reach, or is still chasing
    if (this.target && Math.hypot(this.target[0] - this.phi[0], this.target[1] - this.phi[1]) > 1e-9) {
      const [x, y] = this._toPx(this.target);
      ctx.strokeStyle = this.blocked ? C.muted : C.accent;
      ctx.lineWidth = 1.3;
      ctx.beginPath(); ctx.arc(x, y, 6, 0, 2 * Math.PI); ctx.stroke();
    }
    const [x, y] = this._toPx(this.phi);
    ctx.save();
    ctx.beginPath(); ctx.rect(x0 - 6, y0 - s - 6, s + 12, s + 12); ctx.clip();
    ctx.fillStyle = C.accent;
    ctx.beginPath(); ctx.arc(x, y, this.drag ? 6 : 5, 0, 2 * Math.PI); ctx.fill();
    ctx.restore();
  }

  _drawTriangle(ctx, C, FONT, on) {
    const { s, x0, y0 } = this._geom();
    // the half of the square beyond the triangle phi_1 + phi_2 < pi
    ctx.fillStyle = C.shade;
    ctx.beginPath(); ctx.moveTo(x0 + s, y0); ctx.lineTo(x0 + s, y0 - s); ctx.lineTo(x0, y0 - s); ctx.closePath(); ctx.fill();

    // rational grid: fainter for larger denominators
    for (const u of rationals()) {
      if (u === 0 || u === 1) continue;
      const q = denominator(u);
      ctx.strokeStyle = C.rule;
      ctx.lineWidth = q <= 2 ? 1 : 0.6;
      ctx.globalAlpha = (on ? 1 : 0.45) * (q <= 3 ? 0.9 : 0.5);
      const x = x0 + u * s, y = y0 - u * s;
      ctx.beginPath(); ctx.moveTo(x, y0); ctx.lineTo(x, y0 - s); ctx.moveTo(x0, y); ctx.lineTo(x0 + s, y); ctx.stroke();
    }
    ctx.globalAlpha = on ? 1 : 0.45;

    // frame, hypotenuse, Wente diagonal
    ctx.strokeStyle = C.ink;
    ctx.lineWidth = 1.2;
    ctx.strokeRect(x0, y0 - s, s, s);
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(x0 + s, y0); ctx.lineTo(x0, y0 - s); ctx.stroke();
    ctx.strokeStyle = C.muted;
    ctx.setLineDash([4, 4]);
    ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x0 + s / 2, y0 - s / 2); ctx.stroke();
    ctx.setLineDash([]);

    // labels
    ctx.fillStyle = C.muted;
    ctx.font = `13px ${FONT}`;
    ctx.textBaseline = 'top'; ctx.textAlign = 'center';
    ctx.fillText('0', x0, y0 + 4); ctx.fillText('π', x0 + s, y0 + 4);
    subLabel(ctx, FONT, '1', x0 + s / 2, y0 + 10, 'center');
    ctx.textBaseline = 'middle'; ctx.textAlign = 'right';
    ctx.fillText('π', x0 - 6, y0 - s);
    subLabel(ctx, FONT, '2', x0 - 6, y0 - s / 2, 'right');
    ctx.save();
    ctx.font = `italic 12px ${FONT}`;
    ctx.translate(x0 + s * 0.3, y0 - s * 0.3);
    ctx.rotate(-Math.PI / 4);
    ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
    ctx.fillText('Wente', 0, -3);
    ctx.restore();
  }

  _drawFree(ctx, C, FONT, on) {
    const { s, x0, y0 } = this._geom(), [a, b, c, d] = this.view;
    // the axes phi_l = 0 and a light grid at a round step (in units of pi)
    const step = niceStep((b - a) / PI / 4) * PI;
    ctx.strokeStyle = C.rule;
    ctx.lineWidth = 0.6;
    ctx.beginPath();
    for (let v = Math.ceil(a / step) * step; v < b; v += step) { const [x] = this._toPx([v, c]); ctx.moveTo(x, y0); ctx.lineTo(x, y0 - s); }
    for (let v = Math.ceil(c / step) * step; v < d; v += step) { const [, y] = this._toPx([a, v]); ctx.moveTo(x0, y); ctx.lineTo(x0 + s, y); }
    ctx.stroke();
    ctx.strokeStyle = C.muted;
    ctx.lineWidth = 0.8;
    ctx.beginPath();
    if (a < 0 && b > 0) { const [x] = this._toPx([0, c]); ctx.moveTo(x, y0); ctx.lineTo(x, y0 - s); }
    if (c < 0 && d > 0) { const [, y] = this._toPx([a, 0]); ctx.moveTo(x0, y); ctx.lineTo(x0 + s, y); }
    ctx.stroke();
    ctx.strokeStyle = C.ink;
    ctx.lineWidth = 1.2;
    ctx.strokeRect(x0, y0 - s, s, s);

    ctx.save();
    ctx.beginPath(); ctx.rect(x0, y0 - s, s, s); ctx.clip();
    // where the flow has been, and where it stopped
    ctx.fillStyle = C.muted;
    ctx.globalAlpha = (on ? 1 : 0.45) * 0.35;
    for (const p of this.trail) { const [x, y] = this._toPx(p); ctx.fillRect(x - 1, y - 1, 2, 2); }
    ctx.globalAlpha = on ? 1 : 0.45;
    ctx.strokeStyle = C.ink;
    ctx.lineWidth = 1.1;
    for (const p of this.marks) {
      const [x, y] = this._toPx(p);
      ctx.beginPath(); ctx.moveTo(x - 3, y - 3); ctx.lineTo(x + 3, y + 3); ctx.moveTo(x + 3, y - 3); ctx.lineTo(x - 3, y + 3); ctx.stroke();
    }
    ctx.restore();

    // the window's extent, in units of pi
    const dig = Math.max(2, Math.min(5, 1 - Math.floor(Math.log10(step / PI))));
    const f = (v) => { const t = (Math.abs(v / PI) < 10 ** -dig ? 0 : v / PI).toFixed(dig); return (t.startsWith('-') ? '−' + t.slice(1) : t) + 'π'; };
    ctx.fillStyle = C.muted;
    ctx.font = `12px ${FONT}`;
    ctx.textBaseline = 'top';
    ctx.textAlign = 'left'; ctx.fillText(f(a), x0, y0 + 4);
    ctx.textAlign = 'right'; ctx.fillText(f(b), x0 + s, y0 + 4);
    subLabel(ctx, FONT, '1', x0 + s / 2, y0 + 12, 'center');
    ctx.save();
    ctx.translate(x0 - 5, y0);
    ctx.rotate(-Math.PI / 2);
    ctx.textBaseline = 'bottom';
    ctx.textAlign = 'left'; ctx.fillText(f(c), 0, 0);
    ctx.textAlign = 'right'; ctx.fillText(f(d), s, 0);
    ctx.restore();
    subLabel(ctx, FONT, '2', x0 - 6, y0 - s / 2, 'right');
  }
}

/** 1, 2 or 5 times a power of ten, at least x. */
function niceStep(x) {
  const p = 10 ** Math.floor(Math.log10(x)), m = x / p;
  return (m <= 1 ? 1 : m <= 2 ? 2 : m <= 5 ? 5 : 10) * p;
}

/** phi with a subscript, the pair aligned at (x, y) (vertically centred), anchored left/center/right. */
function subLabel(ctx, FONT, sub, x, y, align) {
  ctx.font = `italic 14px ${FONT}`;
  const wb = ctx.measureText('φ').width;
  ctx.font = `10px ${FONT}`;
  const ws = ctx.measureText(sub).width;
  const left = align === 'center' ? x - (wb + ws) / 2 : align === 'right' ? x - wb - ws : x;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'middle';
  ctx.font = `italic 14px ${FONT}`;
  ctx.fillText('φ', left, y);
  ctx.font = `10px ${FONT}`;
  ctx.fillText(sub, left + wb, y + 4);
}

/** The rationals in [0, 1] with denominator <= Q_MAX. */
let RAT = null;
function rationals() {
  if (RAT) return RAT;
  const set = new Map();
  for (let q = 1; q <= Q_MAX; q++) for (let p = 0; p <= q; p++) set.set((p / q).toFixed(12), p / q);
  RAT = [...set.values()].sort((a, b) => a - b);
  return RAT;
}
export function denominator(u, qmax = Q_MAX, tol = 1e-9) {
  for (let q = 1; q <= qmax; q++) if (Math.abs(u * q - Math.round(u * q)) < tol) return q;
  return Infinity;
}
