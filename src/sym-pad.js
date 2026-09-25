// The pad for the root-preserving flow: the square [0, pi]^2 of Sym integrals (phi_1, phi_2).
//
//  - the current point: a rust dot, draggable; dragging sets a target that the flow chases
//  - faint lines at rational multiples of pi with small denominators; where two cross, the data close
//    up into a torus, and a dragged target snaps onto such a crossing (hold alt to avoid)
//  - the diagonal phi_1 = phi_2 (dashed) is the Wente family; for genus 2 the flow fills the triangle
//    phi_1 + phi_2 < pi (CKKS), so the far half of the square is shaded
//  - where the flow couldn't follow a target (an end of S^2), a grey ring marks the target

const Q_MAX = 6; // denominators of the grid lines and snap points

export class SymPad {
  /** onTarget(phi, dragging) when the user drags; phi in [0, pi]^2 */
  constructor(canvas, onTarget) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.onTarget = onTarget;
    this.phi = null;       // current point, or null (disabled)
    this.target = null;    // point being chased
    this.blocked = false;  // the flow stopped short of the target
    this.drag = false;
    this._bind();
    new ResizeObserver(() => this.draw()).observe(canvas);
    if (document.fonts) document.fonts.load("14px 'EB Garamond'", 'φπ01').then(() => this.draw(), () => {});
  }

  set(opts) {
    Object.assign(this, opts);
    this.draw();
  }

  _geom() {
    const w = this.canvas.clientWidth, h = this.canvas.clientHeight;
    const pad = { l: 24, r: 10, t: 10, b: 22 };
    const s = Math.min(w - pad.l - pad.r, h - pad.t - pad.b);
    return { w, h, s, x0: pad.l, y0: pad.t + s };
  }
  _toPx([p1, p2]) {
    const { s, x0, y0 } = this._geom();
    return [x0 + (p1 / Math.PI) * s, y0 - (p2 / Math.PI) * s];
  }
  _fromPx(x, y) {
    const { s, x0, y0 } = this._geom();
    const c = (v) => Math.min(Math.PI, Math.max(0, v));
    return [c(((x - x0) / s) * Math.PI), c(((y0 - y) / s) * Math.PI)];
  }

  _snap(phi) {
    const q = this._toPx(phi);
    let best = null, bd = 7;
    for (const u of rationals()) for (const v of rationals()) {
      const p = this._toPx([u * Math.PI, v * Math.PI]);
      const d = Math.hypot(p[0] - q[0], p[1] - q[1]);
      if (d < bd) { bd = d; best = [u * Math.PI, v * Math.PI]; }
    }
    return best || phi;
  }

  _bind() {
    const c = this.canvas;
    const pos = (e) => { const r = c.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; };
    const emit = (e, dragging) => {
      let phi = this._fromPx(...pos(e));
      if (!e.altKey) phi = this._snap(phi);
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
    const ink = col('--ink', '#27231E'), muted = col('--muted', '#6A6257'), rule = col('--rule', '#CFC5B0');
    const accent = col('--accent', '#9E3F24'), shade = col('--shade', '#EAE4D6');
    const { w, h, s, x0, y0 } = this._geom();
    ctx.clearRect(0, 0, w, h);
    const on = !!this.phi;
    ctx.globalAlpha = on ? 1 : 0.45;

    // the half of the square beyond the triangle phi_1 + phi_2 < pi
    ctx.fillStyle = shade;
    ctx.beginPath(); ctx.moveTo(x0 + s, y0); ctx.lineTo(x0 + s, y0 - s); ctx.lineTo(x0, y0 - s); ctx.closePath(); ctx.fill();

    // rational grid: fainter for larger denominators
    for (const u of rationals()) {
      if (u === 0 || u === 1) continue;
      const q = denominator(u);
      ctx.strokeStyle = rule;
      ctx.lineWidth = q <= 2 ? 1 : 0.6;
      ctx.globalAlpha = (on ? 1 : 0.45) * (q <= 3 ? 0.9 : 0.5);
      const x = x0 + u * s, y = y0 - u * s;
      ctx.beginPath(); ctx.moveTo(x, y0); ctx.lineTo(x, y0 - s); ctx.moveTo(x0, y); ctx.lineTo(x0 + s, y); ctx.stroke();
    }
    ctx.globalAlpha = on ? 1 : 0.45;

    // frame, hypotenuse, Wente diagonal
    ctx.strokeStyle = ink;
    ctx.lineWidth = 1.2;
    ctx.strokeRect(x0, y0 - s, s, s);
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(x0 + s, y0); ctx.lineTo(x0, y0 - s); ctx.stroke();
    ctx.strokeStyle = muted;
    ctx.setLineDash([4, 4]);
    ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x0 + s / 2, y0 - s / 2); ctx.stroke();
    ctx.setLineDash([]);

    // labels
    const FONT = "'EB Garamond', Georgia, serif";
    ctx.fillStyle = muted;
    ctx.font = `13px ${FONT}`;
    ctx.textBaseline = 'top'; ctx.textAlign = 'center';
    ctx.fillText('0', x0, y0 + 4); ctx.fillText('π', x0 + s, y0 + 4);
    subLabel(ctx, FONT, '1', x0 + s / 2, y0 + 10, 'center');
    ctx.textBaseline = 'middle'; ctx.textAlign = 'right';
    ctx.fillText('π', x0 - 6, y0 - s);
    subLabel(ctx, FONT, '2', x0 - 6, y0 - s / 2, 'right');
    ctx.font = `13px ${FONT}`;
    ctx.save();
    ctx.font = `italic 12px ${FONT}`;
    ctx.translate(x0 + s * 0.3, y0 - s * 0.3);
    ctx.rotate(-Math.PI / 4);
    ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
    ctx.fillText('Wente', 0, -3);
    ctx.restore();
    ctx.globalAlpha = 1;
    if (!on) return;

    // target the flow couldn't reach, or is still chasing
    if (this.target && Math.hypot(this.target[0] - this.phi[0], this.target[1] - this.phi[1]) > 1e-6) {
      const [x, y] = this._toPx(this.target);
      ctx.strokeStyle = this.blocked ? muted : accent;
      ctx.lineWidth = 1.3;
      ctx.beginPath(); ctx.arc(x, y, 6, 0, 2 * Math.PI); ctx.stroke();
    }
    const [x, y] = this._toPx(this.phi);
    ctx.fillStyle = accent;
    ctx.beginPath(); ctx.arc(x, y, this.drag ? 6 : 5, 0, 2 * Math.PI); ctx.fill();
  }
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
