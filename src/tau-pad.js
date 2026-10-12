// The τ pad: a square window on the upper half plane, for moving τ along a leaf with ψ = Σ logit|α_j| fixed (leaf.js).
//
// Always drawn: the fundamental domain F_6 of Gamma(2) (|x| <= 1, outside the circles |τ ± 1/2| = 1/2),
// cut into its six copies of F_1 by the unit circle, the circles |τ ± 1| = 1 and the lines x = ±1/2.
// Each copy is labelled by the order of the classes [w_1] = 1, [w_2] = 2, [w_1 + w_2] = 3 of
// Lambda / 2 Lambda by the length of their shortest vectors (classOrder in leaf.js).
//
// Modes: 'relative' (no lattice: τ is relative to the frame at the anchor and F_6 is only for reference),
// or 'F6', 'F3', 'F1' (a lattice basis: the parts of the plane outside the paper's domain for the
// invariant are shaded). As in the free mode of the root-flow pad: the window zooms with the wheel, pans
// with shift-drag and widens by itself when the point nears its edge; a faint trail shows where the
// flow has been and small crosses where it stopped. Double-click shows all of F_6.
//
// Genus 2 (g2: true): the domain is F_3 = { |x| <= 1, |τ| >= 1 }, drawn alone, with its three copies of F_1
// (cut by the lines x = ±1/2 and the arcs |τ ∓ 1| = 1), and the view starts on all of it. curves: the
// boundary of S^2 in F_3 (genus2.js), at the nodal end (solid) and at the open end (dashed).

import { classOrder } from './cmc/leaf.js';

const F6_VIEW = [-1.25, 1.25, -0.1, 2.4];
const LABELS = [[0, 1.6], [0.78, 1.6], [-0.78, 1.6], [0, 0.75], [0.3, 0.62], [-0.3, 0.62], [0.75, 0.82], [-0.75, 0.82], [0.8, 0.45], [-0.8, 0.45]]
  .map((p) => ({ p, text: classOrder(p) }));
const F3_LABELS = LABELS.filter(({ p: [x, y] }) => x * x + y * y > 1);

export class TauPad {
  /** onTarget(tau, dragging) when the user drags */
  constructor(canvas, onTarget) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.onTarget = onTarget;
    this.tau = null;       // current point, or null (disabled)
    this.target = null;
    this.blocked = false;
    this.drag = false;
    this.mode = 'relative';
    this.view = F6_VIEW.slice();
    this.trail = [];
    this.marks = [];
    this.g2 = false;       // genus 2: draw F_3 alone
    this.curves = null;    // genus 2: { nodal: [pts, ...], open: [pts, ...] }, the boundary of S^2
    this._bind();
    new ResizeObserver(() => this.draw()).observe(canvas);
    if (document.fonts) document.fonts.load("14px 'EB Garamond'", 'τ0123').then(() => this.draw(), () => {});
  }

  set(opts) {
    Object.assign(this, opts);
    if (this.tau) this._fit(this.tau);
    this.draw();
  }

  /** The window on all of F_6 (or F_3). */
  showAll() { this.view = F6_VIEW.slice(); }

  /** A window of half-width R centred on tau (by default a local one, scaled with Im tau). */
  centre(tau, R = 0.2 * tau[1]) {
    this.view = [tau[0] - R, tau[0] + R, tau[1] - R, tau[1] + R];
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
    const pad = { l: 24, r: 10, t: 10, b: 30 };
    const s = Math.min(w - pad.l - pad.r, h - pad.t - pad.b);
    return { w, h, s, x0: pad.l, y0: pad.t + s };
  }
  _toPx([x, y]) {
    const { s, x0, y0 } = this._geom(), [a, b, c, d] = this.view;
    return [x0 + ((x - a) / (b - a)) * s, y0 - ((y - c) / (d - c)) * s];
  }
  _fromPx(px, py) {
    const { s, x0, y0 } = this._geom(), [a, b, c, d] = this.view;
    const cl = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
    const x = cl(a + ((px - x0) / s) * (b - a), a, b), y = cl(c + ((y0 - py) / s) * (d - c), c, d);
    return [x, Math.max(y, 1e-3 * (d - c), 1e-6)]; // stay in the upper half plane
  }
  _scale() { const { s } = this._geom(); return s / (this.view[1] - this.view[0]); }

  _bind() {
    const c = this.canvas;
    const pos = (e) => { const r = c.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; };
    const emit = (e, dragging) => {
      const tau = this._fromPx(...pos(e));
      this.target = tau;
      this.blocked = false;
      this.draw();
      this.onTarget(tau, dragging);
    };
    c.addEventListener('pointerdown', (e) => {
      if (!this.tau) return;
      if (e.shiftKey) {
        this.pan = { at: pos(e), view: this.view.slice() };
        c.setPointerCapture(e.pointerId);
        e.preventDefault();
        return;
      }
      this.drag = true;
      c.setPointerCapture(e.pointerId);
      e.preventDefault();
      emit(e, true);
    });
    c.addEventListener('pointermove', (e) => {
      if (!this.tau) { c.style.cursor = 'default'; return; }
      if (this.pan) {
        const [x, y] = pos(e), { s } = this._geom(), [a, b, cc, d] = this.pan.view;
        const dx = ((x - this.pan.at[0]) / s) * (b - a), dy = ((y - this.pan.at[1]) / s) * (d - cc);
        this.view = [a - dx, b - dx, cc + dy, d + dy];
        c.style.cursor = 'move';
        this.draw();
        return;
      }
      c.style.cursor = this.drag ? 'grabbing' : e.shiftKey ? 'move' : 'crosshair';
      if (this.drag) emit(e, true);
    });
    const up = (e) => {
      if (this.pan) { this.pan = null; return; }
      if (!this.drag) return;
      this.drag = false;
      emit(e, false);
    };
    c.addEventListener('pointerup', up);
    c.addEventListener('pointercancel', up);
    c.addEventListener('wheel', (e) => {
      if (!this.tau) return;
      e.preventDefault();
      const f = Math.exp((e.deltaMode === 1 ? e.deltaY * 16 : e.deltaY) * 0.0015);
      const [a, b, cc, d] = this.view, [px, py] = this.tau;
      this.view = [px + (a - px) * f, px + (b - px) * f, py + (cc - py) * f, py + (d - py) * f];
      this.draw();
    }, { passive: false });
    // double-click: all of F_6, or back to a window around the point
    c.addEventListener('dblclick', (e) => {
      if (!this.tau) return;
      e.preventDefault();
      const all = this.view.every((v, i) => Math.abs(v - F6_VIEW[i]) < 1e-9);
      if (all) this.centre(this.tau); else this.view = F6_VIEW.slice();
      this._fit(this.tau);
      this.draw();
    });
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
      accent: col('--accent', '#9E3F24'), shade: col('--shade', '#EAE4D6'), paper: col('--paper', '#F2EEE3'),
    };
    const FONT = "'EB Garamond', Georgia, serif";
    const { w, h, s, x0, y0 } = this._geom();
    ctx.clearRect(0, 0, w, h);
    const on = !!this.tau;
    const A = on ? 1 : 0.45;
    ctx.globalAlpha = A;

    ctx.save();
    ctx.beginPath(); ctx.rect(x0, y0 - s, s, s); ctx.clip();
    const k = this._scale();
    const P = (p) => this._toPx(p);
    const arc = (cx, r, from, to) => { const [px, py] = P([cx, 0]); ctx.arc(px, py, r * k, from, to); };
    const PI = Math.PI;

    // below the real axis: not part of the upper half plane
    const [, yAxis] = P([0, 0]);
    ctx.fillStyle = C.shade;
    if (yAxis < y0) ctx.fillRect(x0, yAxis, s, y0 - yAxis);

    // with a lattice: shade the plane outside the domain for the invariant (F_6, F_3 or F_1)
    if (this.mode !== 'relative') {
      const [, yTop] = P([0, this.view[3] + 1]);
      ctx.fillRect(x0, y0 - s, s, Math.max(0, Math.min(yAxis, y0) - (y0 - s)));
      ctx.fillStyle = C.paper;
      ctx.beginPath();
      if (this.mode === 'F6') {
        ctx.moveTo(...P([-1, 0]));
        arc(-0.5, 0.5, PI, 2 * PI);
        arc(0.5, 0.5, PI, 2 * PI);
        ctx.lineTo(P([1, 0])[0], yTop);
        ctx.lineTo(P([-1, 0])[0], yTop);
      } else if (this.mode === 'F3') {
        ctx.moveTo(...P([-1, 0]));
        arc(0, 1, PI, 2 * PI);
        ctx.lineTo(P([1, 0])[0], yTop);
        ctx.lineTo(P([-1, 0])[0], yTop);
      } else {
        ctx.moveTo(...P([-0.5, Math.sqrt(3) / 2]));
        arc(0, 1, (4 * PI) / 3, (5 * PI) / 3);
        ctx.lineTo(P([0.5, 0])[0], yTop);
        ctx.lineTo(P([-0.5, 0])[0], yTop);
      }
      ctx.closePath();
      ctx.fill();
    }

    const [, yHi] = P([0, this.view[3] + 1]);
    const r3 = Math.sqrt(3) / 2;
    if (this.g2) {
      // the three copies of F_1 in F_3
      ctx.strokeStyle = C.rule;
      ctx.lineWidth = 0.9;
      ctx.beginPath();
      for (const x of [-0.5, 0.5]) { ctx.moveTo(...P([x, r3])); ctx.lineTo(P([x, 0])[0], yHi); } // x = ±1/2
      ctx.moveTo(...P([0.5, r3])); arc(1, 1, (4 * PI) / 3, 1.5 * PI);  // |τ − 1| = 1, from e^{iπ/3} to 1 + i
      ctx.moveTo(...P([-1, 1])); arc(-1, 1, 1.5 * PI, (5 * PI) / 3);   // |τ + 1| = 1, from −1 + i to e^{2iπ/3}
      ctx.stroke();
      // the boundary of F_3
      ctx.strokeStyle = this.mode === 'relative' ? C.muted : C.ink;
      ctx.lineWidth = 1.1;
      ctx.beginPath();
      for (const x of [-1, 1]) { ctx.moveTo(...P([x, 0])); ctx.lineTo(P([x, 0])[0], yHi); }
      ctx.moveTo(...P([-1, 0])); arc(0, 1, PI, 2 * PI);
      ctx.stroke();
    } else {
      // the six copies of F_1 in F_6
      ctx.strokeStyle = C.rule;
      ctx.lineWidth = 0.9;
      ctx.beginPath();
      arc(0, 1, PI, 2 * PI);                                           // |τ| = 1
      ctx.moveTo(...P([0, 0])); arc(1, 1, PI, 1.5 * PI);               // |τ − 1| = 1, from 0 to 1 + i
      ctx.moveTo(...P([-1, 1])); arc(-1, 1, 1.5 * PI, 2 * PI);         // |τ + 1| = 1, from −1 + i to 0
      for (const x of [-0.5, 0.5]) { ctx.moveTo(...P([x, 0.5])); ctx.lineTo(P([x, 0])[0], yHi); } // x = ±1/2
      ctx.stroke();
      // the boundary of F_6
      ctx.strokeStyle = this.mode === 'relative' ? C.muted : C.ink;
      ctx.lineWidth = 1.1;
      ctx.beginPath();
      for (const x of [-1, 1]) { ctx.moveTo(...P([x, 0])); ctx.lineTo(P([x, 0])[0], yHi); }
      ctx.moveTo(...P([-1, 0])); arc(-0.5, 0.5, PI, 2 * PI); arc(0.5, 0.5, PI, 2 * PI);
      ctx.stroke();
    }
    // genus 2: the boundary of S^2, at the nodal end (solid) and the open end (dashed)
    if (this.g2 && this.curves) {
      ctx.strokeStyle = col('--blue', '#3D5A80');
      ctx.lineWidth = 1.4;
      for (const [key, dash] of [['nodal', []], ['open', [5, 3]]]) {
        ctx.setLineDash(dash);
        for (const pts of this.curves[key]) {
          ctx.beginPath();
          // (a jump means the curve left F_3 and re-entered elsewhere: don't join across it)
          pts.forEach((p, i) => {
            const [px, py] = P(p), jump = i && Math.hypot(p[0] - pts[i - 1][0], p[1] - pts[i - 1][1]) > 0.5;
            if (i && !jump) ctx.lineTo(px, py); else ctx.moveTo(px, py);
          });
          ctx.stroke();
        }
      }
      ctx.setLineDash([]);
    }
    // the real axis
    ctx.strokeStyle = C.muted;
    ctx.lineWidth = 0.8;
    ctx.beginPath(); ctx.moveTo(x0, yAxis); ctx.lineTo(x0 + s, yAxis); ctx.stroke();

    // region labels, where the regions are big enough on screen
    if (k > 45) {
      ctx.fillStyle = C.muted;
      ctx.font = `italic ${k > 110 ? 12 : 10.5}px ${FONT}`;
      ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      for (const { p, text } of this.g2 ? F3_LABELS : LABELS) { const [px, py] = P(p); ctx.fillText(text, px, py); }
    }

    // where the flow has been, and where it stopped
    ctx.fillStyle = C.muted;
    ctx.globalAlpha = A * 0.35;
    for (const p of this.trail) { const [px, py] = P(p); ctx.fillRect(px - 1, py - 1, 2, 2); }
    ctx.globalAlpha = A;
    ctx.strokeStyle = C.ink;
    ctx.lineWidth = 1.1;
    for (const p of this.marks) {
      const [px, py] = P(p);
      ctx.beginPath(); ctx.moveTo(px - 3, py - 3); ctx.lineTo(px + 3, py + 3); ctx.moveTo(px + 3, py - 3); ctx.lineTo(px - 3, py + 3); ctx.stroke();
    }
    ctx.restore();

    ctx.strokeStyle = C.ink;
    ctx.lineWidth = 1.2;
    ctx.strokeRect(x0, y0 - s, s, s);

    // the window's extent
    const [a, b, cc, d] = this.view;
    const dig = Math.max(2, Math.min(5, 2 - Math.floor(Math.log10(b - a))));
    const f = (v) => { const t = (Math.abs(v) < 10 ** -dig ? 0 : v).toFixed(dig); return t.startsWith('-') ? '−' + t.slice(1) : t; };
    ctx.fillStyle = C.muted;
    ctx.font = `12px ${FONT}`;
    ctx.textBaseline = 'top';
    ctx.textAlign = 'left'; ctx.fillText(f(a), x0, y0 + 4);
    ctx.textAlign = 'right'; ctx.fillText(f(b), x0 + s, y0 + 4);
    ctx.textAlign = 'center';
    ctx.font = `italic 14px ${FONT}`; ctx.fillText('Re τ', x0 + s / 2, y0 + 8);
    ctx.save();
    ctx.translate(x0 - 5, y0);
    ctx.rotate(-Math.PI / 2);
    ctx.textBaseline = 'bottom';
    ctx.font = `12px ${FONT}`;
    ctx.textAlign = 'left'; ctx.fillText(f(cc), 0, 0);
    ctx.textAlign = 'right'; ctx.fillText(f(d), s, 0);
    ctx.textAlign = 'center';
    ctx.font = `italic 14px ${FONT}`; ctx.fillText('Im τ', s / 2, 0);
    ctx.restore();
    ctx.globalAlpha = 1;
    if (!on) return;

    if (this.target && Math.hypot(this.target[0] - this.tau[0], this.target[1] - this.tau[1]) > 1e-9) {
      const [px, py] = this._toPx(this.target);
      ctx.strokeStyle = this.blocked ? C.muted : C.accent;
      ctx.lineWidth = 1.3;
      ctx.beginPath(); ctx.arc(px, py, 6, 0, 2 * Math.PI); ctx.stroke();
    }
    const [px, py] = this._toPx(this.tau);
    ctx.save();
    ctx.beginPath(); ctx.rect(x0 - 6, y0 - s - 6, s + 12, s + 12); ctx.clip();
    ctx.fillStyle = C.accent;
    ctx.beginPath(); ctx.arc(px, py, this.drag ? 6 : 5, 0, 2 * Math.PI); ctx.fill();
    ctx.restore();
  }
}
