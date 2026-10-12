// The plane of differentials B_a, drawn in the displacement coordinate w (the z-plane: Theta_w has
// principal part omega = p_w(0) = -a(0) w / (2 kappa_0), so this is the omega-plane up to a fixed complex
// factor). Shown:
//  - the current frame w_1, w_2 (arrows) and the domain of the surface (dashed), for orientation
//  - integer period vectors m in Z^g near the period plane, each at the w whose periods L(w) are the
//    projection of m onto the plane; the darker, the smaller the predicted move of the branch points to
//    lock m in; diamonds: m is a period (in the plane), so w is a period of the solution
//    (triangles: longer vectors found by LLL reduction, much closer to the plane; zoom out to see them)
//  - locked points (rust rings): the plane is held through them
// Clicking a point calls onPick(point); hovering calls onHover(point or null). The wheel zooms, shift-drag
// pans, double-click fits the view.

export class PeriodPad {
  constructor(canvas, { onPick, onHover }) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.onPick = onPick;
    this.onHover = onHover;
    this.points = [];   // [{ w, m, dist, move }]
    this.moveMax = 0.4; // the predicted move shaded lightest
    this.locked = [];   // [w]
    this.frame = null;  // [w1, w2]
    this.domain = null; // [corner, ...] in w
    this.hover = null;
    this.view = [-1, 1, -1, 1];
    this.enabled = false;
    this._bind();
    new ResizeObserver(() => this.draw()).observe(canvas);
  }

  set(opts) {
    Object.assign(this, opts);
    this.draw();
  }

  /** A view centred on 0 that shows the frame, the domain and the locked points. */
  fit() {
    let R = 0;
    for (const w of [...(this.frame || []), ...(this.domain || []), ...this.locked]) R = Math.max(R, Math.hypot(...w));
    const near = this.points.filter((p) => !p.lll).slice(0, 12).map((p) => Math.hypot(...p.w));
    if (near.length) R = Math.max(R, ...near);
    R = 1.15 * (R || 1);
    this.view = [-R, R, -R, R];
  }

  _geom() {
    const w = this.canvas.clientWidth, h = this.canvas.clientHeight;
    const pad = { l: 10, r: 10, t: 10, b: 22 };
    const s = Math.min(w - pad.l - pad.r, h - pad.t - pad.b);
    return { w, h, s, x0: pad.l + (w - pad.l - pad.r - s) / 2, y0: pad.t + s };
  }
  _toPx([x, y]) {
    const { s, x0, y0 } = this._geom(), [a, b, c, d] = this.view;
    return [x0 + ((x - a) / (b - a)) * s, y0 - ((y - c) / (d - c)) * s];
  }

  _nearest(px, py) {
    let best = null, bd = 9;
    for (const p of this.points) {
      const [x, y] = this._toPx(p.w), d = Math.hypot(x - px, y - py);
      if (d < bd) { bd = d; best = p; }
    }
    return best;
  }

  _bind() {
    const c = this.canvas;
    const pos = (e) => { const r = c.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; };
    c.addEventListener('pointerdown', (e) => {
      if (!this.enabled) return;
      if (e.shiftKey) {
        this.pan = { at: pos(e), view: this.view.slice() };
        c.setPointerCapture(e.pointerId);
        e.preventDefault();
        return;
      }
      const p = this._nearest(...pos(e));
      if (p) this.onPick(p);
    });
    c.addEventListener('pointermove', (e) => {
      if (!this.enabled) { c.style.cursor = 'default'; return; }
      if (this.pan) {
        const [x, y] = pos(e), { s } = this._geom(), [a, b, cc, d] = this.pan.view;
        const dx = ((x - this.pan.at[0]) / s) * (b - a), dy = ((y - this.pan.at[1]) / s) * (d - cc);
        this.view = [a - dx, b - dx, cc + dy, d + dy];
        c.style.cursor = 'move';
        this.draw();
        return;
      }
      const p = this._nearest(...pos(e));
      c.style.cursor = p ? 'pointer' : e.shiftKey ? 'move' : 'default';
      if (p !== this.hover) { this.hover = p; this.draw(); this.onHover(p); }
    });
    const up = () => { this.pan = null; };
    c.addEventListener('pointerup', up);
    c.addEventListener('pointercancel', up);
    c.addEventListener('pointerleave', () => { if (this.hover) { this.hover = null; this.draw(); this.onHover(null); } });
    c.addEventListener('wheel', (e) => {
      if (!this.enabled) return;
      e.preventDefault();
      const f = Math.exp((e.deltaMode === 1 ? e.deltaY * 16 : e.deltaY) * 0.0015);
      const [a, b, cc, d] = this.view, mx = (a + b) / 2, my = (cc + d) / 2;
      this.view = [mx + (a - mx) * f, mx + (b - mx) * f, my + (cc - my) * f, my + (d - my) * f];
      this.draw();
    }, { passive: false });
    c.addEventListener('dblclick', (e) => { if (this.enabled) { e.preventDefault(); this.fit(); this.draw(); } });
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
      accent: col('--accent', '#9E3F24'), blue: col('--blue', '#3D5A80'),
    };
    const FONT = "'EB Garamond', Georgia, serif";
    const { w, h, s, x0, y0 } = this._geom();
    ctx.clearRect(0, 0, w, h);
    ctx.globalAlpha = this.enabled ? 1 : 0.45;
    ctx.save();
    ctx.beginPath(); ctx.rect(x0, y0 - s, s, s); ctx.clip();
    const P = (p) => this._toPx(p);

    // axes through 0
    const [ox, oy] = P([0, 0]);
    ctx.strokeStyle = C.rule;
    ctx.lineWidth = 0.8;
    ctx.beginPath(); ctx.moveTo(x0, oy); ctx.lineTo(x0 + s, oy); ctx.moveTo(ox, y0 - s); ctx.lineTo(ox, y0); ctx.stroke();

    // the domain of the surface
    if (this.domain) {
      ctx.strokeStyle = C.muted;
      ctx.setLineDash([4, 3]);
      ctx.beginPath();
      this.domain.forEach((q, i) => { const [x, y] = P(q); if (i) ctx.lineTo(x, y); else ctx.moveTo(x, y); });
      ctx.closePath(); ctx.stroke();
      ctx.setLineDash([]);
    }

    // the frame (under the points, so a lattice point at its tip stays visible)
    const arrow = (wv, label) => {
      const [x, y] = P(wv), L = Math.hypot(x - ox, y - oy);
      if (L < 4) return;
      const ux = (x - ox) / L, uy = (y - oy) / L;
      ctx.beginPath(); ctx.moveTo(ox, oy); ctx.lineTo(x, y);
      ctx.moveTo(x, y); ctx.lineTo(x - 7 * ux + 3.5 * uy, y - 7 * uy - 3.5 * ux);
      ctx.moveTo(x, y); ctx.lineTo(x - 7 * ux - 3.5 * uy, y - 7 * uy + 3.5 * ux);
      ctx.stroke();
      ctx.font = `italic 14px ${FONT}`;
      ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.fillText('w', x + 11 * ux, y + 11 * uy);
      ctx.font = `10px ${FONT}`;
      ctx.fillText(label, x + 11 * ux + 6, y + 11 * uy + 4);
    };
    if (this.frame) {
      ctx.strokeStyle = C.blue; ctx.fillStyle = C.blue; ctx.lineWidth = 1.4;
      arrow(this.frame[0], '1');
      arrow(this.frame[1], '2');
    }

    // integer period vectors, darker the smaller the predicted move to lock them in (log scale from 0.005
    // to moveMax; without a prediction, by the distance to the plane)
    for (const p of this.points) {
      const [x, y] = P(p.w);
      const a = p.dist < 1e-8 ? 1
        : p.move !== undefined ? Math.min(0.85, Math.max(0.12, 0.85 - 0.73 * Math.log(Math.max(p.move, 0.005) / 0.005) / Math.log(this.moveMax / 0.005)))
          : Math.min(0.85, Math.max(0.12, -Math.log10(p.dist) / 4));
      ctx.globalAlpha = (this.enabled ? 1 : 0.45) * a;
      ctx.fillStyle = C.ink;
      ctx.beginPath();
      if (p.dist < 1e-8) { ctx.moveTo(x, y - 4.2); ctx.lineTo(x + 4.2, y); ctx.lineTo(x, y + 4.2); ctx.lineTo(x - 4.2, y); ctx.closePath(); }
      else if (p.lll) { ctx.moveTo(x, y - 3.6); ctx.lineTo(x + 3.2, y + 2.2); ctx.lineTo(x - 3.2, y + 2.2); ctx.closePath(); }
      else ctx.arc(x, y, 2.6, 0, 2 * Math.PI);
      ctx.fill();
    }
    ctx.globalAlpha = this.enabled ? 1 : 0.45;

    // locked points and the hovered one
    ctx.strokeStyle = C.accent;
    ctx.lineWidth = 1.6;
    for (const q of this.locked) { const [x, y] = P(q); ctx.beginPath(); ctx.arc(x, y, 6.5, 0, 2 * Math.PI); ctx.stroke(); }
    if (this.hover) {
      const [x, y] = P(this.hover.w);
      ctx.lineWidth = 1;
      ctx.beginPath(); ctx.arc(x, y, 5, 0, 2 * Math.PI); ctx.stroke();
    }
    ctx.restore();

    ctx.strokeStyle = C.ink;
    ctx.lineWidth = 1.2;
    ctx.strokeRect(x0, y0 - s, s, s);
    // the window's half-width
    const [a, b] = this.view;
    ctx.fillStyle = C.muted;
    ctx.font = `12px ${FONT}`;
    ctx.textAlign = 'right'; ctx.textBaseline = 'top';
    ctx.fillText(`width ${(b - a).toPrecision(3)}`, x0 + s, y0 + 4);
    ctx.textAlign = 'left';
    ctx.font = `italic 13px ${FONT}`;
    ctx.fillText('w-plane', x0, y0 + 4);
    ctx.globalAlpha = 1;
  }
}
