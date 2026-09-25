// Runs the surface computation off the main thread.
// Message in:  { id, params }   (params as for computeSurface), { id, kind: 'close', params },
//              or { id, kind: 'family', params: { alphas, smax, maxSteps } } (the Whitham curve)
// Message out: { id, result } or { id, error }

import { computeSurface, closingInfo } from './cmc/cmc.js';
import { WhithamCurve } from './cmc/whitham.js';

self.onmessage = (e) => {
  const { id, params, kind } = e.data;
  if (kind === 'family') {
    try {
      const t0 = performance.now();
      const tr = new WhithamCurve(params.alphas).trace(params.smax, params.maxSteps);
      self.postMessage({ id, result: { ...tr, ms: performance.now() - t0 } });
    } catch (err) {
      self.postMessage({ id, error: String(err && err.message ? err.message : err) });
    }
    return;
  }
  if (kind === 'close') {
    try {
      self.postMessage({ id, result: { s: closingInfo(params, 's'), t: closingInfo(params, 't') } });
    } catch (err) {
      self.postMessage({ id, error: String(err && err.message ? err.message : err) });
    }
    return;
  }
  try {
    const r = computeSurface(params);
    self.postMessage({ id, result: r }, [r.pos.buffer, r.nrm.buffer, r.u.buffer, r.st.buffer]);
  } catch (err) {
    self.postMessage({ id, error: String(err && err.message ? err.message : err) });
  }
};
