// Runs the surface computation off the main thread.
// Message in:  { id, params }   (params as for computeSurface)
// Message out: { id, result } or { id, error }

import { computeSurface } from './cmc/cmc.js';

self.onmessage = (e) => {
  const { id, params } = e.data;
  try {
    const r = computeSurface(params);
    self.postMessage({ id, result: r }, [r.pos.buffer, r.nrm.buffer, r.u.buffer, r.st.buffer]);
  } catch (err) {
    self.postMessage({ id, error: String(err && err.message ? err.message : err) });
  }
};
