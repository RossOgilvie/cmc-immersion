# CMC surfaces from spectral data

An interactive viewer for constant-mean-curvature (H = ½) immersions of the plane in R³,
computed from spectral data (branch points αᵢ, Sym point λ₀, isospectral times τ) by integrating
polynomial Killing fields. The mathematics is in `cmc_spectral_to_immersion.md`; the plan in `PLAN.md`.

## Running

No build step. Serve the folder statically and open `index.html`:

    python3 -m http.server 8000      # then visit http://localhost:8000

(Opening the file directly doesn't work: browsers block module workers on `file://`.)

## Controls

- **Viewer:** drag rotates, shift+drag rolls, wheel or `+`/`−` zooms, `WASD` moves in the screen
  plane, `Q`/`E` toward/away, `F` frames the surface, `R` resets. Touch: drag, pinch, twist.
- **λ-plane:** drag branch points (dots) and the Sym point (diamond); shift+drag turns everything;
  double-click adds or removes a branch point. Blue rings are the eigenline divisor at the domain centre.
- The state lives in the URL hash, so links reproduce the surface.

## Layout

    src/cmc/cmc.js         Killing fields, isospectral flows, frame integration, Sym formula
    src/cmc/poly.js        complex polynomials and roots
    src/worker.js          computation off the main thread
    src/viewer.js          three.js scene, controls, shader (two-sided, curvature lines, colour maps)
    src/spectral-widget.js λ-plane editor
    src/app.js             panel, scheduling (coarse preview while dragging, then full), presets, export
    test/                  node tests against cmc_reference.py:  python3 test/fixtures.py && node test/run.mjs
    vendor/                three.js r186
    old/                   the earlier theta-function attempt
