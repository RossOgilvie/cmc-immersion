# Implementation plan: interactive CMC immersions from spectral data

Goal: a single web page. Most of the screen is a 3D viewer of the CMC ($H=\tfrac12$) immersion
$f:\Omega\subset\mathbb C\to\mathbb R^3$. A panel on the right edits the spectral data, and the surface
recomputes live. The math follows `cmc_spectral_to_immersion.md` (§§1–8); `cmc_reference.py` is the ground truth.

---

## 1. Architecture

No build step: plain ES modules, with three.js vendored into `vendor/` (works offline, no CDN at runtime).
Serve with any static server (`python -m http.server`).

```
index.html              layout, panel markup, styles
src/
  cmc/complex.js        tiny helpers for complex numbers / 2x2 complex matrices on Float64Array
  cmc/killing.js        a(λ), κ0, base Killing field (§4), Fejér–Riesz (optional A≠0), divisor readout
  cmc/flows.js          generator G_{k,c}, flow RHS [ζ,G], isospectral motion (§5)
  cmc/surface.js        grid integration of (ζ, F, G), Sym formula, N, u, K (§§6–7)
  worker.js             runs surface.js off the main thread; latest request wins
  viewer.js             three.js scene, trackball/keyboard controls, mesh + shader material
  spectral-widget.js    interactive λ-plane editor (canvas/SVG)
  app.js                glue: state, panel controls, progressive recompute
test/
  fixtures.py           dumps JSON test data from cmc_reference.py
  run.mjs               node test runner (fixtures + geometric invariants + Delaunay radii)
```

## 2. Math core (port of §§3–7)

- **Storage.** A Laurent polynomial of 2×2 complex matrices is one flat `Float64Array` (8 doubles per coefficient),
  the same layout as the Python (`X[i] = X_{lo+i}`). No allocation in the inner loop: RK4 uses preallocated scratch buffers.
- **State per integration path:** ζ ($8(g+2)$ doubles) plus $F$ and $G=\partial_\lambda F$ at $\lambda_0$ (16 doubles).
- **Direction-general stepping.** A step in direction $e^{i\phi}$ uses the generator $\cos\phi\,(U+V)+\sin\phi\,i(U-V)$.
  So the parameter grid can be rotated by $\phi=-\tfrac12\arg Q$, which makes the grid lines **curvature lines**
  (on by default, toggleable).
- **Grid order.** Integrate from the domain centre outward: along the "horizontal" axis both ways, then each column
  up and down. This halves the path length versus corner-start, and so halves the accumulated error.
- **Step size decoupled from mesh resolution.** The RK4 step is $h\le h_{\max}$ (default 0.02, "accuracy" setting).
  The mesh samples every $m$-th point, so a 200×200 mesh over a large $z$-domain stays accurate.
- **Invariant hygiene.** After each output point: re-impose reality $\zeta_k\leftarrow\tfrac12(\zeta_k-\zeta_{g-1-k}^\dagger)$
  and re-unitarise $F$. Monitor $\det\zeta+a/\lambda$ and report the max drift in a diagnostics line.
- **Outputs per vertex:** position $f$, normal $N$ (exact, from the frame, not finite differences), $u$, $K$, the $z$-coordinate (for grid lines).
- **Isospectral parameters.** For every $g$, the two $k=0$ flows *are* the $x/y$ translations. So the UI exposes them as
  "domain centre $z_0$", and exposes only the $g-2$ genuine shape times $\tau$ as sliders (none for $g\le2$).
  $\zeta^{\rm init}$ is the base point flowed by these (§5).

## 3. Responsiveness

- All computation runs in a Web Worker. Each request carries an id, and stale results are dropped.
- **Progressive refinement.** While a handle is being dragged, compute a coarse mesh (~64×64) for immediate feedback.
  When input stops, compute the full mesh (default ~256×256, adjustable).
- Cache $\zeta^{\rm init}$ (the isospectral flow) so that moving only $\lambda_0$ or $z_0$ skips it. The spec notes that
  $\lambda_0$ changes could reuse a cached ζ-grid. I'll do this only if timing shows it's needed, because re-integrating
  everything is probably fast enough.
- Rough budget: 256² mesh, $L\approx 20$, $h=0.02$, $g=4$ is about $10^6$ RHS evaluations, which should take a few hundred ms in JS.
  Coarse previews should take a few tens of ms.

## 4. 3D viewer

- **Controls** (as in the example page): a quaternion trackball. Drag rotates about the screen axis perpendicular to
  the drag, shift+drag rolls about the view axis. Wheel or `+`/`−` zooms. `WASDQE` translates (see open question).
  `R` resets the view, and `F` re-frames the camera to the bounding box. Touch: one finger rotates, pinch zooms.
- **Framing.** Surfaces can be huge or tiny (§7, "domain scale"). The camera auto-frames the bounding box on genus
  change, preset load, or `F`. It does *not* auto-frame on every slider move, so a moving surface doesn't jitter.
- **Material** (custom `ShaderMaterial`): two-sided with distinct front/back colours (so the mean-curvature side is
  visible), Blinn–Phong lighting, and anti-aliased grid lines drawn in the fragment shader from the $z$-coordinate
  (`fwidth`). These are the curvature lines when the grid is rotated.
- **Colour modes:** plain two-sided, conformal factor $u$, Gauss curvature $K$ (diverging map centred at 0).
- **Extras:** PNG screenshot, OBJ export of the current mesh.

## 5. Spectral-data panel

- **λ-plane widget** (the main interaction): a canvas showing the unit circle and a region a little beyond it.
  - Branch points $\alpha_i$ are draggable dots in the punctured disc, clamped to $\varepsilon\le|\alpha|\le1-\varepsilon$.
    Their reflections $1/\bar\alpha_i$ are drawn faintly when on screen.
  - The Sym point $\lambda_0$ is a draggable dot on the circle.
  - The divisor readout $\mu_j$ (§4) is drawn as hollow markers, as a read-only view of where on the isospectral torus we are.
    Points outside the view are clamped to the edge with an arrow.
  - Add / remove branch point buttons (genus 0…6). Hold shift while dragging to rotate all $\alpha_i$ and $\lambda_0$ together
    (a pure $z$-rotation, per §2).
- **Numeric fields** for each $\alpha_i$ ($r,\theta$) and for $\lambda_0$ ($\theta$), for exact values.
- **Shape times** $\tau_1..\tau_{g-2}$: sliders in $[-2\pi,2\pi]$.
- **Domain:** centre $z_0$, width/height in $z$, curvature-line alignment toggle, mesh resolution, accuracy.
- **Animation:** "play" buttons that sweep $\lambda_0$ around the circle (the associated family) or sweep one $\tau$.
- **Presets:** round cylinder (g=0), unduloid, nodoid, twizzler (g=1), and a few hand-picked g=2,3 examples.
- **Diagnostics line:** compute time, invariant drift, and the $Q$ check.
- **Shareable state:** parameters are encoded in the URL hash.

## 6. Verification

1. `test/fixtures.py` dumps, for several $(g,\alpha,A,\tau,\lambda_0)$: $\zeta^{\rm base}$, $\zeta^{\rm init}$, and a small grid of $f,N,u$.
   `node test/run.mjs` compares the JS port to these values (target ~1e-10 for ζ, ~1e-8 for the grid).
2. The same test runner re-implements §8's finite-difference checks in JS: conformality, $H=\tfrac12$, $Q$, the Delaunay neck and bulge radii,
   and the base-point divisor. It also checks commutativity on a closed rectangle.
3. Visual check: headless Chromium screenshots of the presets, which I inspect myself before handing over.

## 7. Milestones

1. Math core + node tests passing against the Python reference.
2. Worker + viewer with trackball/keyboard controls, showing a hard-coded Delaunay surface.
3. λ-plane widget + panel wired up; progressive recompute.
4. Colour modes, curvature lines, presets, URL state, export, animation.
5. Polish pass: screenshots, performance tuning, mobile/narrow layout, README.

## 8. Status of later additions

- Done: H slider (a homothety in R³); genus-1 Delaunay locus in the λ-plane; bubbletons (double branch
  point at a resonance point of the cylinder); the three-lobed Wente torus (McIntosh's data, refined by
  Newton); `src/cmc/periods.js` computes the periods of Θ_w and the lattice Γ (§9.1–9.3); `closingInfo`
  does the §9.4 shooting check; metric-adapted meshes; per-step adaptive RK4 step.

## 9. Out of scope for now (noted for later)

- General torus search (closing conditions in the UI), and a "close up" that uses the lattice Γ rather
  than the grid axes (so non-rectangular fundamental domains work).
- Whitham family in the λ-plane: done (traced in a worker, chord Newton, W from the residue formula).
  Not done: the algebraic Whitham ODE of the lemma (not needed for speed yet: 60–170 ms per family).
- Whitham: done as a slider (continuation keeping [ell] in CP^{g-1} and arg a(0) fixed). Open: use the
  Willmore functional as the flow parameter, and follow the common root of B_a as a moving Sym point.
- The $S^3$ version (two Sym points), and a pen-and-ink hatching style like the example page.
- Direct dragging of divisor points (this needs an inverse map from divisor to $\tau$).
