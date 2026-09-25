# Implementation plan: interactive CMC immersions from spectral data

Goal: a single web page. Most of the screen is a 3D viewer of the CMC ($H=\tfrac12$) immersion
$f:\Omega\subset\mathbb C\to\mathbb R^3$. A panel on the right edits the spectral data, and the surface
recomputes live. The math follows `cmc_spectral_to_immersion.md` (§§1–8); the values of the (retired) Python reference are frozen in `test/fixtures.json`.

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
- **Texture modes:** plain two-sided, conformal factor $u$, Gauss curvature $K$ (diverging map centred at 0), wireframe (the parameter lines alone). A transparency slider fades the fill; then the lines are drawn first as an opaque, depth-writing pass and the fill blended over them, so lines never mis-sort against the fill.
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

1. `test/fixtures.json` holds, for several $(g,\alpha,A,\tau,\lambda_0)$: $\zeta^{\rm base}$, $\zeta^{\rm init}$, and a small grid of $f,N,u$.
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

## 9. Future work

Done so far, for reference: the Whitham family in the λ-plane (traced in a worker, chord Newton, 𝒲 from
the residue formula); the Whitham slider, parametrised by 𝒲, with its ends at the obstacles and
"Flow to common root"; the genus-2 root-preserving flow as a pad of the Sym integrals (phi_1, phi_2),
solved by Gauss–Newton with the common root pinned at lam0 (src/cmc/rootflow.js).

Open:

- **Video export of the τ-Whitham family.** Candidate formats: video (WebM/MP4 via MediaRecorder on the
  canvas), GIF, an SVG animation (SMIL or CSS keyframes over per-frame line drawings), or something else.
- **Doubly periodic case: τ and the fundamental domain.** Compute the conformal type τ of the lattice Γ
  and display the F_1, F_3 or F_6 domain.
- **Projection to the doubly periodic case.** Doubly periodic data are dense in H^g, so "the nearest"
  isn't well defined. One idea: for the current finite domain in C, give each point a distance to the
  nearest integer periods; pick the two independent points with the lowest value; step in H^g in the
  direction that moves those periods towards integral periods. Other schemes are possible, and the
  choice of distance to the integer points is open to experiment.
- General torus search (closing conditions in the UI), and a "close up" that uses the lattice Γ rather
  than the grid axes (so non-rectangular fundamental domains work; for the root-preserving flow, the
  torus lattice at rational phi).
- Root-preserving flow: done for even g (g >= 4 with the period plane as 2(g - 2) extra equations, a basis
  orthonormal in the plane, and a pad that marks where the flow stops). phi folds over within ~0.02 pi of
  the anchor in the examples tried, so the pad can't pass folds: continuing past one needs path-following
  on S^g itself (e.g. arclength along the fold, or a pad in other coordinates). Open: odd g (the definition
  of phi); an integer basis of the period plane when it is rational (so grid points are tori again); the
  flow off the main thread (a g = 4 step takes ~100 ms); an extension off S^g, ideally to all of H^g.
- Blow-up view of the α → 0 corner (Carberry–Klein–Schmidt, arXiv:2110.01574): zoom into one bubble
  and rescale z and f, once the normalisation is settled.
- Bubbleton presets: the 𝒲 tag shows NaN at a double branch point (show "–" and disable the slider).
- The algebraic Whitham ODE of the lemma (not needed for speed yet: 60–170 ms per family).
- The $S^3$ version (two Sym points), and a pen-and-ink hatching style like the example page.
- Direct dragging of divisor points (this needs an inverse map from divisor to $\tau$).
