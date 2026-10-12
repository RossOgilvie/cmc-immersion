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
- Done (src/cmc/leaf.js): **the leaf and its frame.** The app keeps a frame (w_1, w_2) of the period plane
  at the current data (leaf.w): (1, i) where the leaf is entered, or an adapted lattice basis when the plane
  is rational (always for g = 2; recognised by LLL). The 𝒲 slider and the τ pad move along the leaf and
  carry the frame; any other edit starts a new leaf.
- Done: **τ pad** (tau-pad.js) under the 𝒲 slider: τ = w_2/w_1 with psi = sum_j logit|alpha_j| fixed, by
  square Newton in (alpha, w_1, w_2) with the frame's periods, gauge, τ and psi as equations. F_6 always
  drawn, its six copies of F_1 labelled by the order of the classes of Λ/2Λ by length. With a lattice
  basis, τ is reduced into F_6 / F_3 / F_1 on release. (First version held 𝒲 fixed, Q'(0) = 0: see below.)
- Done: **Doubly periodic section** (period-pad.js): the plane in the w-plane with nearby integer period
  vectors (a box, plus LLL's vectors as triangles; periods as diamonds), the frame and the domain. Click to lock a period in
  (PlaneFlow: continuation, min-norm Gauss–Newton in alpha), a second to make the data doubly periodic
  (saturated lattice, adapted basis, invariants). Verified by shooting: zeta(z + w) = zeta(z) to 1e-10.
- Open, doubly periodic: edit the domain in the period picture (it is drawn there read-only), e.g. snap
  its sides to lattice vectors, which would merge with the z-domain controls and with "Close up"; the
  closing conditions (iv), (v) on top of a locked lattice; the nodes (insert a node with rational nodal
  B-periods, open it along the leaf: every invariant triple, Remark "bper nodes" of HOS); continuing the
  τ pad and the τ-flow through a node on S^1 (the paper shows the leaf passes through U^g).
- Done: the lock flow stops when two branch points come within 1% of |alpha| of each other, and when one
  step takes over 250 ms (near a collision the period quadrature refines without bound: a genus-5 lock
  used to spend 5 s per frame there); the note says why it stopped.
- Why psi and not 𝒲 for the τ pad: the slice 𝒲 = const meets each τ-fibre (an interval) only where 𝒲's
  narrow, often non-monotone range allows, so the pad soon reached a fibre's end. Measured (genus 2, 14
  random fibres): psi is strictly monotone on every fibre, from the end at 0 to the end at the circle, so
  psi = const is a global slice and (τ, psi) global coordinates; holding psi reached targets 0.3–0.8 away
  where 𝒲 stalled after 0.002–0.16. In genus 3, 4 most fibres end at the circle at both ends and psi
  turns once (2–4x better than 𝒲, but it still stops, at a fold of the slice).
- **Future: a better control function for the τ pad** (there might not exist one). Wanted: a function on
  the leaf that is monotone along every τ-fibre with range all of R, or some other global slice. Candidate
  for g >= 3: the middle of each fibre, the maximum along the fibre of a barrier
  B = sum log r_j + sum log(1 - r_j) + sum log|alpha_i - alpha_j| (it tends to -inf at every kind of end,
  so the maximum exists; uniqueness unknown; the fibre direction comes from the Newton Jacobian). Also:
  starting near one end of a fibre puts the whole psi-slice near that end.
- Done: points of the period picture are shaded by the **predicted move** |Δalpha| of one linearised
  min-norm step of the lock flow (predictMoves: one linearisation of the plane, 2g period evaluations,
  then a small least-squares problem per point). Measured in genus 3–5: within a few percent of the
  actual moves; the locks that failed (collision, 0) were predicted at 0.56–0.64, all successes <= 0.38.
  Points predicted over MOVE_MAX = 0.4 (app.js) are not shown: to be tuned with experience.
- Open: the τ pad stops at folds; past one needs path-following on the slice (as for the root-preserving
  pad). The paper's formula for dτ/dt has a sign slip (c_2 b_1 - b_2 c_1 = -aQ).
- Done: panel order Spectral data, Plane of differentials (was: Doubly periodic), Whitham deformation: all
  configurations, then a period condition, then the deformations that preserve it. The Domain section is
  gone: width, height, curvature-line alignment, Close up and Reset domain are in the Plane of
  differentials; x0, y0 (the translations) in Isospectral deformation. Reset domain became a reset
  button on each domain slider (and the angle). Open: revisit the behaviour of Align grid.
- Done: **genus 2 τ pad** draws F_3 alone (three copies of F_1), starts on all of F_3, and draws the
  boundary of S^2 (genus2.js; HOS genus 2 paper): at the nodal end t = t_0(k) with the paper's Re/Im τ,
  at the open end the real locus of f_0 = i (K - E)/E'. Complete elliptic integrals of complex modulus
  by the AGM. Checked against the app's periods: nodal end to 7e-10 (so the paper's Re τ integral is
  right; its TODO doubting the part along the real axis can go); open end to 1e-3 at s = 1e-3 (converging),
  with real δ_2/δ_1 on the locus. Note for the paper: the open-end τ is τ_0/(1 + τ_0), τ_0 = iK(α)/K(α'),
  for the basis dual to B_1^+, B_2^+ (the two differentials written down there have periods on B_1^+ and
  B_1^+ - B_2^+); iK/K' itself is another sheet of F_3 over F_1.
- Done: "Close up" (shooting along the two grid lines for a period of zeta, then the Euclidean motion,
  closing after q <= 24 periods) replaced by "Rescale domain to lattice": the smallest rectangle (same
  centre and grid angle) enclosing the parallelogram of the frame w_1, w_2. closingInfo stays in cmc.js
  (and its tests).
- **Future: a lattice-aware close up.** With a lattice Γ: take a reduced basis, check the closing
  conditions on small multiples of the basis vectors (no translational period, p_w(λ_0) = 0; no
  rotational period, μ_w(λ_0) = ±1, via the Sym integrals), find the torus sublattice Γ' ⊂ Γ, and set the
  domain to its fundamental parallelogram. Needs non-rectangular (parallelogram) domains in the mesh; keep
  the shooting along the grid lines (closingInfo) as the fallback without a lattice; for the
  root-preserving flow, the torus lattice at rational phi. Also: a general torus search (all the closing
  conditions in the UI).
- Root-preserving flow: done for even g (g >= 4 with the period plane as 2(g - 2) extra equations, a basis
  orthonormal in the plane, and a pad that marks where the flow stops). phi folds over within ~0.02 pi of
  the anchor in the examples tried, so the pad can't pass folds: continuing past one needs path-following
  on S^g itself (e.g. arclength along the fold, or a pad in other coordinates). Open: an integer basis of the period plane when it is rational (so grid points are tori again); the
  flow off the main thread (a g = 4 step takes ~100 ms); an extension off S^g, ideally to all of H^g.
- **Root-preserving flow for odd g.** Over the unit circle Σ_R is two circles swapped by σ, so the
  primitive q of Lemma 1.2 (CKKS) with σ*q = −q is not determined: there is no path in Σ_R from σ(y) to y,
  the integration constant on each circle is free, and each circle is a cycle with a possibly nonzero
  period. Going through a branch point (μ = ±1) gives φ modulo periods: an angle mod π when the periods
  are integral (rational plane, integer basis), nothing canonical otherwise. Needs a new normalisation.
- Done: the differentials b_1, b_2 (for 𝒲 and the arrows in the plane) follow the locked periods: one
  lock w gives (w, i w); two that are a basis of Λ give themselves (the second negated for orientation if
  needed); otherwise a locked period primitive in Λ (or the primitive part of the first) and the
  lattice vector completing it to a basis, with |Re b_2/b_1| <= 1/2 (first version: (w_1, τ w_1), which
  is not a lattice basis). With a lattice they are a constant real combination A of the adapted
  basis, which the τ pad keeps, so its τ stays the invariant; 𝒲 scales by det A.
- **Better handling of branch points near λ = 0.** The open ends of the τ-fibres (all alpha -> 0) are
  half of the genus-2 picture, but the flows stop at |alpha| = 0.02 (R_MIN in leaf.js, whitham.js), the
  period quadrature refines without bound there, and the surfaces bubble (u unbounded). Ideas: work in the
  rescaled coordinate lambda^+ = lambda / t, t = sqrt(alpha_1 alpha_2 ...) of the genus 2 paper (the
  roots of a^+ stay of size 1, and the differentials rescale by t^{1/2}), so that the open end H^{2+} is
  a boundary at t = 0 that the flows can reach and the periods stay cheap; then the rescaled surface
  (the blow-up view below).
- Blow-up view of the α → 0 corner (Carberry–Klein–Schmidt, arXiv:2110.01574): zoom into one bubble
  and rescale z and f, once the normalisation is settled.
- Bubbleton presets: the 𝒲 tag shows NaN at a double branch point (show "–" and disable the slider).
- **Soul conjecture.** First version done: "Set H = 1 / min |alpha - lam0|" in the Mean curvature section,
  with the domain enlarged by the same factor. Open: other scalings, and more mesh resolution as the
  domain grows. Original note: rescale the mean curvature automatically by a function of the distance from the
  branch points to lam0. Hypothesis: as one or more alpha -> lam0 the surface degenerates to a cylinder,
  but with H scaled correctly the cylinder's width shrinks and the limit is an interesting space curve.
- **Soul conjecture II.** The same in a coordinate kappa adapted to the Sym point (from our preliminary
  work): lam0 at kappa = 0 and the unit circle mapped to R. Rescaling kappa as the branch points come to
  the Sym point should give the spectral data of a curve gamma: R -> R^3, as for the NLS system.
- The algebraic Whitham ODE of the lemma (not needed for speed yet: 60–170 ms per family).
- The $S^3$ version (two Sym points), and a pen-and-ink hatching style like the example page.
- Direct dragging of divisor points (this needs an inverse map from divisor to $\tau$).
