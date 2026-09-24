# From spectral data to a CMC immersion of the plane

Mathematical specification for an interactive renderer. The target is a map $f:\Omega\subset\mathbb C\to\mathbb R^3$ with constant mean curvature $H=\tfrac12$, computed from a real hyperelliptic spectral curve, a point of its isospectral set, and a Sym point. All formulas in Sections 1–8 have been checked numerically with the accompanying `cmc_reference.py`, and Section 9 has been partially checked (see there). That script is ground truth for porting.

Method: polynomial Killing fields (Pinkall–Sterling 1989; in the $\lambda$-conventions of Kilian–Schmidt and Burstall–Ferus–Pedit–Pinkall), integrated as commuting ODEs. There are no theta functions and no period matrices, except for the closing conditions in Section 9, which are not needed for rendering planes.

---

## 0. Pipeline summary

1. **Input.**
   - Genus $g\ge 0$.
   - Branch points $\alpha_1,\dots,\alpha_g$ in the punctured unit disc.
   - Sym point $\lambda_0\in S^1$.
   - Isospectral times $\tau\in\mathbb R^g$.
2. **Base Killing field.** Build $a(\lambda)$ and the base Killing field $\zeta^{\rm base}$ in closed form (Section 4).
3. **Isospectral motion.** Flow $\zeta^{\rm base}$ by the $g$ isospectral flows for times $\tau$ to get $\zeta^{\rm init}$ (Section 5).
4. **Integrate over the grid.** Over a grid in $z=x+iy$, integrate simultaneously:
   - $\zeta$ (the $x$- and $y$-flows),
   - the frame $F$ at $\lambda_0$,
   - $G=\partial_\lambda F$ at $\lambda_0$ (Section 6).
5. **Output.** Evaluate the Sym formula for $f$, the normal $N$, and the scalar fields $u$, $K$, $\kappa_\pm$ (Section 6).

---

## 1. Conventions

- **Coordinates.** $z=x+iy$, $\partial_z=\tfrac12(\partial_x-i\partial_y)$, $\Delta=\partial_x^2+\partial_y^2$.
- **Surface normalization.** $H=\tfrac12$, conformal metric $4e^{2u}|dz|^2$, no umbilics.
- **Gauss equation** (elliptic sinh-Gordon, Pinkall–Sterling's (1.18) with $u=\omega$):
$$u_{z\bar z}+\tfrac12\sinh 2u=0\quad\Longleftrightarrow\quad \Delta u+2\sinh 2u=0 .$$
- **Hopf differential.** $Q=\langle f_{zz},N\rangle$ is constant with $|Q|=1$. It equals $-E$ in Pinkall–Sterling's notation.
- **$\mathbb R^3\cong\mathfrak{su}(2)$.** Use the orthonormal basis
$$e_1=\begin{pmatrix}0&i\\ i&0\end{pmatrix},\quad e_2=\begin{pmatrix}0&1\\ -1&0\end{pmatrix},\quad e_3=\begin{pmatrix}i&0\\ 0&-i\end{pmatrix},$$
with inner product $\langle X,Y\rangle=-\tfrac12\operatorname{tr}(XY)$. Coordinates are $x_k(X)=-\tfrac12\operatorname{tr}(Xe_k)$; take the real part numerically.
- **Matrix notation.** $E_{12}=\begin{pmatrix}0&1\\0&0\end{pmatrix}$, $E_{21}=E_{12}^T$, and ${}^\dagger$ denotes conjugate transpose.
- **Reflection of Laurent polynomials.** For $X(\lambda)=\sum_p X_p\lambda^p$ with $X_p\in\mathfrak{sl}(2,\mathbb C)$:
$$R(X)(\lambda):=-X(1/\bar\lambda)^\dagger=\sum_p(-X_p^\dagger)\lambda^{-p}.$$
  Consequently $X$ is $\mathfrak{su}(2)$-valued on $S^1$ iff $R(X)=X$.

---

## 2. Spectral data

**Branch points.** $\alpha_1,\dots,\alpha_g$ satisfy $0<|\alpha_i|<1$ and are pairwise distinct for a smooth curve (see Section 10 for degenerations). Define
$$a(\lambda)=-\prod_{i=1}^g(\lambda-\alpha_i)(1-\bar\alpha_i\lambda),\qquad \deg a=2g.$$
Write $a_0=a(0)=-\prod(-\alpha_i)$ and $a_{2g}$ for the leading coefficient, so $a_{2g}=\bar a_0$ and $|a_0|=\prod|\alpha_i|$.

Reality: $\overline{a(1/\bar\lambda)}=\lambda^{-2g}a(\lambda)$, and $\lambda^{-g}a(\lambda)=-\prod|\lambda-\alpha_i|^2<0$ on $S^1$. For $g=0$, $a\equiv-1$.

**Spectral curve.**
$$\Sigma:\ \nu^2=\lambda\,a(\lambda),$$
branched over $0,\infty,\alpha_i,1/\bar\alpha_i$, of genus $g$.

**Normalization constant.**
$$\kappa_0:=2\sqrt{|a_0|}=2\prod_i|\alpha_i|^{1/2}.$$

**Sym point.** $\lambda_0\in S^1$, a free parameter. Rotating all $\alpha_i$ and $\lambda_0$ by a common phase gives the same surface up to a rotation of the $z$-plane, so one real parameter is redundant. Moving $\lambda_0$ alone traverses the associated family.

**Eigenline divisor.** A point of the isospectral set (a real $g$-torus) is parametrized as in Section 5.

---

## 3. Polynomial Killing fields

A **polynomial Killing field** is
$$\zeta(\lambda)=\sum_{k=-1}^{g}\zeta_k\lambda^k=\begin{pmatrix}A&B\\ C&-A\end{pmatrix},$$
with the following structure.

**Degrees.**
- $A=\sum_{k=0}^{g-1}A_k\lambda^k$ (absent when $g=0$),
- $B=\sum_{k=-1}^{g-1}B_k\lambda^k$,
- $C=\sum_{k=0}^{g}C_k\lambda^k$.

So $\zeta_{-1}=B_{-1}E_{12}$ and $\zeta_g=C_gE_{21}$ are nilpotent.

**Reality.**
$$\zeta_k=-\zeta_{g-1-k}^\dagger\qquad(k=-1,\dots,g),$$
equivalently $A_k=-\bar A_{g-1-k}$ and $C_k=-\bar B_{g-1-k}$. This says $\lambda^{-(g-1)/2}\zeta\in\mathfrak{su}(2)$ on $S^1$.

**Spectral curve.**
$$-\det\zeta=A^2+BC=\frac{a(\lambda)}{\lambda}.$$
The eigenvalues of $\zeta(\lambda)$ are $\pm\nu/\lambda$, and a point $(\lambda,\nu)\in\Sigma$ corresponds to the eigenline $\ker(\zeta(\lambda)-\nu/\lambda)$.

**Lax pair.** The Lax pair is read off from $\zeta$. With $\pi_0(M)=\begin{pmatrix}\tfrac12M_{11}&0\\ M_{21}&\tfrac12M_{22}\end{pmatrix}$:
$$U(\lambda)=\kappa_0^{-1}\big(\lambda^{-1}\zeta_{-1}+\pi_0(\zeta_0)\big)=\kappa_0^{-1}\begin{pmatrix}\tfrac12A_0&\lambda^{-1}B_{-1}\\ C_0&-\tfrac12A_0\end{pmatrix},$$
$$V(\lambda)=R(U)(\lambda)=\kappa_0^{-1}\begin{pmatrix}\tfrac12A_{g-1}&B_{g-1}\\ \lambda C_g&-\tfrac12A_{g-1}\end{pmatrix}.$$
For $g=0$, drop the $A$ terms and read $B_{g-1}=B_{-1}$, $C_g=C_0$.

**Equations of motion.**
$$\zeta_z=[\zeta,U],\qquad \zeta_{\bar z}=[\zeta,V],$$
or in real form
$$\zeta_x=[\zeta,\,U+V],\qquad \zeta_y=[\zeta,\,i(U-V)].$$
These vector fields preserve the degree structure, the reality condition and $\det\zeta$, and they commute (Adler–Kostant–Symes). A frame $F(z,\lambda)\in SL(2,\mathbb C)$ solving
$$F_x=F\,(U+V),\qquad F_y=F\,i(U-V),\qquad F(0)=I$$
exists and is $SU(2)$-valued for $|\lambda|=1$.

**Recovered geometry.**
$$e^{u}=\frac{2|B_{-1}|}{\kappa_0}\qquad\text{(equivalently } e^{2u}=|B_{-1}/C_0|\text{)} .$$
This $u$ solves $\Delta u+2\sinh 2u=0$. Only $\zeta_{-1},\zeta_0$ (and, by reality, $\zeta_{g-1},\zeta_g$) enter $U,V$; the higher coefficients are what make the system closed.

---

## 4. Constructing an initial Killing field (Fejér–Riesz)

Write $b(\lambda):=\lambda B(\lambda)$, a polynomial of degree $g$. The reality condition forces $C=-b^\sharp$, where
$$b^\sharp(\lambda):=\lambda^g\,\overline{b(1/\bar\lambda)}$$
(reverse and conjugate the coefficients). The determinant condition becomes
$$b\,b^\sharp=q,\qquad q(\lambda):=\lambda A(\lambda)^2-a(\lambda),$$
an identity between polynomials of degree $\le 2g$.

**Algorithm.**

1. **Choose $A$.** Pick $A$ with $A_k=-\bar A_{g-1-k}$; this is $g$ real parameters. It must satisfy
   $$|A(\lambda)|\le\prod_i|\lambda-\alpha_i|\quad\text{on }|\lambda|=1,$$
   which is exactly $\lambda^{-g}q\ge0$ on $S^1$. With strict inequality, $q$ has no roots on $S^1$.
2. **Roots of $q$.** They come in pairs $\rho_j,\,1/\bar\rho_j$ with $|\rho_j|<1$, for $j=1,\dots,g$.
3. **Form $b$.** Set $b(\lambda)=\beta\prod_j(\lambda-\rho_j)$, where
   - $|\beta|^2=q_{2g}\prod_j(-1/\bar\rho_j)$, which is automatically real and positive;
   - $\arg\beta$ is chosen so that $b(0)>0$. This is a gauge choice: a constant diagonal $U(1)$ conjugation, which rotates the surface in $\mathbb R^3$.
4. **Assemble** $\zeta=\begin{pmatrix}A&b/\lambda\\ -b^\sharp&-A\end{pmatrix}$.

Choosing the roots outside the disc instead gives other points of the isospectral set; they are not needed.

**Canonical base point, $A=0$.** Then $q=-a$, $\rho_j=\alpha_j$, and
$$\zeta^{\rm base}=\begin{pmatrix}0&\beta\lambda^{-1}\prod_j(\lambda-\alpha_j)\\ -\bar\beta\prod_j(1-\bar\alpha_j\lambda)&0\end{pmatrix},\qquad \beta=e^{-i\arg\prod_j(-\alpha_j)} .$$
Its eigenline divisor is the set of branch points $\{\alpha_j\}$ (see below).

**Dimension check.** $A$ carries $g$ real parameters, and the isospectral set is a real $g$-torus. The Fejér–Riesz picture covers it by $2^g$ copies of the region $\{|A|\le\prod|\lambda-\alpha_i|\}$ (one per inside/outside root choice), glued along boundaries where roots cross $S^1$. This is awkward as a UI parametrization, hence Section 5.

**Divisor readout.** For display, the eigenline divisor is defined as the points where the eigenline of $\zeta$ is $\mathbb C\,(0,1)^T$. These are
$$D=\{(\mu_j,\nu_j)\}_{j=1}^g,\qquad b(\mu_j)=0,\qquad \nu_j=-\mu_jA(\mu_j).$$
Check: $\nu_j^2=\mu_j a(\mu_j)$ because $\mu_jA(\mu_j)^2=a(\mu_j)$ when $b(\mu_j)=0$. The $\mu_j$ lie in $\mathbb C^*$ and are generally not in the unit disc, so the display must handle $|\mu_j|>1$, for example with the full plane or a sphere. $\nu_j$ says which sheet the point is on. At the base point, $A=0$ and $D=\{\alpha_j\}$.

---

## 5. The isospectral flows (parametrizing the real $g$-torus)

For $k\ge0$ and $c\in\mathbb C$, define
$$W_k=\kappa_0^{-1}\Big(\sum_{p=-1-k}^{-1}\zeta_{p+k}\lambda^{p}+\pi_0(\zeta_k)\Big)\quad(\text{the “}\pi\text{-projection” of }\kappa_0^{-1}\lambda^{-k}\zeta),$$
$$G_{k,c}=c\,W_k+R(c\,W_k),\qquad \frac{d\zeta}{dt}=[\zeta,G_{k,c}] .$$
The generator $G_{k,c}$ has powers $-1-k,\dots,1+k$, but the commutator stays within powers $-1,\dots,g$.

**The $x$ and $y$ flows.** $k=0$ reproduces them: $G_{0,1}=U+V$ gives $\partial_x$, and $G_{0,i}=i(U-V)$ gives $\partial_y$.

**Properties** (all verified):
- all $G_{k,c}$-flows commute with each other and with $x,y$;
- they preserve the degree structure, reality and $\det\zeta$;
- $G_{k,c}$ and $G_{g-1-k,c'}$ span the same 2-plane;
- for odd $g$, the flow $G_{(g-1)/2,\,1}\equiv0$. In particular, for $g=1$ the $x$-flow is trivial and $u=u(y)$.

**A basis of $g$ independent real flows:**
$$\mathcal B=\{(k,1),(k,i):0\le k<g/2\}\ \cup\ \{((g-1)/2,\,i)\ \text{if } g \text{ odd}\}.$$

**Isospectral parametrization.** Given $\tau\in\mathbb R^g$,
$$\zeta^{\rm init}=\Phi^{\mathcal B_g}_{\tau_g}\circ\cdots\circ\Phi^{\mathcal B_1}_{\tau_1}(\zeta^{\rm base}).$$
Integrate the flows one after another; order is irrelevant because they commute.

- **Redundant directions.** For $g\ge2$, the two $k=0$ directions only translate the domain (same surface, shifted origin). The shape-changing parameters are the other $g-2$. For $g=1$ the single flow is the $y$-translation, so there are no shape parameters. For $g=0$ there is no flow.
- **Range of $\tau$.** The flows are linear on the torus but generically not periodic in $\tau$, so the $\tau$ are "times", not angles. Bounded ranges such as $|\tau_j|\le 2\pi$ are a reasonable default.
- **Alternative parametrization.** Use $A$ directly (Section 4, root choice inside). This covers only part of the torus and has a constrained domain.

---

## 6. Frame, Sym formula, and derived quantities

**What to integrate.** Fix $\lambda_0\in S^1$. Along the grid, integrate the state $(\zeta,F,G)$ with $F=F(\cdot,\lambda_0)\in SU(2)$ and $G=\partial_\lambda F(\cdot,\lambda_0)$. With $X$ the generator ($U+V$ for $x$-steps, $i(U-V)$ for $y$-steps):
$$\partial\zeta=[\zeta,X],\qquad \partial F=F\,X(\lambda_0),\qquad \partial G=G\,X(\lambda_0)+F\,X'(\lambda_0).$$
The initial values are $\zeta(0)=\zeta^{\rm init}$, $F(0)=I$, $G(0)=0$. The $\lambda$-derivatives are
$$\partial_\lambda U=-\lambda^{-2}\kappa_0^{-1}B_{-1}E_{12},\qquad \partial_\lambda V=\kappa_0^{-1}C_gE_{21}.$$

**Sym–Bobenko formula** ($\lambda=e^{i\theta}$):
$$f=4\,(\partial_\theta F)F^{-1}\big|_{\lambda_0}=4i\lambda_0\,G\,F^{-1}\ \in\mathfrak{su}(2)\cong\mathbb R^3,$$
$$N=-F\,e_3\,F^{-1}.$$
This $N$ agrees with $f_x\times f_y/|f_x\times f_y|$ in the coordinates of Section 1, and $H=+\tfrac12$ with respect to it. A different initial frame $F(0)$ changes $f$ by a rigid motion.

**Derived fields** (useful for coloring and diagnostics):
- **Conformal factor:** $e^u=2|B_{-1}|/\kappa_0$; the metric is $4e^{2u}|dz|^2$.
- **Gauss curvature:** $K=\tfrac14(1-e^{-4u})$.
- **Principal curvatures:** $\kappa_\pm=\tfrac12(1\pm e^{-2u})$.
- **Hopf differential:** $Q=\langle f_{zz},N\rangle=-e^{i\arg a_0}/\lambda_0$, constant (Pinkall–Sterling torsion invariant $E=-Q$).
- **Curvature lines.** In the $z$-plane they are the straight lines in the directions $e^{i\phi}$ with $\phi=-\tfrac12\arg Q$ and $\phi+\pi/2$. Rotating the parameter grid by $\phi$ gives a curvature-line parametrization, as in Pinkall–Sterling's pictures. The $x,y$ lines are themselves curvature lines iff $\lambda_0=\pm e^{i\arg a_0}$.

---

## 7. Grid integration

- **Integration order.** The flows commute, so the result is path independent. Integrate along the $x$-axis from $0$, then along each vertical line. Other orders are fine too, for example outward from the origin to the four quadrants, or sharing work across a tree.
- **Integrator.** RK4 with step about $0.01$–$0.03$ is ample for moderate $|\alpha_i|$ and domain size. The system is small: $2(g+2)$ matrices for $\zeta$ plus two $2\times2$ matrices for $F$ and $G$.
- **Invariants** to monitor or project back onto:
  - $\det\zeta=-a/\lambda$;
  - the reality condition $\zeta_k=-\zeta_{g-1-k}^\dagger$, which can be re-imposed by $\zeta_k\leftarrow\tfrac12(\zeta_k-\zeta_{g-1-k}^\dagger)$;
  - $F\in SU(2)$.
- **Domain scale.**
  - $\Omega$ carries metric $4e^{2u}|dz|^2$, where $u$ is bounded but can be large when some $|\alpha_i|$ is small. A fixed $z$-domain can then give very large or small surfaces; normalizing the view by the bounding box is advisable.
  - Typical features have size $O(1)$–$O(10)$ in $z$. The Delaunay period in $y$ is about $3$ for $|\alpha|=0.5$.
- **Changing the Sym point.** When $\lambda_0$ changes, $\zeta$ along the grid is unchanged; only $F$ and $G$ must be recomputed. Caching $\zeta$ on the grid is therefore worthwhile.
- **Changing $\tau$ or the $\alpha$.** This requires recomputing everything.

---

## 8. Test cases (with expected results)

1. **$g=0$ (round cylinder).** $a=-1$, $\zeta=\begin{pmatrix}0&\lambda^{-1}\\-1&0\end{pmatrix}$, $\kappa_0=2$, $u\equiv0$. The surface is a cylinder of radius $1$ for every $\lambda_0$.
2. **$g=1$ (Delaunay).** Take $\alpha=r e^{i\varphi}$ and the base point.
   - **Generic Sym point.** $u=u(y)$ only, the $x$-curves are helices, and the surface is a helicoidal ("twizzler") CMC surface.
   - **$\lambda_0=-\alpha/|\alpha|$: unduloid.** Neck and bulge radii $\dfrac{2r}{1+r},\ \dfrac{2}{1+r}$.
   - **$\lambda_0=+\alpha/|\alpha|$: nodoid.** Radii $\dfrac{2r}{1-r},\ \dfrac{2}{1-r}$.

   The $x$-curves are the parallels. As $r\to1$ the unduloid tends to the unit cylinder; as $r\to0$ it tends to a chain of spheres.
3. **Generic invariants for any $g$ and data.** Finite-difference checks on the output mesh:
   - $|f_x|=|f_y|=2e^u$ and $f_x\cdot f_y=0$;
   - $H=\tfrac12$;
   - $Q=-e^{i\arg a_0}/\lambda_0$;
   - $\Delta u+2\sinh 2u=0$.
4. **Commutativity.** Integrating around a closed rectangle in $z$ returns $\zeta$ and $F$ to their start values up to integrator error.
5. **Base-point divisor.** At $A=0$, the divisor readout returns exactly $\{\alpha_j\}$.

`cmc_reference.py` runs checks 1–3 and 5.

---

## 9. Closing conditions (for later: cylinders and tori)

Not needed for rendering planes. Stated here in the same conventions.

### 9.1 The differentials $\Theta_w$

For a displacement $w\in\mathbb C$ in the $z$-plane, let $\Theta_w$ be the meromorphic differential on $\Sigma$ defined by the following:

- **Form.** $\Theta_w=\dfrac{p_w(\lambda)\,d\lambda}{\lambda\,\nu}$ with $\deg p_w\le g+1$. Such differentials are odd under $\nu\mapsto-\nu$, holomorphic away from $0,\infty$, and have double poles (no residues) there.
- **Principal parts.**
  $$\Theta_w\sim d\!\left(\frac{w}{\kappa_0}\frac{\nu}{\lambda}\right)\ \text{at }\lambda=0,\qquad \Theta_w\sim d\!\left(\frac{\bar w}{\kappa_0}\frac{\nu}{\lambda^{g}}\right)\ \text{at }\lambda=\infty .$$
  Equivalently, $p_w(0)=-\dfrac{w\,a_0}{2\kappa_0}$ and $[\lambda^{g+1}]p_w=\dfrac{\bar w\,a_{2g}}{2\kappa_0}$.
- **Real normalization.** $\operatorname{Re}\oint_\gamma\Theta_w=0$ for all $\gamma\in H_1(\Sigma,\mathbb Z)$. This is $2g$ real linear conditions on the $g$ complex middle coefficients of $p_w$, so $\Theta_w$ is unique.

**Meaning.** If $w$ is a period of $\zeta$, then $\Theta_w=d\log\mu_w$, where $\mu_w$ is the eigenvalue of the monodromy $M_w(\lambda)=F(z+w,\lambda)F(z,\lambda)^{-1}$ on the eigenline corresponding to the point of $\Sigma$. The principal parts come from the eigenvalue of $U\approx\kappa_0^{-1}\zeta$ near $0$ and of $V\approx\kappa_0^{-1}\lambda^{1-g}\zeta$ near $\infty$.

**Linearity.** $w\mapsto\Theta_w$ is $\mathbb R$-linear: $p_w=\operatorname{Re}(w)\,p^{(1)}+\operatorname{Im}(w)\,p^{(i)}$.

### 9.2 Conditions

1. **Intrinsic periods.** $w$ is a period of $\zeta$ (hence of $u$) iff $\oint_\gamma\Theta_w\in2\pi i\mathbb Z$ for all cycles $\gamma$. The set $\Gamma$ of such $w$ is a subgroup of $\mathbb C$; $u$ is doubly periodic iff $\Gamma$ is a lattice. This is Pinkall–Sterling's "lattice plane" condition, Theorem 7.5.
2. **Closing in $\mathbb R^3$ with Sym point $\lambda_0$.** For generators $w_1,w_2$ of $\Gamma$, both of the following must hold for $j=1,2$:
   - (a) $p_{w_j}(\lambda_0)=0$, i.e. $\Theta_{w_j}$ vanishes over $\lambda_0$. This makes $\partial_\lambda M_{w_j}(\lambda_0)=0$: no translational period. By linearity, (a) for two independent $w$ is equivalent to $p^{(1)}(\lambda_0)=p^{(i)}(\lambda_0)=0$.
   - (b) $\mu_{w_j}(\lambda_0)=\pm1$, i.e. $M_{w_j}(\lambda_0)=\pm I$. Compute it as
     $$\log\mu_{w}(\lambda_0)=\int_{\alpha_1}^{\lambda_0}\Theta_{w}\ \pmod{\pi i},$$
     using that $\mu=\pm1$ at branch points. Substitute $\lambda=\alpha_1+s^2$ near the endpoint to remove the square-root singularity.

   For a closed cylinder, impose (a) and (b) for a single $w$.

   Consequence: for tori, the Sym point is no longer free; it must be a common zero on $S^1$ of $p^{(1)},p^{(i)}$. Parameter count: the Pinkall–Sterling §8 count carries over.

### 9.3 Numerics for periods

- **Cycles.** Use closed loops in the $\lambda$-plane, each encircling exactly two branch points. The loops around $\{0,\alpha_j\}$ and $\{\alpha_j,1/\bar\alpha_j\}$ for $j=1,\dots,g$ span $H_1$.
- **Continuation of $\nu$.** Track $\nu=\sqrt{\lambda a(\lambda)}$ by continuity along each loop: flip sign whenever $|\nu_n+\nu_{n-1}|<|\nu_n-\nu_{n-1}|$.
- **Quadrature.** Use the trapezoid rule on these smooth closed loops; it converges geometrically.

### 9.4 Independent check

For a candidate $w$, integrate $\zeta$, $F$ and $G$ along the straight segment from $0$ to $w$. Closing holds iff all three of the following hold:
- $\zeta(w)=\zeta(0)$,
- $F(w)=\pm I$,
- $f(w)=0$.

This uses only the ODEs of Sections 3 and 6 and should be used to validate any implementation of 9.1–9.2.

### 9.5 Verification status

In $g=1$ (Delaunay, $w=iT$ with $T$ the numerically determined $y$-period):
- the real-normalized $\Theta_w$ has periods $0$ and $-2\pi i$;
- $\Theta_w$ agrees pointwise with $d\log\mu_w$ computed from the numerical monodromy;
- the Delaunay Sym points $\pm\alpha/|\alpha|$ are exactly the zeros of $p_{w}$ for real $w$.

Not yet independently verified for $g\ge2$.

---

## 10. Degenerations and pitfalls

- **$|\alpha_i|\to1$.** Branch points collide with their reflections on $S^1$, giving a double point of $\Sigma$. The construction still works (for example, Delaunay tends to the cylinder), but the Fejér–Riesz step becomes ill-conditioned if $A\ne0$. Keep $|\alpha_i|\le1-\varepsilon$ in the UI.
- **$\alpha_i\to0$.** Then $a_0\to0$, $\kappa_0\to0$ and $u$ is unbounded; the surfaces degenerate towards sphere chains or bubbletons. Keep $|\alpha_i|\ge\varepsilon$.
- **$\alpha_i=\alpha_j$.** The curve is singular, but the base-point formulas are unaffected.
- **Branch of $\sqrt{\ }$.** Nothing in Sections 3–7 requires a square root. Only Section 9 does.
- **$\lambda_0\notin S^1$.** Not allowed for $\mathbb R^3$: $F$ would not be $SU(2)$-valued and $f$ would leave $\mathfrak{su}(2)$. Two distinct Sym points $\lambda_1\ne\lambda_2$ on $S^1$ would instead give the $S^3$ version:
  $$f=F(\lambda_1)F(\lambda_2)^{-1},\qquad H=i\,\frac{\lambda_1+\lambda_2}{\lambda_1-\lambda_2}.$$
  This is not treated here.

---

## 11. Relation to the sources

- **Pinkall–Sterling.**
  - Their $\omega$ is our $u$ (same metric $4e^{2\omega}|dz|^2$ and same Gauss equation).
  - Their torsion invariant is $E=-Q=e^{i\arg a_0}/\lambda_0$.
  - Their canonical Killing field $J$ (§6), a Laurent polynomial in $\sqrt E$, plays the role of $\zeta$; its norm $\langle J,J\rangle$ produces their polynomial $P$, whose roots correspond to our branch points under $E\leftrightarrow\lambda^{-1}$ up to a phase.
  - Their type $n$ gives spectral genus $2n$ generically.
  - Their jet-space ODE (§5, (18)) is an equivalent but less convenient coordinate system for the same flows.
- **Bobenko 1991 and Ercolani–Knörrer–Trubowitz 1993** give the same surfaces in closed form via theta functions; EKT likewise use $H=\tfrac12$. The Killing-field ODE approach here avoids computing Riemann matrices and theta functions. Those formulas could serve as a second independent check but are not transcribed here.
