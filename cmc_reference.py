"""
Reference implementation (numpy) of the construction in cmc_spectral_to_immersion.md.
Not optimized; intended as ground truth for testing a JS/threejs port.
Run `python3 cmc_reference.py` to execute the self-tests.

Storage convention: a Laurent polynomial sum_{p=lo}^{hi} X_p lam^p of 2x2 matrices
is an array X[0..hi-lo] with X[i] = X_{lo+i}.  The Killing field zeta has lo=-1, hi=g.
"""
import numpy as np

E1 = np.array([[0, 1j], [1j, 0]])
E2 = np.array([[0, 1], [-1, 0]], dtype=complex)
E3 = np.array([[1j, 0], [0, -1j]])


def a_poly(alphas):
    """a(lam) = -prod (lam-alpha)(1-conj(alpha) lam); coefficients by increasing power."""
    a = np.array([-1.0 + 0j])
    for al in alphas:
        a = np.convolve(a, np.convolve([-al, 1], [1, -np.conj(al)]))
    return a


def sharp(p, deg):
    """p^#(lam) = lam^deg * conj(p(1/conj lam))."""
    p = np.concatenate([p, np.zeros(deg + 1 - len(p))])
    return np.conj(p[::-1])


def killing_field(alphas, A=None):
    """Initial polynomial Killing field (Section 4). A=None -> base point A=0."""
    g = len(alphas)
    a = a_poly(alphas)
    if g == 0:
        b = np.array([1 + 0j])
        A = np.zeros(0, complex)
    else:
        A = np.zeros(g, complex) if A is None else np.asarray(A, complex)
        q = -a.copy()
        lA2 = np.concatenate([[0], np.convolve(A, A)])
        q[: len(lA2)] += lA2
        r = np.roots(q[::-1])
        rin = r[np.abs(r) < 1]
        if len(rin) != g:
            raise ValueError("A violates |A| <= prod|lam-alpha| on the unit circle")
        b = np.array([1 + 0j])
        for rho in rin:
            b = np.convolve(b, [-rho, 1])
        s = q[-1] / b[-1] / np.conj(b[0])          # = |beta|^2 (real, > 0)
        b = b * np.sqrt(s.real)
        b = b * np.exp(-1j * np.angle(b[0]))       # gauge: b(0) > 0
    Z = np.zeros((g + 2, 2, 2), complex)
    for k in range(len(A)):                        # A: powers 0..g-1
        Z[k + 1, 0, 0] += A[k]
        Z[k + 1, 1, 1] -= A[k]
    for k in range(g + 1):                         # B = b/lam: powers -1..g-1
        Z[k, 0, 1] = b[k]
    bs = sharp(b, g)
    for k in range(g + 1):                         # C = -b^#: powers 0..g
        Z[k + 1, 1, 0] = -bs[k]
    return Z, a


def generator(Z, kap, k, c):
    """G_{k,c} = c W_k + R(c W_k) as Laurent coeffs with powers -1-k .. 1+k."""
    W = Z[: k + 2] / kap                           # powers -1-k .. 0
    M = W[-1]
    W = W.copy()
    W[-1] = np.array([[M[0, 0] / 2, 0], [M[1, 0], M[1, 1] / 2]])
    lo = -1 - k
    G = np.zeros((2 * (k + 2) - 1, 2, 2), complex)
    for i in range(k + 2):
        p = lo + i
        G[i] += c * W[i]
        G[-p - lo] += -np.conj(c) * W[i].conj().T
    return G, lo


def flow_rhs(Z, kap, k, c):
    """d zeta / dt = [zeta, G_{k,c}]; result truncated to powers -1..g (the rest is 0)."""
    G, lo = generator(Z, kap, k, c)
    n = len(Z)
    R = np.zeros((n + len(G) - 1, 2, 2), complex)
    for i in range(n):
        for j in range(len(G)):
            R[i + j] += Z[i] @ G[j] - G[j] @ Z[i]
    return R[-lo: -lo + n]


def laurent_eval(X, lo, lam, deriv=False):
    if deriv:
        return sum((lo + i) * X[i] * lam ** (lo + i - 1) for i in range(len(X)))
    return sum(X[i] * lam ** (lo + i) for i in range(len(X)))


def rk4(f, y, h):
    k1 = f(y); k2 = f(y + h / 2 * k1); k3 = f(y + h / 2 * k2); k4 = f(y + h * k3)
    return y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def isospectral_basis(g):
    """Independent real flows (k, c) spanning the isospectral torus."""
    B = [(k, c) for k in range(g // 2) for c in (1, 1j)]
    if g % 2 == 1:
        B.append(((g - 1) // 2, 1j))
    return B


def move_on_isospectral_torus(Z, kap, times, h=0.01):
    g = len(Z) - 2
    for (k, c), t in zip(isospectral_basis(g), times):
        n = max(1, int(np.ceil(abs(t) / h)))
        for _ in range(n):
            Z = rk4(lambda Y: flow_rhs(Y, kap, k, c), Z, t / n)
    return Z


def su2_to_r3(X):
    return np.array([-0.5 * np.trace(X @ e).real for e in (E1, E2, E3)])


def surface(Z0, a, lam0, xs, ys, substeps=4):
    """Integrate zeta, F, dF/dlam over the grid xs x ys (xs[0]=ys[0]=0). Returns f, N, u."""
    kap = 2 * np.sqrt(abs(a[0]))
    ns = Z0.size

    def rhs(c):
        def f(y):
            Z = y[:ns].reshape(Z0.shape)
            F = y[ns:ns + 4].reshape(2, 2)
            G = y[ns + 4:].reshape(2, 2)
            W, lo = generator(Z, kap, 0, c)
            Wl = laurent_eval(W, lo, lam0)
            dWl = laurent_eval(W, lo, lam0, deriv=True)
            return np.concatenate([flow_rhs(Z, kap, 0, c).ravel(), (F @ Wl).ravel(),
                                   (G @ Wl + F @ dWl).ravel()])
        return f

    def step(y, c, h):
        for _ in range(substeps):
            y = rk4(rhs(c), y, h / substeps)
        return y

    nx, ny = len(xs), len(ys)
    f = np.zeros((nx, ny, 3)); N = np.zeros((nx, ny, 3)); u = np.zeros((nx, ny))
    y = np.concatenate([Z0.ravel(), np.eye(2, dtype=complex).ravel(), np.zeros(4, complex)])
    for i in range(nx):
        if i > 0:
            y = step(y, 1, xs[i] - xs[i - 1])
        yy = y
        for j in range(ny):
            if j > 0:
                yy = step(yy, 1j, ys[j] - ys[j - 1])
            Z = yy[:ns].reshape(Z0.shape)
            F = yy[ns:ns + 4].reshape(2, 2)
            G = yy[ns + 4:].reshape(2, 2)
            Fi = np.linalg.inv(F)
            f[i, j] = su2_to_r3(4j * lam0 * G @ Fi)
            N[i, j] = su2_to_r3(-F @ E3 @ Fi)
            u[i, j] = np.log(2 * abs(Z[0, 0, 1]) / kap)
    return f, N, u


def divisor(Z):
    """Points (mu_j, nu_j) on nu^2 = lam a(lam) where the eigenline of zeta is [0:1]."""
    g = len(Z) - 2
    if g == 0:
        return []
    b = Z[:, 0, 1][: g + 1]                        # b(lam) = lam*B(lam)
    A = Z[1: g + 1, 0, 0]
    mus = np.roots(b[::-1])
    return [(m, -m * np.polyval(A[::-1], m)) for m in mus]


# ----------------------------------------------------------------------------- tests
def _geom_check(alphas, A, lam0, times=None, L=1.0, n=41):
    Z, a = killing_field(alphas, A)
    kap = 2 * np.sqrt(abs(a[0]))
    if times is not None:
        Z = move_on_isospectral_torus(Z, kap, times)
    xs = np.linspace(0, L, n); h = xs[1]
    P, N, u = surface(Z, a, lam0, xs, xs)
    Px = (P[2:, 1:-1] - P[:-2, 1:-1]) / (2 * h); Py = (P[1:-1, 2:] - P[1:-1, :-2]) / (2 * h)
    Pxx = (P[2:, 1:-1] - 2 * P[1:-1, 1:-1] + P[:-2, 1:-1]) / h ** 2
    Pyy = (P[1:-1, 2:] - 2 * P[1:-1, 1:-1] + P[1:-1, :-2]) / h ** 2
    Pxy = (P[2:, 2:] - P[2:, :-2] - P[:-2, 2:] + P[:-2, :-2]) / (4 * h * h)
    E = (Px * Px).sum(-1); G = (Py * Py).sum(-1); F = (Px * Py).sum(-1)
    Nn = N[1:-1, 1:-1]
    H = ((Pxx + Pyy) * Nn).sum(-1) / (2 * E)
    Q = (0.25 * (Pxx - Pyy - 2j * Pxy) * Nn).sum(-1)
    uu = u[1:-1, 1:-1]
    print(f"  g={len(alphas)}: |E-G|/E={np.abs(E-G).max()/E.mean():.1e}  F/E={np.abs(F).max()/E.mean():.1e}"
          f"  E/(4e^2u)-1={np.abs(E/(4*np.exp(2*uu))-1).max():.1e}  H in [{H.min():.5f},{H.max():.5f}]"
          f"  Q err={np.abs(Q + np.exp(1j*np.angle(a[0]))/lam0).max():.1e}")


def _delaunay_check(r, sign):
    al = r * np.exp(0.9j)
    lam0 = sign * al / abs(al)
    Z, a = killing_field([al])
    ys = np.linspace(0, 12, 400)
    P, _, _ = surface(Z, a, lam0, np.linspace(0, 6, 15), ys[:1])   # x-curve (circle)
    X = P[:, 0]
    Acirc = np.c_[2 * X, np.ones(len(X))]
    s = np.linalg.lstsq(Acirc, (X ** 2).sum(1), rcond=None)[0]
    Y = X - X.mean(0); nrm = np.linalg.eigh(Y.T @ Y)[1][:, 0]
    cc = s[:3] - ((s[:3] - X.mean(0)) @ nrm) * nrm
    P2, _, _ = surface(Z, a, lam0, ys[:1], ys)                      # meridian x=0
    prof = P2[0] - cc
    ht = prof @ nrm
    rad = np.linalg.norm(prof - np.outer(ht, nrm), axis=1)
    exp = (2 * r / (1 + r), 2 / (1 + r)) if sign < 0 else (2 * r / (1 - r), 2 / (1 - r))
    print(f"  |alpha|={r}, lam0={'-' if sign<0 else '+'}alpha/|alpha|: radii [{rad.min():.4f},{rad.max():.4f}]"
          f"  expected [{exp[0]:.4f},{exp[1]:.4f}]")


if __name__ == "__main__":
    print("geometry checks (finite-difference errors are O(h^2)):")
    _geom_check([], None, np.exp(0.3j))
    _geom_check([0.4 + 0.2j], None, np.exp(0.7j))
    _geom_check([0.5 + 0.1j, -0.3 + 0.6j], [0.03 + 0.02j, -0.03 + 0.02j], np.exp(1.3j))
    _geom_check([0.3, 0.5j, -0.6 - 0.2j], None, 1j, times=[0.2, -0.3, 0.4])
    print("Delaunay checks (g=1):")
    for r in (0.3, 0.7):
        for sgn in (-1, 1):
            _delaunay_check(r, sgn)
    print("divisor at base point = branch points:",
          np.round(sorted([d[0] for d in divisor(killing_field([0.3, 0.5j])[0])], key=np.angle), 6))
