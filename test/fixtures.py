"""Dump reference values from cmc_reference.py to test/fixtures.json for the JS port.

Run from the repository root:  python3 test/fixtures.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import cmc_reference as ref  # noqa: E402


def cx(z):
    return [float(np.real(z)), float(np.imag(z))]


def flat(Z):
    """(n,2,2) complex -> flat [re, im, ...] in the JS slot layout."""
    return [float(v) for m in Z for z in m.reshape(4) for v in (z.real, z.imag)]


def point(Z0, a, lam0, x, y, h=0.005):
    """f, N, u at z = x + iy, integrating along x then y from F(0) = I."""
    kap = 2 * np.sqrt(abs(a[0]))
    ns = Z0.size

    def rhs(c):
        def f(v):
            Z = v[:ns].reshape(Z0.shape)
            F = v[ns:ns + 4].reshape(2, 2)
            G = v[ns + 4:].reshape(2, 2)
            W, lo = ref.generator(Z, kap, 0, c)
            Wl = ref.laurent_eval(W, lo, lam0)
            dWl = ref.laurent_eval(W, lo, lam0, deriv=True)
            return np.concatenate([ref.flow_rhs(Z, kap, 0, c).ravel(), (F @ Wl).ravel(),
                                   (G @ Wl + F @ dWl).ravel()])
        return f

    v = np.concatenate([Z0.ravel(), np.eye(2, dtype=complex).ravel(), np.zeros(4, complex)])
    for c, d in ((1, x), (1j, y)):
        n = max(1, int(np.ceil(abs(d) / h)))
        for _ in range(n):
            v = ref.rk4(rhs(c), v, d / n)
    F = v[ns:ns + 4].reshape(2, 2)
    G = v[ns + 4:].reshape(2, 2)
    Fi = np.linalg.inv(F)
    Z = v[:ns].reshape(Z0.shape)
    return {
        "z": [x, y],
        "f": ref.su2_to_r3(4j * lam0 * G @ Fi).tolist(),
        "N": ref.su2_to_r3(-F @ ref.E3 @ Fi).tolist(),
        "u": float(np.log(2 * abs(Z[0, 0, 1]) / kap)),
    }


CASES = [
    dict(alphas=[], A=None, theta0=0.3, times=None),
    dict(alphas=[0.4 + 0.2j], A=None, theta0=0.7, times=None),
    dict(alphas=[0.5 + 0.1j, -0.3 + 0.6j], A=[0.03 + 0.02j, -0.03 + 0.02j], theta0=1.3, times=None),
    dict(alphas=[0.3, 0.5j, -0.6 - 0.2j], A=None, theta0=np.pi / 2, times=[0.2, -0.3, 0.4]),
    dict(alphas=[0.6, 0.4j, -0.5 - 0.3j, 0.2 - 0.7j], A=None, theta0=2.0, times=[0.1, 0.2, -0.3, 0.25]),
]

out = []
for c in CASES:
    Zb, a = ref.killing_field(c["alphas"], c["A"])
    kap = 2 * np.sqrt(abs(a[0]))
    Zi = Zb if c["times"] is None else ref.move_on_isospectral_torus(Zb, kap, c["times"], h=0.002)
    lam0 = np.exp(1j * c["theta0"])
    pts = [point(Zi, a, lam0, x, y) for (x, y) in ((0.7, 0.0), (0.5, -0.8), (-1.2, 1.5))]
    out.append({
        "alphas": [cx(z) for z in c["alphas"]],
        "A": None if c["A"] is None else [cx(z) for z in c["A"]],
        "theta0": c["theta0"],
        "times": c["times"],
        "a": [v for z in a for v in cx(z)],
        "Zbase": flat(Zb),
        "Zinit": flat(Zi),
        "divisor": [cx(d[0]) for d in ref.divisor(Zi)],
        "points": pts,
    })
    print(f"g={len(c['alphas'])} done", flush=True)

path = os.path.join(os.path.dirname(__file__), "fixtures.json")
with open(path, "w") as fh:
    json.dump(out, fh, indent=1)
print("wrote", path)
