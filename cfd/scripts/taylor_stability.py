"""Linear stability of circular Couette flow (outer cylinder fixed) to axisymmetric disturbances.

Chebyshev collocation in primitive variables (u, v, w, p), lengths scaled by the gap d and velocity
by the inner-wall speed, so Re = Omega1 R1 d / nu. Returns the critical Re and axial wavenumber.
"""
import json
import os
import sys

import numpy as np
from scipy.linalg import eig
from scipy.optimize import brentq, minimize_scalar


def cheb(n):
    x = np.cos(np.pi * np.arange(n + 1) / n)
    c = np.hstack([2, np.ones(n - 1), 2]) * (-1) ** np.arange(n + 1)
    X = np.tile(x, (n + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1 / c) / (dX + np.eye(n + 1))
    D -= np.diag(D.sum(axis=1))
    return D, x


def growth(re, k, eta, n=40):
    r1, r2 = eta / (1 - eta), 1 / (1 - eta)
    D, x = cheb(n)
    r = r1 + (x + 1) * (r2 - r1) / 2
    D = D * 2 / (r2 - r1)
    D2 = D @ D
    I = np.eye(n + 1)
    R = np.diag(r)
    Ri = np.diag(1 / r)
    Ri2 = np.diag(1 / r**2)

    # base flow V = a r + b/r with V(r1) = 1, V(r2) = 0
    a = -r1 / (r2**2 - r1**2)
    b = r1 * r2**2 / (r2**2 - r1**2)
    V = a * r + b / r
    VV = np.diag(V)
    Zeta = np.diag(np.full_like(r, 2 * a))  # dV/dr + V/r

    lap_v = (D2 + Ri @ D - Ri2 - k**2 * I) / re
    lap_w = (D2 + Ri @ D - k**2 * I) / re
    Z = np.zeros_like(I)
    N = n + 1
    A = np.block([
        [lap_v, 2 * Ri @ VV, Z, -D],
        [-Zeta, lap_v, Z, Z],
        [Z, Z, lap_w, -1j * k * I],
        [D + Ri, Z, 1j * k * I, Z],
    ]).astype(complex)
    B = np.block([
        [I, Z, Z, Z],
        [Z, I, Z, Z],
        [Z, Z, I, Z],
        [Z, Z, Z, Z],
    ]).astype(complex)
    for blk in range(3):
        for row in (blk * N, blk * N + n):
            A[row, :] = 0
            B[row, :] = 0
            A[row, row] = 1
    w = eig(A, B, right=False)
    w = w[np.isfinite(w)]
    w = w[np.abs(w) < 1e6]
    return float(np.max(w.real))


def critical_re(eta, n=40):
    def re_neutral(k):
        f = lambda re: growth(re, k, eta, n)
        return brentq(f, 20, 2000, xtol=1e-6)

    res = minimize_scalar(re_neutral, bracket=(2.0, 3.1, 4.5), tol=1e-6)
    return res.fun, res.x


if __name__ == "__main__":
    out = {}
    for eta, label in ((0.5, "check (literature 68.19, k=3.16)"), (0.3, "this mug: R1/R2 = 0.3")):
        rec, kc = critical_re(eta)
        out[f"eta={eta}"] = {"Re_c": rec, "k_c": kc, "wavelength_over_d": 2 * np.pi / kc}
        print(f"eta={eta}: Re_c={rec:.3f}, k_c={kc:.3f}, axial wavelength={2*np.pi/kc:.3f} d   [{label}]")
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "post"))
    json.dump(out, open(os.path.join(root, "taylor_linear_stability.json"), "w"), indent=2)
