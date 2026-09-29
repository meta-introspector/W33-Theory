#!/usr/bin/env python3
"""Pass 11162: psi* = sqrt(7/15)|000> + sqrt(2/15) sum_{a=1,2} |a>(|0a> + |a0>) is a strict, nondegenerate local maximum of
N_AB + N_AC = 4/sqrt15 modulo local unitaries, and 96 further random restarts find nothing higher.

(1) Smoothness.  At psi* the partial transposes of rho_AB and rho_AC have spectrum
        {-0.1915 (x2), -2/15, 2/15 (x3), 0.3249 (x2), 7/15},
    no eigenvalue within 2/15 of zero, so the negativity (minus the sum of the negative eigenvalues) is real-analytic in a
    neighbourhood of psi* and a second-order test is legitimate.
(2) Hessian.  On the 52-dimensional real tangent space of the unit sphere modulo phase at psi*, the gradient vanishes
    (|grad| ~ 3e-9, finite differences) and the Hessian has 32 strictly negative eigenvalues (-4.58 ... -0.72) and exactly
    20 zero eigenvalues -- the dimension of the local-unitary orbit through psi* (3 x su(3) = 24 minus a 4-dimensional
    stabiliser), computed independently from the orbit's tangent vectors.  So psi* is a Morse-Bott maximum: strict modulo
    local unitaries, with a margin 0.72 against finite-difference errors ~1e-7.
(3) Global evidence.  96 restarts (analysis/w33_pass11162_scan_restarts.py; random scales 0.2-3) all converge to
    4/sqrt15 to 1e-15; none exceeds it (plus 16 restarts in Pass 11148 and the Pass 11154 scan).
Scope: a certified local maximum and strong numerical evidence for the global one -- NOT a proof of global optimality.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESTARTS = ROOT / "data" / "w33_pass11162_restarts.json"
OUT = ROOT / "data" / "w33_pass11162_polygamy_local_maximum.json"


def pair(v, keep):
    T = v.reshape(3, 3, 3)
    M = T.reshape(9, 3) if keep == (0, 1) else T.transpose(0, 2, 1).reshape(9, 3)
    return M @ M.conj().T


def pt_spectrum(rho):
    pt = rho.reshape(3, 3, 3, 3).transpose(0, 3, 2, 1).reshape(9, 9)
    return np.linalg.eigvalsh((pt + pt.conj().T) / 2)


def f(v):
    v = v / np.linalg.norm(v)
    return sum(float(-e[e < 0].sum()) for e in (pt_spectrum(pair(v, (0, 1))), pt_spectrum(pair(v, (0, 2)))))


def psi_star():
    x = 7 / 15
    psi = np.zeros(27, complex)
    psi[0] = np.sqrt(x)
    b = np.sqrt((1 - x) / 4)
    for a in (1, 2):
        psi[a * 9 + a] = b
        psi[a * 9 + a * 3] = b
    return psi


def tangent_basis(psi):
    basis = []
    for k in range(27):
        for c in (1, 1j):
            e = np.zeros(27, complex)
            e[k] = c
            basis.append(e - (psi.conj() @ e) * psi)
    Bm = np.array(basis).T
    R = np.vstack([Bm.real, Bm.imag])
    Q, s, _ = np.linalg.svd(R, full_matrices=False)
    Q = Q[:, s > 1e-9]
    return Q[:27] + 1j * Q[27:]


def lu_orbit_dim(psi):
    gm = []
    for a in range(3):
        for b in range(3):
            E = np.zeros((3, 3), complex)
            E[a, b] = 1
            gm += [E + E.conj().T, 1j * (E - E.conj().T)]
    tang = []
    for site in range(3):
        for M in gm:
            ops = [np.eye(3)] * 3
            ops[site] = M
            t = 1j * np.kron(np.kron(ops[0], ops[1]), ops[2]) @ psi
            t = t - (psi.conj() @ t) * psi
            tang.append(np.concatenate([t.real, t.imag]))
    return int(np.linalg.matrix_rank(np.array(tang), tol=1e-8))


def summarize():
    psi = psi_star()
    f0 = f(psi)
    spec = pt_spectrum(pair(psi, (0, 1)))
    T = tangent_basis(psi)
    n, h = T.shape[1], 1e-4
    H = np.zeros((n, n))
    for i in range(n):
        H[i, i] = (f(psi + h * T[:, i]) - 2 * f0 + f(psi - h * T[:, i])) / h ** 2
        for j in range(i + 1, n):
            H[i, j] = H[j, i] = (f(psi + h * (T[:, i] + T[:, j])) - f(psi + h * (T[:, i] - T[:, j]))
                                 - f(psi - h * (T[:, i] - T[:, j])) + f(psi - h * (T[:, i] + T[:, j]))) / (4 * h * h)
    grad = np.array([(f(psi + h * T[:, i]) - f(psi - h * T[:, i])) / (2 * h) for i in range(n)])
    ev = np.linalg.eigvalsh(H)
    rs = json.loads(RESTARTS.read_text())
    res = dict(pass_id=11162, f_at_psi=f0, four_over_sqrt15=4 / np.sqrt(15), pt_spectrum=[round(float(e), 6) for e in spec],
               pt_gap_from_zero=float(np.abs(spec).min()), tangent_dim=n, grad_max=float(np.abs(grad).max()),
               hessian_negative=int((ev < -1e-3).sum()), hessian_zero=int((np.abs(ev) < 1e-3).sum()),
               hessian_positive=int((ev > 1e-3).sum()), hessian_most_negative=float(ev.min()),
               hessian_least_negative=float(ev[ev < -1e-3].max()), lu_orbit_dim=lu_orbit_dim(psi),
               restarts=len(rs['values']), restarts_best=rs['best'],
               restarts_above_target=sum(v > rs['target'] + 1e-9 for v in rs['values']),
               restarts_at_target=sum(abs(v - rs['target']) < 1e-9 for v in rs['values']))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
