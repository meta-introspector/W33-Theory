#!/usr/bin/env python3
"""Pass 11157: why ninth roots of unity -- one operator identity explains both cubics of the Pauli-only CGLMP Bell operator.

Pass 11152 found char(B) = (x^3 - 6x - 2)(x^3 - 3x + 1)^2 for the Bell operator B of the optimal Pauli strategy (Z and X
bases on both qutrits) and left the ninth roots of unity unexplained.  Here:
  * WEYL FORM: B is a sum of eight two-qutrit Weyl operators, Z(x)Z^2, Z(x)X^2, X(x)Z^2, X(x)X^2 and their adjoints, each
    with coefficient of modulus 1/sqrt3 and phase +-pi/6 or +-5pi/6;
  * SYMMETRY: the only nontrivial Pauli symmetry is U = (X Z^2)(x)(X Z^2), a joint relabelling of both parties' outcomes
    (U^3 = I); B splits into three 3-dimensional charge sectors U = 1, w, w^2;
  * THE IDENTITY (checked at 40 digits, exact because every entry lies in (1/9)Z[w]):
        B^3 - 3B + 1 = (B + 1)(I + U + U^2).
    On the charged sectors (U = w, w^2) the right side vanishes: B^3 - 3B + 1 = 0, i.e. with B = 2cos(theta),
    2cos(3 theta) = -1 = 2cos(2pi/3) -- the Chebyshev trisection of the symmetry's angle 2pi/3: theta = 2pi m/9,
    m = 1, 2, 4.  That is where the ninth roots come from.
    On the neutral sector (U = 1): B^3 - 3B + 1 = 3(B + 1), i.e. B^3 - 6B - 2 = 0, whose largest root 2.6017 is the
    Pauli-only CGLMP maximum.
So both cubics are one identity read in different charge sectors of the relabelling symmetry.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11152_pauli_bell_charpoly as P  # noqa: E402

OUT = ROOT / "data" / "w33_pass11157_ninth_roots_explained.json"


def summarize():
    _, B = P.bell_charpoly(40)
    mp.mp.dps = 40
    w = mp.mpc(-0.5, mp.sqrt(3) / 2)
    X = mp.matrix([[1 if i == (j + 1) % 3 else 0 for j in range(3)] for i in range(3)])
    Z = mp.matrix([[w ** i if i == j else 0 for j in range(3)] for i in range(3)])
    kron = lambda A, Bm: mp.matrix([[A[i // 3, j // 3] * Bm[i % 3, j % 3] for j in range(9)] for i in range(9)])
    U = kron(X * Z * Z, X * Z * Z)
    I = mp.eye(9)
    dev = max(abs(x) for x in (B * B * B - 3 * B + I - (B + I) * (I + U + U * U)))
    comm = max(abs(x) for x in (U * B - B * U))
    # Weyl support (numerical)
    Bn = np.array(B.tolist(), dtype=complex)
    W3 = np.exp(2j * np.pi / 3)
    Xn, Zn = np.roll(np.eye(3), 1, axis=0), np.diag([1, W3, W3 * W3])
    pa = lambda u: np.linalg.matrix_power(Xn, u[0]) @ np.linalg.matrix_power(Zn, u[1])
    supp = {}
    for u in itertools.product(range(3), repeat=2):
        for v in itertools.product(range(3), repeat=2):
            c = np.trace(np.kron(pa(u), pa(v)).conj().T @ Bn) / 9
            if abs(c) > 1e-9:
                supp[f"{u}|{v}"] = [round(abs(c), 9), round(float(np.angle(c) / np.pi), 9)]
    sym = [f"{u}|{v}" for u in itertools.product(range(3), repeat=2) for v in itertools.product(range(3), repeat=2)
           if np.abs(np.kron(pa(u), pa(v)) @ Bn - Bn @ np.kron(pa(u), pa(v))).max() < 1e-9]
    res = dict(pass_id=11157, identity_deviation=float(dev), commutator=float(comm), weyl_support=supp,
               weyl_terms=len(supp), all_moduli_one_over_sqrt3=all(abs(m - 1 / np.sqrt(3)) < 1e-8 for m, _ in supp.values()),
               pauli_symmetries=sym, ninth_root_eigs=[float(2 * np.cos(2 * np.pi * m / 9)) for m in (1, 2, 4)])
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    print(json.dumps({k: v for k, v in r.items() if k != 'weyl_support'}, indent=1))
