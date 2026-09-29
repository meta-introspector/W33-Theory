#!/usr/bin/env python3
"""Pass 11169: the perfect three-qutrit Clifford gates fall into exactly six local orbits, labelled by a permutation pi in
S3 -- which future receives each past ORIENTATION-REVERSED -- and into a single orbit once parties may be relabelled.

For a perfect S in Sp(6,3) every column of block determinants is (1, 1, -1) in some order (Pass 11158), and so is every
row (row law); so the positions of the -1's form a permutation matrix: det S_ij = -1 iff i = pi(j).  A block of
determinant -1 is an orientation-reversing (anti-symplectic) map of the qutrit phase space -- a local time reversal
(Pass 11156: the anti-unitary law) -- so pi says in which future each past appears time-reversed.
Local Clifford operations S -> L S R (L, R in SL(2,3)^3) preserve every block determinant, so pi is a local invariant.
Result: for a representative of each of the six classes, the stabiliser {R : S R S^{-1} local} has order 4, so its orbit
under SL(2,3)^6 (order 24^6 = 191 102 976) has 47 775 744 elements = exactly 1/6 of the 286 654 464 perfect gates
(Pass 11163).  Hence each pi-class is ONE local orbit: pi is a complete local invariant.  Relabelling output and input
parties (pi -> sigma pi tau^{-1}) joins the six, so up to local operations and relabelling the perfect three-qutrit gate is
unique -- as the perfect two-qutrit gate is (Pass 11156, stabiliser order 24 recomputed here).
The F9-linear perfect gates (Pass 11164; block det = norm of the entry) are spread evenly: 2048 in each class.
"""
from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11165_three_qutrit_arrow as A3  # noqa: E402

OUT = ROOT / "data" / "w33_pass11169_perfect_gate_local_orbits.json"


def sl23():
    return [np.array([[a, b], [c, d]]) for a, b, c, d in itertools.product(range(3), repeat=4) if (a * d - b * c) % 3 == 1]


def dets(S, n):
    return [[int(S[2 * i, 2 * j] * S[2 * i + 1, 2 * j + 1] - S[2 * i, 2 * j + 1] * S[2 * i + 1, 2 * j]) % 3 for j in range(n)] for i in range(n)]


def pattern(S, n=3):
    D = dets(S, n)
    if any(D[i][j] == 0 for i in range(n) for j in range(n)):
        return None
    return tuple(next(i for i in range(n) if D[i][j] == 2) for j in range(n))      # pi(j) = row of the -1 in column j


def inv_mod3(S):
    """inverse of a symplectic matrix: S^{-1} = -J S^T J"""
    n = S.shape[0] // 2
    J = np.zeros((2 * n, 2 * n), int)
    for k in range(n):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, -1
    return (-J @ S.T @ J) % 3


def local_batch(n):
    """all R in SL(2,3)^n as a (24^n, 2n, 2n) array"""
    G = np.array(sl23())
    out = []
    for combo in itertools.product(range(24), repeat=n):
        R = np.zeros((2 * n, 2 * n), int)
        for k, g in enumerate(combo):
            R[2 * k:2 * k + 2, 2 * k:2 * k + 2] = G[g]
        out.append(R)
    return np.array(out)


def stabiliser_order(S, Rs):
    n = S.shape[0] // 2
    Si = inv_mod3(S)
    T = np.einsum('ij,rjk,kl->ril', S, Rs, Si) % 3
    mask = np.ones((2 * n, 2 * n), bool)
    for k in range(n):
        mask[2 * k:2 * k + 2, 2 * k:2 * k + 2] = False
    return int((T[:, mask] == 0).all(axis=1).sum())


def summarize():
    rng = np.random.default_rng(11169)
    reps = {}
    while len(reps) < 6:
        M, _ = A3.random_symplectic(rng)
        p = pattern(M)
        if p is not None and p not in reps:
            reps[p] = M
    R3 = local_batch(3)
    stabs = {str(p): stabiliser_order(S, R3) for p, S in sorted(reps.items())}
    orbit = {p: 24 ** 6 // s for p, s in stabs.items()}
    # two-qutrit comparison: a perfect gate (both determinants -1)
    R2 = local_batch(2)
    two = None
    while two is None:
        M = np.eye(4, dtype=int)
        J = np.zeros((4, 4), int)
        J[0, 1], J[1, 0], J[2, 3], J[3, 2] = 1, -1, 1, -1
        for _ in range(40):
            v = rng.integers(0, 3, 4)
            M = (M + rng.integers(1, 3) * np.outer(v, (J @ v) @ M)) % 3
        D = dets(M, 2)
        if all(D[i][j] == 2 for i in range(2) for j in range(2)):
            two = M
    stab2 = stabiliser_order(two, R2)
    # F9-linear gates: enumerate U(3,F9) perfect ones directly and record their pi
    pats = Counter()
    nz = [(a, b) for a in range(3) for b in range(3) if (a, b) != (0, 0)]
    norm = {x: (x[0] ** 2 + x[1] ** 2) % 3 for x in nz}

    def herm(u, v):
        re = sum(x[0] * y[0] + x[1] * y[1] for x, y in zip(u, v)) % 3
        im = sum(x[0] * y[1] - x[1] * y[0] for x, y in zip(u, v)) % 3
        return re, im
    rows = [r for r in itertools.product(nz, repeat=3) if sum(norm[x] for x in r) % 3 == 1]
    for r1 in rows:
        for r2 in rows:
            if herm(r1, r2) != (0, 0):
                continue
            for r3 in rows:
                if herm(r1, r3) != (0, 0) or herm(r2, r3) != (0, 0):
                    continue
                K = [r1, r2, r3]
                ok = True
                for (i, j) in itertools.combinations(range(3), 2):
                    for (k, l) in itertools.combinations(range(3), 2):
                        a, b, c, e = K[i][k], K[i][l], K[j][k], K[j][l]
                        m = ((a[0] * e[0] - a[1] * e[1] - b[0] * c[0] + b[1] * c[1]) % 3,
                             (a[0] * e[1] + a[1] * e[0] - b[0] * c[1] - b[1] * c[0]) % 3)
                        ok &= m != (0, 0)
                if ok:
                    pats[tuple(next(i for i in range(3) if norm[K[i][j]] == 2) for j in range(3))] += 1
    total = 286654464
    res = dict(pass_id=11169, stabiliser_orders=stabs, orbit_sizes=orbit, class_size=total // 6,
               each_class_one_orbit=all(o == total // 6 for o in orbit.values()), n_classes=len(reps),
               two_qutrit_stabiliser=stab2, two_qutrit_orbit=24 ** 4 // stab2,
               f9_pattern_counts={str(k): v for k, v in sorted(pats.items())}, f9_total=sum(pats.values()))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
