"""Pass 11204: a second law from counting -- arrow-free reversible dynamics becomes rare as qutrits are added.

A Clifford tick S of n qutrits has intrinsic arrow A(S) = 0 (Pass 11183) iff in some stabilizer tensor factorisation it
exports nothing, i.e. iff F3^{2n} is an orthogonal direct sum of n S-invariant nondegenerate planes (S is then a product
of independent one-qutrit ticks in that split).  Exact fractions (Pass 11188): n = 1: 1; n = 2: 133/360; n = 3:
55241/884520.  This pass tests the decomposition exactly for uniformly random S in Sp(2n,3) (random symplectic
bases) and estimates the arrow-free fraction for n = 4, 5, validating the sampler on n = 2, 3 against the exact values.

Test.  Invariant planes are span(v, Sv) with S^2 v in span(v, Sv) (v not an eigenvector), or spans of two eigenvectors
with the same eigenvalue; keep the nondegenerate ones and search (backtracking) for n of them that are pairwise
orthogonal -- they then span F3^{2n}.
"""

from __future__ import annotations

import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11204_second_law_counting.json"


def form(n):
    J = np.zeros((2 * n, 2 * n), np.int64)
    for k in range(n):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, -1
    return J


def random_sp(n, rng):
    """uniform element of Sp(2n,3): columns = a random symplectic basis"""
    J = form(n)
    D = 2 * n
    cols = []
    while len(cols) < D:
        basis = np.array(cols).T if cols else np.zeros((D, 0), np.int64)
        while True:
            u = rng.integers(0, 3, D)
            if not u.any():
                continue
            if basis.shape[1] and ((u @ J @ basis) % 3).any():
                continue
            break
        while True:
            v = rng.integers(0, 3, D)
            if basis.shape[1] and ((v @ J @ basis) % 3).any():
                continue
            w = (u @ J @ v) % 3
            if w:
                v = (v * pow(int(w), -1, 3)) % 3
                break
        cols += [u, v]
    S = np.array(cols).T % 3
    assert np.array_equal((S.T @ J @ S) % 3, J % 3)
    return S


def all_vectors(n):
    V = np.array(list(itertools.product(range(3), repeat=2 * n)), np.int64)[1:]
    first = np.argmax(V != 0, axis=1)
    lead = V[np.arange(len(V)), first]
    return V[lead == 1]                      # one representative per projective point


def rank_le2(a, b, c):
    """rows a, b, c (m x D) span a space of dim <= 2 (all 3x3 minors of [a b c] vanish), vectorised"""
    D = a.shape[1]
    ok = np.ones(len(a), bool)
    for i, j, k in itertools.combinations(range(D), 3):
        M = np.stack([a[:, [i, j, k]], b[:, [i, j, k]], c[:, [i, j, k]]], 1)
        det = (M[:, 0, 0] * (M[:, 1, 1] * M[:, 2, 2] - M[:, 1, 2] * M[:, 2, 1])
               - M[:, 0, 1] * (M[:, 1, 0] * M[:, 2, 2] - M[:, 1, 2] * M[:, 2, 0])
               + M[:, 0, 2] * (M[:, 1, 0] * M[:, 2, 1] - M[:, 1, 1] * M[:, 2, 0])) % 3
        ok &= det == 0
        if not ok.any():
            break
    return ok


def invariant_planes(S, V, J):
    SV = (V @ S.T) % 3
    S2V = (SV @ S.T) % 3
    eig = {}
    for lam in (1, 2):
        eig[lam] = V[((SV - lam * V) % 3 == 0).all(1)]
    is_eig = ((SV - V) % 3 == 0).all(1) | ((SV - 2 * V) % 3 == 0).all(1)
    cyc = ~is_eig & rank_le2(V, SV, S2V)
    planes = []
    for v, w in zip(V[cyc], SV[cyc]):
        if (v @ J @ w) % 3:
            planes.append((v, w))
    for lam in (1, 2):
        E = eig[lam]
        for a, b in itertools.combinations(range(len(E)), 2):
            if (E[a] @ J @ E[b]) % 3:
                planes.append((E[a], E[b]))
    return planes


def decomposes(planes, n, J):
    """is there a set of n pairwise orthogonal planes among the candidates?"""
    P = [np.stack(p) for p in planes]
    m = len(P)
    orth = np.zeros((m, m), bool)
    for i in range(m):
        for j in range(i + 1, m):
            orth[i, j] = orth[j, i] = not ((P[i] @ J @ P[j].T) % 3).any()

    def rec(chosen, cand):
        if len(chosen) == n:
            return True
        for k in cand:
            if rec(chosen + [k], [c for c in cand if c > k and orth[k, c]]):
                return True
        return False
    return rec([], list(range(m)))


def arrow_free_fraction(n, samples, seed):
    rng = np.random.default_rng(seed)
    V = all_vectors(n)
    J = form(n)
    hits = 0
    for _ in range(samples):
        S = random_sp(n, rng)
        hits += decomposes(invariant_planes(S, V, J), n, J)
    return hits


def run():
    res = dict(pass_id=11204, exact={"1": "1", "2": "133/360", "3": "55241/884520"})
    plan = {2: 20000, 3: 20000, 4: 20000, 5: 3000}
    est = {}
    for n, N in plan.items():
        h = arrow_free_fraction(n, N, 11204 + n)
        p = h / N
        est[n] = dict(samples=N, arrow_free=h, fraction=p, stderr=float(np.sqrt(max(p * (1 - p), 1e-12) / N)))
    res["estimates"] = {str(k): v for k, v in est.items()}
    for n, fr in ((2, Fraction(133, 360)), (3, Fraction(55241, 884520))):
        e = est[n]
        res[f"validation_n{n}"] = abs(e["fraction"] - float(fr)) < 4 * e["stderr"] + 1e-3
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
