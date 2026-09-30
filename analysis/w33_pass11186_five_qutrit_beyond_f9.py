#!/usr/bin/env python3
"""Pass 11186: perfect five-qutrit gates beyond F9 -- AME(10,3) graph states by tabu search, and their orientation
patterns.

Pass 11175: the determinant law allows five orientation patterns for perfect five-qutrit gates ('C10', 'C4+C6',
'cross+perm', 'double cross', 'all'); the 2 642 411 520 F9-linear perfect gates realise only 'C10' and 'cross+perm', and a
plain hill climb found no AME(10,3) graph state.  Here: tabu search (analysis/w33_pass11186_scan_ame10_tabu.py; frozen
data/w33_pass11186_ame10_tabu.json) for weighted qutrit graphs on 10 vertices with all 126 balanced blocks invertible;
each graph found gives, for every choice of five input legs, a perfect gate S = M T (Pass 11175's conversion) whose
pattern is recorded.
RESULTS.  Tabu search found 71 AME(10,3) graph states in 300 restarts (the other 229 stalled at 8 singular blocks); the
71 x 252 = 17 892 perfect gates they give are symplectic and perfect, and realise EXACTLY THE SAME TWO patterns as the F9
census: 'C10' (5112) and 'cross+perm' (12 780).  No 'C4+C6', 'double cross' or 'all' in either source -- conjecture: only
these two orientation patterns occur for perfect five-qutrit gates.
THE DETERMINANT LAW OF EVERY ORDER (verified here on random symplectic matrices, n = 2..5).  The symplectic form is a sum
of local forms, omega = sum_i omega_i, and S preserves omega^k / k! = sum_{|I| = k} wedge_{i in I} omega_i; evaluating on
the 2k-dimensional input space of a party set J gives
        sum_{|I| = k} det S[I, J] = 1  (mod 3)   for every k and every set J of k input parties
(k = 1 is Pass 11158's law, k = n is det S = 1).  Jacobi's identity with S^-1 = -J S^T J gives the complement symmetry
det S[I, J] = det S[I^c, J^c] exactly.  These laws constrain the higher-order patterns but do not by themselves exclude
the three unrealised first-order patterns; that exclusion remains open.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11175_five_qutrit_patterns as P  # noqa: E402

SCAN = ROOT / "data" / "w33_pass11186_ame10_tabu.json"
OUT = ROOT / "data" / "w33_pass11186_five_qutrit_beyond_f9.json"


def determinant_laws(seed=1, reps=6):
    import itertools
    import numpy as np

    def rsym(n, rng, steps=120):
        D = 2 * n
        J = np.zeros((D, D), np.int64)
        for k in range(n):
            J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, -1
        M = np.eye(D, dtype=np.int64)
        for _ in range(steps):
            v = rng.integers(0, 3, D)
            c = int(rng.integers(1, 3))
            M = (M + c * np.outer(v, (v @ J) @ M)) % 3
        return M

    def minor(S, I, J):
        r = [2 * i + t for i in I for t in range(2)]
        c = [2 * j + t for j in J for t in range(2)]
        return int(round(np.linalg.det(S[np.ix_(r, c)].astype(float)))) % 3
    rng = np.random.default_rng(seed)
    law = jac = True
    for n in (2, 3, 4, 5):
        for _ in range(reps):
            S = rsym(n, rng)
            for k in range(1, n + 1):
                for Js in itertools.combinations(range(n), k):
                    law &= sum(minor(S, I, Js) for I in itertools.combinations(range(n), k)) % 3 == 1
                    if k < n:
                        for I in itertools.combinations(range(n), k):
                            Ic = [x for x in range(n) if x not in I]
                            Jc = [x for x in range(n) if x not in Js]
                            jac &= minor(S, I, Js) == minor(S, Ic, Jc)
    return bool(law), bool(jac)


def summarize():
    d = json.loads(SCAN.read_text())
    named = Counter()
    for k, v in d['gate_patterns'].items():
        named[P.name(k)] += v
    f9 = {'C10', 'cross+perm'}
    law, jac = determinant_laws()
    res = dict(pass_id=11186, determinant_law_every_order=law, jacobi_complement_symmetry=jac, restarts=d['restarts'], best_distribution=d['distribution'], graphs_found=len(d['graphs']),
               gates_all_perfect_symplectic=d['gates_all_perfect_symplectic'], patterns=dict(named),
               new_beyond_f9=sorted(set(named) - f9), all_five_realised=set(named) | f9 == {'C10', 'C4+C6', 'cross+perm', 'double cross', 'all'})
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
