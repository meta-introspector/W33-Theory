"""Pass 11312: one-qutrit T-violation versus the number of cubic gates -- exact through four gates.

Pass 11266 (exhaustive): a uniformly random one-qutrit Clifford word C_k T ... C_1 T violates substrate time reversal
with probability 1/12, 13/32, 421/768 for k = 1, 2, 3.

REDUCTION.  The nine diagonal Cliffords D = Z^a S^b (mod phase) commute with T.  Writing each C_i = R_i D_i with R_i one
of the 24 right-coset representatives of Cl/Diag, C_i T = R_i T D_i, and D_i is absorbed into C_{i-1} (still uniform).
So the uniform distribution over 216^k words equals the uniform distribution over words R_k T ... R_2 T C_1 T with
R_i over 24 coset representatives and C_1 over all 216: 216 * 24^(k-1) words, each of equal weight.  Checked against
the exhaustive k = 1, 2, 3 values; k = 4 (2,985,984 words) is new.  k = 5..8 by uniform sampling.
"""

from __future__ import annotations

import itertools
import json
import sys
from fractions import Fraction
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11213_cubic_t_violation as P1  # noqa: E402
import w33_pass11252_exact_reversibility as R  # noqa: E402

OUT = ROOT / "data" / "w33_pass11312_depth_law.json"
_S = {}


def coset_reps(Cf):
    """24 right-coset representatives of Cl / Diag (C = R D, D diagonal)"""
    diag = [C for C in Cf if np.allclose(C, np.diag(np.diag(C)))]
    assert len(diag) == 9
    reps, seen = [], []
    for C in Cf:
        if any(P1.key(C @ np.linalg.inv(D)) in seen for D in diag):
            continue
        reps.append(C)
        seen.extend(P1.key(C @ np.linalg.inv(D)) for D in diag)
    assert len(reps) == 24
    return reps


def _init():
    R.WEYL[1] = R.Weyl(1)
    Cf = P1.clifford1()
    _S["Cf"] = Cf
    _S["R"] = coset_reps(Cf)


def _block(args):
    k, prefix = args
    Cf, Rr, T = _S["Cf"], _S["R"], P1.T1
    rng = np.random.default_rng(0)
    v = 0
    for rest in itertools.product(range(24), repeat=k - 1 - len(prefix)):
        word = list(prefix) + list(rest)                     # R_k ... R_2 (indices)
        W = np.eye(3, dtype=complex)
        for r in word:
            W = W @ Rr[r] @ T
        for C in Cf:                                         # C_1
            if R.decide(W @ C @ T, 1, rng)[0] is False:
                v += 1
    return v


def _sample(args):
    k, seed, n = args
    Cf, T = _S["Cf"], P1.T1
    rng = np.random.default_rng(seed)
    v = 0
    for _ in range(n):
        U = np.eye(3, dtype=complex)
        for _ in range(k):
            U = Cf[rng.integers(216)] @ T @ U
        if R.decide(U, 1, rng)[0] is False:
            v += 1
    return v


def run():
    res = dict(pass_id=11312, exact={}, sampled={})
    with Pool(11, initializer=_init) as pool:
        for k in (1, 2, 3, 4):
            if k == 1:
                jobs = [(1, ())]
            else:
                jobs = [(k, (i,)) for i in range(24)] if k <= 3 else [(k, (i, j)) for i in range(24) for j in range(24)]
            v = sum(pool.map(_block, jobs))
            tot = 216 * 24 ** (k - 1)
            res["exact"][str(k)] = dict(violating=v, words=tot, probability=str(Fraction(v, tot)), float=v / tot)
            print(k, res["exact"][str(k)], flush=True)
        for k in (5, 6, 8, 12):
            n = 66000
            v = sum(pool.map(_sample, [(k, 100 * k + j, n // 11) for j in range(11)]))
            p = v / n
            res["sampled"][str(k)] = dict(violating=v, samples=n, probability=p, stderr=float(np.sqrt(p * (1 - p) / n)))
            print(k, res["sampled"][str(k)], flush=True)
    ex = res["exact"]
    res["matches_pass_11266"] = (ex["1"]["probability"] == "1/12" and ex["2"]["probability"] == "13/32"
                                 and ex["3"]["probability"] == "421/768")
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
