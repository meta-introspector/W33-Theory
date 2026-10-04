"""Pass 11422: the exact one-qutrit reversible fraction at depth 5, and the 0.72 coincidence.

Pass 11312 gave the exact violating fractions 1/12, 13/32, 421/768, 1425/2048 for 1..4 cubic gates and sampled values
beyond; the reversible fraction decays by about 0.72 per gate.  Pass 11372 found rho_(8,8) = 0.72276 (root of
2592x^3 - 1332x^2 - 585x + 140) and showed that sup rho > 0.72, so the match is unexplained.

Decider (exact, Pass 11355 Theorem 1): U is reversible iff U^T = lambda C U C^dag for some Clifford C, i.e. iff
max_C |tr((C U C^dag)^dag U^T)| = 3.  The 216 Cliffords are checked in one batched contraction; the gap between 3 and
the next value is reported (no threshold is guessed).  Validated at depth 4 against the exact count 2,077,650.

Words: R_k T ... R_2 T C_1 T (24 coset representatives per later gate, Pass 11312) -- the same multiset of verdicts
as the uniform walk.  Depth 5: 216 * 24^4 = 71,663,616 words.  Deeper depths are sampled with the same decider.
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
import w33_pass11312_depth_law as DL  # noqa: E402

OUT = ROOT / "data" / "w33_pass11422_depth5_exact.json"
_S = {}


def _init():
    CF = np.array(P1.clifford1())
    _S["CF"], _S["CFd"] = CF, CF.conj().transpose(0, 2, 1)
    _S["R"] = np.array(DL.coset_reps(list(CF)))
    _S["T"] = P1.T1
    _S["CT"] = np.einsum('cij,jk->cik', CF, P1.T1)


def overlaps(U):
    """max over Cliffords of |tr((C U C^dag)^dag U^T)| for a batch U (W,3,3)"""
    CF, CFd = _S["CF"], _S["CFd"]
    X = np.einsum('cij,wjk,ckl->wcil', CF, U, CFd)                 # C U C^dag
    return np.abs(np.einsum('wcij,wji->wc', X.conj(), U)).max(axis=1)   # sum_ij conj(X)_ij (U^T)_ij


def _block(prefix):
    """all words P R_c T R_d T ... C T with the given prefix (indices into R), last factor over all 216 C T"""
    R, T, CT = _S["R"], _S["T"], _S["CT"]
    k_rest = _S["k"] - 1 - len(prefix)
    Pm = np.eye(3, dtype=complex)
    for r in prefix:
        Pm = Pm @ R[r] @ T
    rev = 0
    gap_lo, gap_hi = 3.0, 0.0
    for rest in itertools.product(range(24), repeat=k_rest):
        W = Pm.copy()
        for r in rest:
            W = W @ R[r] @ T
        U = np.einsum('ij,cjk->cik', W, CT)
        ov = overlaps(U)
        rv = ov > 2.999
        rev += int(rv.sum())
        if (~rv).any():
            gap_hi = max(gap_hi, float(ov[~rv].max()))
        if rv.any():
            gap_lo = min(gap_lo, float(ov[rv].min()))
    return rev, gap_lo, gap_hi


def exact(k, nproc=8):
    _S["k"] = k
    prefixes = [()] if k <= 2 else [tuple(p) for p in itertools.product(range(24), repeat=min(2, k - 2))]
    with Pool(nproc, initializer=_init_k, initargs=(k,)) as pool:
        out = pool.map(_block, prefixes, chunksize=1)
    words = 216 * 24 ** (k - 1)
    rev = sum(o[0] for o in out)
    return dict(words=words, reversible=rev, violating=words - rev,
                reversible_fraction=str(Fraction(rev, words)), violating_fraction=str(Fraction(words - rev, words)),
                min_overlap_reversible=min(o[1] for o in out), max_overlap_violating=max(o[2] for o in out))


def _init_k(k):
    _init()
    _S["k"] = k


def _sample(args):
    k, seed, n = args
    _init()
    CT, rng = _S["CT"], np.random.default_rng(seed)
    rev = 0
    for _ in range(n // 512):
        idx = rng.integers(216, size=(512, k))
        U = np.broadcast_to(np.eye(3, dtype=complex), (512, 3, 3)).copy()
        for j in range(k):
            U = np.einsum('wij,wjk->wik', CT[idx[:, j]], U)
        rev += int((overlaps(U) > 2.999).sum())
    return rev, (n // 512) * 512


def sampled(ks, per_k, nproc=8):
    res = {}
    with Pool(nproc) as pool:
        for k in ks:
            parts = pool.map(_sample, [(k, 114220000 + 1000 * k + i, per_k // nproc) for i in range(nproc)])
            rev, n = sum(p[0] for p in parts), sum(p[1] for p in parts)
            p = rev / n
            res[str(k)] = dict(samples=n, reversible=rev, reversible_fraction=p, stderr=float(np.sqrt(p * (1 - p) / n)))
            print(k, res[str(k)], flush=True)
    return res


def run():
    res = dict(pass_id=11422)
    known = {1: Fraction(11, 12), 2: Fraction(19, 32), 3: Fraction(347, 768), 4: Fraction(623, 2048)}
    ex = {}
    for k in (1, 2, 3, 4, 5):
        ex[str(k)] = exact(k)
        if k in known:
            assert Fraction(ex[str(k)]["reversible_fraction"]) == known[k], (k, ex[str(k)])
        print(k, ex[str(k)], flush=True)
    res["exact"] = ex
    fr = [Fraction(ex[str(k)]["reversible_fraction"]) for k in range(1, 6)]
    res["exact_ratios"] = [float(fr[i + 1] / fr[i]) for i in range(4)]
    res["sampled"] = sampled([6, 8, 10, 12, 14, 16, 20], 4_000_000)
    rho88 = float(max(np.roots([2592, -1332, -585, 140]).real))
    res["rho_8_8"] = rho88
    # geometric fit of the reversible fraction on k >= 8 (log-linear, weighted)
    ks = np.array([8, 10, 12, 14, 16, 20])
    p = np.array([res["sampled"][str(k)]["reversible_fraction"] for k in ks])
    se = np.array([res["sampled"][str(k)]["stderr"] for k in ks])
    w = (p / se) ** 2
    A = np.vstack([np.ones_like(ks, dtype=float), ks.astype(float)]).T
    coef, *_ = np.linalg.lstsq(A * np.sqrt(w)[:, None], np.log(p) * np.sqrt(w), rcond=None)
    cov = np.linalg.inv((A * w[:, None]).T @ A)
    rate, rate_se = float(np.exp(coef[1])), float(np.exp(coef[1]) * np.sqrt(cov[1, 1]))
    res["fitted_decay_rate_k8_20"] = dict(rate=rate, stderr=rate_se, z_vs_rho88=(rate - rho88) / rate_se)
    print(res["fitted_decay_rate_k8_20"], flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
