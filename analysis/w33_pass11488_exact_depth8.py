"""Pass 11488: the exact reversible fraction at eight cubic gates, by weighted deduplication of operators.

The K-reduced walk of Pass 11436: words R_(i_(k-1)) T ... R_(i_1) T C T, C over the 216 Cliffords and R over the 8 coset
representatives of K (|K| = 27), uniform.  P_k = #(reversible words) / (216 * 8^(k-1)), reversibility decided by
Theorem 1 (Pass 11355) as max_C |tr((C U C^dag)^dag U^T)| = 3, cut 3 - 1e-7 (Pass 11459; gap reported).

Two words with the same operator (up to phase) have the same extensions, so the walk can be run on DISTINCT operators
carrying their word multiplicities.  Validated by reproducing the exact P_1..P_7 of Passes 11312/11422/11436; P_7 and P_8
are weighted counts over the 8 and 64 extensions of the depth-6 distinct set (no deduplication at depth 7 needed).

Deduplication keys: entries after fixing the phase of the first entry of modulus > 0.05, rounded at 1e-6; the count of
distinct operators is required to be identical at 1e-5 and 1e-9 roundings (a merge or split would change it).
"""

from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11422_depth5_exact as E  # noqa: E402
import w33_pass11436_depth_reduction as RD  # noqa: E402

OUT = ROOT / "data" / "w33_pass11488_exact_depth8.json"
TIGHT = 3 - 1e-7
KNOWN = {1: Fraction(11, 12), 2: Fraction(19, 32), 3: Fraction(347, 768), 4: Fraction(623, 2048),
         5: Fraction(11003, 49152), 6: Fraction(20327, 131072), 7: Fraction(353927, 3145728)}


def _canon(U, dec):
    f = U.reshape(len(U), 9)
    idx = np.argmax(np.abs(f) > 0.05, axis=1)
    ph = f[np.arange(len(f)), idx]
    g = np.round(f * (np.abs(ph) / ph)[:, None], dec) + 0.0
    return np.ascontiguousarray(g.view(np.float64).reshape(len(g), 18))


def dedupe(U, w):
    keys = _canon(U, 6)
    _, first, inv = np.unique(keys, axis=0, return_index=True, return_inverse=True)
    W = np.bincount(inv.ravel(), weights=w, minlength=len(first))
    n5 = np.unique(_canon(U, 5), axis=0).shape[0]
    n9 = np.unique(_canon(U, 9), axis=0).shape[0]
    assert n5 == n9 == len(first), (n5, n9, len(first))
    return U[first], W.astype(np.int64)


def _fast_overlaps(U):
    """max_C |tr((C U C^dag)^dag U^T)| = max_C |tr(C U^dag C^dag U^T)|, as two batched matmuls (4x faster than the
    three-operand einsum of Pass 11422; agreement to 1e-15 checked in run())"""
    CF, CFd = E._S["CF"], E._S["CFd"]
    T = np.matmul(np.matmul(CF[None], np.conj(np.transpose(U, (0, 2, 1)))[:, None]), CFd[None])
    return np.abs(np.einsum('wcik,wki->wc', T, np.transpose(U, (0, 2, 1)))).max(axis=1)


def overlaps(U, chunk=5000):
    return np.concatenate([_fast_overlaps(U[i:i + chunk]) for i in range(0, len(U), chunk)])


def _stream_job(args):
    """weighted reversible counts and overlap extremes over the 8 and 64 extensions of a block of distinct operators"""
    U, w = args
    if "CF" not in E._S:
        E._init()
        RD._init()
    RT = RD._S["RT"]
    X1 = np.einsum('rij,wjk->wrik', RT, U).reshape(-1, 3, 3)
    w1 = np.repeat(w, 8)
    X2 = np.einsum('rij,wjk->wrik', RT, X1).reshape(-1, 3, 3)
    w2 = np.repeat(w1, 8)
    out = []
    for X, ww in ((X1, w1), (X2, w2)):
        ov = overlaps(X)
        rev = ov > TIGHT
        out.append((int((ww * rev).sum()), float(ov[rev].min()) if rev.any() else 3.0,
                    float(ov[~rev].max()) if (~rev).any() else 0.0))
    return out


def run(kmax=8):
    RD._init()
    E._init()
    CT, RT = RD._S["CT"], RD._S["RT"]
    probe = np.einsum('rij,wjk->rwik', RT, CT).reshape(-1, 3, 3)
    assert np.abs(_fast_overlaps(probe) - E.overlaps(probe)).max() < 1e-12
    U = CT.copy()
    w = np.ones(len(U), dtype=np.int64)
    res = dict(pass_id=11488, levels={})
    gap = [3.0, 0.0]
    for k in range(1, kmax - 1):
        t0 = time.time()
        if k > 1:
            U = np.einsum('rij,wjk->wrik', RT, U).reshape(-1, 3, 3)
            w = np.repeat(w, 8)
            U, w = dedupe(U, w)
        ov = overlaps(U)
        rev = ov > TIGHT
        gap = [min(gap[0], float(ov[rev].min())), max(gap[1], float(ov[~rev].max()))]
        P = Fraction(int(w[rev].sum()), 216 * 8 ** (k - 1))
        assert int(w.sum()) == 216 * 8 ** (k - 1)
        res["levels"][k] = dict(distinct_operators=len(U), P=str(P), matches_known=(P == KNOWN.get(k)) if k in KNOWN else None,
                                seconds=round(time.time() - t0, 1))
        print(k, res["levels"][k], flush=True)
    # depths kmax-1 and kmax stream from the depth (kmax-2) distinct set: 8 and 64 extensions, no deduplication
    from multiprocessing import Pool
    t0 = time.time()
    num1 = num2 = 0
    B = 1000
    blocks = [(U[i:i + B], w[i:i + B]) for i in range(0, len(U), B)]
    with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 6) as pool:
        for n, (r1, r2) in enumerate(pool.imap_unordered(_stream_job, blocks), 1):
            num1 += r1[0]
            num2 += r2[0]
            gap = [min(gap[0], r1[1], r2[1]), max(gap[1], r1[2], r2[2])]
            if n % 50 == 0:
                print("  blocks", n, len(blocks), round(time.time() - t0), flush=True)
    for k, num in ((kmax - 1, num1), (kmax, num2)):
        P = Fraction(num, 216 * 8 ** (k - 1))
        res["levels"][k] = dict(P=str(P), P_float=float(P), matches_known=(P == KNOWN.get(k)) if k in KNOWN else None,
                                seconds=round(time.time() - t0, 1))
        print(k, res["levels"][k], flush=True)
    res["overlap_gap"] = dict(min_reversible=gap[0], max_violator=gap[1])
    Ps = [Fraction(res["levels"][k]["P"]) for k in range(1, kmax + 1)]
    res["ratios"] = [float(Ps[i + 1] / Ps[i]) for i in range(len(Ps) - 1)]
    res["denominators"] = [Ps[i].denominator for i in range(len(Ps))]
    res["denominator_law_2^(3k-1)3^[k odd]"] = [Ps[i].denominator == 2 ** (3 * (i + 1) - 1) * 3 ** ((i + 1) % 2)
                                                 for i in range(len(Ps))]
    print(res["levels"][kmax], res["ratios"], res["denominators"], flush=True)
    return res


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    res = run(kmax)
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
