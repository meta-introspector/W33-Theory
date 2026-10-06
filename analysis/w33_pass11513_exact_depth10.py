"""Pass 11513: the exact reversible fraction at ten cubic gates.

Pass 11500's method one level deeper.  The depth-8 set has ~4.5e7 candidates (8 extensions of 5,567,184 distinct depth-7
operators), too many to deduplicate as matrices; instead each candidate is the pair (depth-7 index, extension r) carrying
two 64-bit hashes of its canonical key (phase-fixed entries rounded at 1e-6 and at 1e-5, integer-encoded, mixed by a
fixed multiplicative hash).  Deduplication = np.unique on the 1e-6 hash; merge check = the 1e-5 hash gives the SAME number
of classes (two distinct operators would have to agree to 1e-5 to merge, and distinct depth-8 operators are >~ 3^-8 apart);
splitting of equal operators by rounding is harmless (copies keep their weights).  Then every distinct depth-8 operator is
re-built from (index, r), its reversibility gives P_8 (validation), its 8 extensions P_9 (validation) and its 64 extensions
P_10.  Reversibility: Pass 11500's exact cubed-phase prefilter + Theorem 1 confirmation.
"""

from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11422_depth5_exact as E  # noqa: E402
import w33_pass11436_depth_reduction as RD  # noqa: E402
import w33_pass11488_exact_depth8 as X  # noqa: E402
import w33_pass11500_exact_depth9 as P9  # noqa: E402

OUT = ROOT / "data" / "w33_pass11513_exact_depth10.json"
KNOWN = dict(P9.KNOWN)
KNOWN[9] = Fraction(11519231, 201326592)
_C = np.random.default_rng(11513).integers(1, 2 ** 62, size=(2, 18), dtype=np.int64) | 1


def hashes(U):
    """two 64-bit hashes of the phase-fixed operator, at roundings 1e-6 and 1e-5"""
    f = U.reshape(len(U), 9)
    idx = np.argmax(np.abs(f) > 0.05, axis=1)
    ph = f[np.arange(len(f)), idx]
    g = f * (np.abs(ph) / ph)[:, None]
    g = np.concatenate([g.real, g.imag], axis=1)                                  # (W, 18)
    out = []
    for i, scale in enumerate((1e6, 1e5)):
        q = np.rint(g * scale).astype(np.int64)
        with np.errstate(over="ignore"):
            h = (q * _C[i]).sum(axis=1)
            h ^= (h >> 29)
            h *= np.int64(0x5DEECE66D)
            h ^= (h >> 31)
        out.append(h)
    return out


def prefilter_fast(U):
    """Pass 11500's cubed-phase prefilter, same decision, computing phases only for the (word, L) pairs whose moduli match"""
    c = np.einsum('pij,wji->wp', P9.WP.conj(), U) / 3
    mod = np.abs(c)
    okmod = (np.abs(mod[:, P9.PERM] - mod[:, None, :]) < 1e-9).all(axis=2)          # (W, 24), real arithmetic only
    out = np.zeros(len(U), bool)
    wi, li = np.nonzero(okmod)
    if len(wi) == 0:
        return out
    cc = c[wi]
    supp = mod[wi] > 1e-9
    cL = cc[np.arange(len(wi))[:, None], P9.PERM[li]]
    r3 = (cL / np.where(supp, cc, 1)) ** 3
    ref = r3[np.arange(len(wi)), np.argmax(supp, axis=1)][:, None]
    ok = np.where(supp, np.abs(r3 - ref), 0.0).max(axis=1) < 1e-7
    out[wi[ok]] = True
    return out


def reversible(U, chunk=20000):
    """exact: fast prefilter, then Theorem 1 on the survivors"""
    rev = np.zeros(len(U), bool)
    lo, hi, kept = 3.0, 0.0, 0
    for i in range(0, len(U), chunk):
        Uc = U[i:i + chunk]
        keep = prefilter_fast(Uc)
        kept += int(keep.sum())
        if keep.any():
            ov = X.overlaps(Uc[keep])
            r = ov > P9.TIGHT
            rev[i + np.flatnonzero(keep)[r]] = True
            if r.any():
                lo = min(lo, float(ov[r].min()))
            if (~r).any():
                hi = max(hi, float(ov[~r].max()))
    return rev, lo, hi, kept


def _stream_job(args):
    U, w = args
    if "CF" not in E._S:
        E._init()
        RD._init()
    RT = RD._S["RT"]
    res = []
    rev, lo, hi, _ = reversible(U)
    res.append((int((w * rev).sum()), lo, hi))
    X1 = np.einsum('rij,wjk->wrik', RT, U).reshape(-1, 3, 3)
    w1 = np.repeat(w, 8)
    rev, lo, hi, _ = reversible(X1)
    res.append((int((w1 * rev).sum()), lo, hi))
    tot, lo2, hi2 = 0, 3.0, 0.0
    for r in range(8):                                                             # 64 extensions, 8 at a time
        X2 = np.einsum('rij,wjk->wrik', RT, X1[r::8]).reshape(-1, 3, 3)
        rev, lo, hi, _ = reversible(X2)
        tot += int((np.repeat(w, 8) * rev).sum())
        lo2, hi2 = min(lo2, lo), max(hi2, hi)
    res.append((tot, lo2, hi2))
    return res


def run(nproc=6):
    RD._init()
    E._init()
    CT, RT = RD._S["CT"], RD._S["RT"]
    t0 = time.time()
    U = CT.copy()
    w = np.ones(len(U), dtype=np.int64)
    res = dict(pass_id=11513, levels={}, dedupe_checks={})
    for k in range(2, 8):                                                          # distinct depth-7 set (Pass 11500)
        U = np.einsum('rij,wjk->wrik', RT, U).reshape(-1, 3, 3)
        w = np.repeat(w, 8)
        U, w = P9.dedupe(U, w)
        res["levels"][k] = dict(distinct_operators=len(U))
        print(k, len(U), round(time.time() - t0), flush=True)
    assert len(U) == 5567184
    # depth 8 by hashes
    n7 = len(U)
    H6 = np.empty(8 * n7, np.int64)
    H5 = np.empty(8 * n7, np.int64)
    B = 200000
    for i in range(0, n7, B):
        Xc = np.einsum('rij,wjk->wrik', RT, U[i:i + B]).reshape(-1, 3, 3)       # order: (operator, r)
        h6, h5 = hashes(Xc)
        H6[8 * i:8 * i + len(Xc)] = h6
        H5[8 * i:8 * i + len(Xc)] = h5
    u6, first, inv = np.unique(H6, return_index=True, return_inverse=True)
    n5 = len(np.unique(H5))
    W8 = np.bincount(inv.ravel(), weights=np.repeat(w, 8).astype(np.float64), minlength=len(u6))
    assert abs(W8.sum() - 216 * 8 ** 7) < 0.5
    W8 = np.rint(W8).astype(np.int64)
    res["dedupe_checks"]["8"] = dict(candidates=int(8 * n7), classes_at_1e6=int(len(u6)), classes_at_1e5=int(n5))
    assert n5 == len(u6), res["dedupe_checks"]["8"]
    print("depth 8 distinct", len(u6), "hash check", n5, round(time.time() - t0), flush=True)
    del H6, H5, inv
    idx7, rr = first // 8, first % 8
    res["levels"][8] = dict(distinct_operators=int(len(u6)))
    # stream
    tot = [0, 0, 0]
    gap = [3.0, 0.0]
    C = 2000
    starts = list(range(0, len(first), C))

    def blocks():
        for s in starts:
            Ub = np.einsum('wij,wjk->wik', RT[rr[s:s + C]], U[idx7[s:s + C]])
            yield Ub, W8[s:s + C]

    with Pool(nproc) as pool:
        for n, out in enumerate(pool.imap(_stream_job, blocks(), chunksize=1), 1):
            for i, (num, lo, hi) in enumerate(out):
                tot[i] += num
                gap = [min(gap[0], lo), max(gap[1], hi)]
            if n % 500 == 0:
                print("  blocks", n, len(starts), round(time.time() - t0), flush=True)
    for i, k in enumerate((8, 9, 10)):
        P = Fraction(tot[i], 216 * 8 ** (k - 1))
        res["levels"].setdefault(k, {})
        res["levels"][k].update(P=str(P), P_float=float(P), matches_known=(P == KNOWN.get(k)) if k in KNOWN else None)
        print(k, res["levels"][k], flush=True)
    assert res["levels"][8]["matches_known"] and res["levels"][9]["matches_known"]
    res["overlap_gap"] = dict(min_reversible=gap[0], max_prefilter_passing_violator=gap[1])
    Ps = [Fraction(v) for v in [*(str(KNOWN[k]) for k in range(1, 10)), res["levels"][10]["P"]]]
    res["denominator_law"] = [Ps[i].denominator == 2 ** (3 * (i + 1) - 1) * 3 ** ((i + 1) % 2) for i in range(10)]
    res["counts_over_18"] = [str(Ps[i] * 216 * 8 ** i / 18) for i in range(10)]
    res["two_step_rates"] = [float(Ps[i + 2] / Ps[i]) ** 0.5 for i in range(8)]
    res["seconds"] = round(time.time() - t0)
    print(json.dumps({k: v for k, v in res.items() if k != "levels"}), flush=True)
    return res


def main():
    res = run(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
