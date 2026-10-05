"""Pass 11500: the exact reversible fraction at nine cubic gates.

Method of Pass 11488 (distinct operators carrying word multiplicities; deduplicated through depth 7, depths 8 and 9 streamed
as the 8 and 64 extensions of the depth-7 set), with one change: a cheap EXACT-NECESSARY prefilter before the overlap test.

PREFILTER (Pass 11252's criterion, cubed).  U reversible => c(Lp) = mu omega^<b,Lp> c(p) for an anti-symplectic L of F3^2
(24 of them), where c(p) = tr(W(p)^dag U)/3.  Cubing removes omega^<b,Lp>, so a reversible U must have |c(Lp)| = |c(p)| and
(c(Lp)/c(p))^3 constant on the support for some L.  Words failing this are violating; words passing it are confirmed by
Theorem 1 (overlap = 3, cut 3 - 1e-7).  Measured on 20 000 depth-8 words: the prefilter passes exactly the reversible ones.
"""

from __future__ import annotations

import itertools
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

OUT = ROOT / "data" / "w33_pass11500_exact_depth9.json"
TIGHT = X.TIGHT
KNOWN = dict(X.KNOWN)
KNOWN[8] = Fraction(665275, 8388608)


DEDUPE_CHECKS = {}


def dedupe(U, w):
    """Pass 11488's weighted deduplication with the merge check made one-sided: keys at 1e-5 and 1e-6 must give the same
    count (no merging of distinct operators, which are >~ 3^-k apart); 1e-9 may only SPLIT equal operators (floating noise
    at depth >= 7), which is harmless because split copies keep their own weights."""
    keys = X._canon(U, 6)
    _, first, inv = np.unique(keys, axis=0, return_index=True, return_inverse=True)
    W = np.bincount(inv.ravel(), weights=w, minlength=len(first))
    n5 = np.unique(X._canon(U, 5), axis=0).shape[0]
    n9 = np.unique(X._canon(U, 9), axis=0).shape[0]
    assert n5 == len(first) <= n9, (n5, len(first), n9)
    DEDUPE_CHECKS[len(U)] = dict(at_1e5=n5, at_1e6=len(first), at_1e9=n9)
    return U[first], W.astype(np.int64)


def _weyl_tables():
    w = np.exp(2j * np.pi / 3)
    Xm = np.roll(np.eye(3), 1, axis=0)
    Z = np.diag([1, w, w * w])
    labs = [(a, b) for a in range(3) for b in range(3)]
    Wp = np.array([w ** ((2 * a * b) % 3) * np.linalg.matrix_power(Xm, a) @ np.linalg.matrix_power(Z, b) for a, b in labs])
    anti = [np.array([[a, b], [c, d]]) for a, b, c, d in itertools.product(range(3), repeat=4) if (a * d - b * c) % 3 == 2]
    P = np.array([[labs.index(tuple(int(x) for x in (S @ np.array(p)) % 3)) for p in labs] for S in anti])
    return Wp, P


WP, PERM = _weyl_tables()


def prefilter(U):
    c = np.einsum('pij,wji->wp', WP.conj(), U) / 3
    mod = np.abs(c)
    cL = c[:, PERM]
    okmod = (np.abs(np.abs(cL) - mod[:, None, :]) < 1e-9).all(axis=2)
    supp = mod > 1e-9
    r3 = (cL / np.where(supp, c, 1)[:, None, :]) ** 3
    ref = r3[np.arange(len(U)), :, np.argmax(supp, axis=1)][:, :, None]          # ratio at the first support point
    dev = np.where(supp[:, None, :], np.abs(r3 - ref), 0.0).max(axis=2)
    return (okmod & (dev < 1e-7)).any(axis=1)


def reversible(U, chunk=20000):
    if len(U) > chunk:
        parts = [reversible(U[i:i + chunk], chunk) for i in range(0, len(U), chunk)]
        return (np.concatenate([p[0] for p in parts]), min(p[1] for p in parts), max(p[2] for p in parts),
                sum(p[3] for p in parts))
    return _reversible(U)


def _reversible(U):
    """exact: prefilter, then Theorem 1 on the survivors; returns (mask, min overlap of reversible, max overlap of
    prefilter-passing violators)"""
    keep = prefilter(U)
    rev = np.zeros(len(U), bool)
    lo, hi = 3.0, 0.0
    if keep.any():
        ov = X.overlaps(U[keep])
        r = ov > TIGHT
        rev[np.flatnonzero(keep)[r]] = True
        if r.any():
            lo = float(ov[r].min())
        if (~r).any():
            hi = float(ov[~r].max())
    return rev, lo, hi, int(keep.sum())


def _stream_job(args):
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
    for Xs, ww in ((X1, w1), (X2, w2)):
        rev, lo, hi, kept = reversible(Xs)
        out.append((int((ww * rev).sum()), lo, hi, kept, len(Xs)))
    return out


def run(kmax=9, nproc=6):
    RD._init()
    E._init()
    CT, RT = RD._S["CT"], RD._S["RT"]
    U = CT.copy()
    w = np.ones(len(U), dtype=np.int64)
    res = dict(pass_id=11500, levels={})
    gap = [3.0, 0.0]
    # control: the prefilter+confirm decider equals Theorem 1 on 20 000 random depth-8 words
    rng = np.random.default_rng(11500)
    V = CT[rng.integers(216, size=20000)]
    for _ in range(7):
        V = np.einsum('wij,wjk->wik', RT[rng.integers(8, size=20000)], V)
    full = X.overlaps(V) > TIGHT
    mine = reversible(V)[0]
    pre = prefilter(V)
    res["control_depth8_words"] = dict(words=20000, agree=bool((full == mine).all()), reversible=int(full.sum()),
                                       prefilter_passing=int(pre.sum()))
    assert res["control_depth8_words"]["agree"]
    print(res["control_depth8_words"], flush=True)
    for k in range(1, kmax - 1):
        t0 = time.time()
        if k > 1:
            U = np.einsum('rij,wjk->wrik', RT, U).reshape(-1, 3, 3)
            w = np.repeat(w, 8)
            U, w = dedupe(U, w)
        rev, lo, hi, _ = reversible(U)
        gap = [min(gap[0], lo), max(gap[1], hi)]
        P = Fraction(int(w[rev].sum()), 216 * 8 ** (k - 1))
        assert int(w.sum()) == 216 * 8 ** (k - 1)
        res["levels"][k] = dict(distinct_operators=len(U), P=str(P), matches_known=(P == KNOWN.get(k)) if k in KNOWN else None,
                                seconds=round(time.time() - t0, 1))
        print(k, res["levels"][k], flush=True)
        if k in KNOWN:
            assert P == KNOWN[k]
    t0 = time.time()
    tot = [0, 0]
    kept = [0, 0]
    B = 1000
    blocks = [(U[i:i + B], w[i:i + B]) for i in range(0, len(U), B)]
    with Pool(nproc) as pool:
        for n, out in enumerate(pool.imap_unordered(_stream_job, blocks), 1):
            for i, (num, lo, hi, kp, _) in enumerate(out):
                tot[i] += num
                kept[i] += kp
                gap = [min(gap[0], lo), max(gap[1], hi)]
            if n % 100 == 0:
                print("  blocks", n, len(blocks), round(time.time() - t0), flush=True)
    for i, k in enumerate((kmax - 1, kmax)):
        P = Fraction(tot[i], 216 * 8 ** (k - 1))
        res["levels"][k] = dict(P=str(P), P_float=float(P), matches_known=(P == KNOWN.get(k)) if k in KNOWN else None,
                                prefilter_survivors=kept[i], seconds=round(time.time() - t0, 1))
        print(k, res["levels"][k], flush=True)
    assert res["levels"][kmax - 1]["matches_known"] in (True, None)
    res["dedupe_checks"] = {str(k): v for k, v in DEDUPE_CHECKS.items()}
    res["overlap_gap"] = dict(min_reversible=gap[0], max_prefilter_passing_violator=gap[1])
    Ps = [Fraction(res["levels"][k]["P"]) for k in range(1, kmax + 1)]
    res["ratios"] = [float(Ps[i + 1] / Ps[i]) for i in range(len(Ps) - 1)]
    res["two_step_rates"] = [float(Ps[i + 2] / Ps[i]) ** 0.5 for i in range(len(Ps) - 2)]
    res["denominators"] = [Ps[i].denominator for i in range(len(Ps))]
    res["denominator_law"] = [Ps[i].denominator == 2 ** (3 * (i + 1) - 1) * 3 ** ((i + 1) % 2) for i in range(len(Ps))]
    res["counts_over_18"] = [str(Ps[i] * 216 * 8 ** i / 18) for i in range(len(Ps))]
    print(json.dumps({k: v for k, v in res.items() if k != "levels"}), flush=True)
    return res


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    res = run(kmax, nproc)
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
