"""Pass 11459: where does the one-qutrit reversible fraction's decay rate go at large depth?

Exact values to depth 7 (Pass 11436): ratios 0.648, 0.761, 0.673, 0.736, 0.693, 0.725 -- an odd/even oscillation damping
toward ~0.71; the sampled local rates of Pass 11422 crept up to 0.727 at depths 14-16.

CORRECTED DECIDER THRESHOLD (found here): reversible iff max overlap > 3 - 1e-7, not > 2.999 (see TIGHT).
Here: 10 million samples per depth at k = 12, 18, 24, 30 using the K-reduction of Pass 11436 (the left-most Clifford
uniform over 216, every inner gate uniform over the 8 coset representatives -- the same law as the full walk) and the
exact decider of Pass 11422.  Local rates are geometric means over 6 gates (even-to-even, so the odd/even oscillation
cancels).  This is a measurement: no transfer operator governs a Haar-null event (Pass 11372); the question is whether
the rate keeps creeping up.
"""

from __future__ import annotations

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11422_depth5_exact as E  # noqa: E402
import w33_pass11436_depth_reduction as RD  # noqa: E402

OUT = ROOT / "data" / "w33_pass11459_deep_decay.json"
# THRESHOLD.  Exactly reversible words have max overlap 3 to rounding (deficit < 1e-10); deep violating words come
# arbitrarily close by density (deficits ~1e-4..1e-3 at depth 24-30).  The loose cut 2.999 used in Pass 11422's sampled
# tail counts such near-misses as reversible (38% false at depth 24, 80% at depth 30) -- see threshold_sensitivity().
TIGHT = 3 - 1e-7


def _sample(args):
    k, seed, n = args[:3]
    lo = args[3] if len(args) > 3 else 0                # lo = 1 excludes the identity coset (no T-merging)
    RD._init()
    CT, RT = RD._S["CT"], RD._S["RT"]
    assert RD._S["reps"][0] in RD._S["K"]               # representative 0 is the identity class K
    rng = np.random.default_rng(seed)
    rev, done = 0, 0
    B = 512
    for _ in range(n // B):
        U = np.broadcast_to(np.eye(3, dtype=complex), (B, 3, 3)).copy()
        idx = rng.integers(lo, 8, size=(B, k - 1))
        for j in range(k - 1):
            U = np.einsum('wij,wjk->wik', RT[idx[:, j]], U)
        U = np.einsum('wij,wjk->wik', CT[rng.integers(216, size=B)], U)
        rev += int((E.overlaps(U) > TIGHT).sum())
        done += B
    return rev, done


def run(ks=(12, 18, 24, 30), per_k=10_000_000, nproc=8):
    res = dict(pass_id=11459, depths={})
    with Pool(nproc) as pool:
        for k in ks:
            parts = pool.map(_sample, [(k, 114590000 + 1000 * k + i, per_k // nproc) for i in range(nproc)])
            rev, n = sum(p[0] for p in parts), sum(p[1] for p in parts)
            p = rev / n
            res["depths"][str(k)] = dict(samples=n, reversible=rev, fraction=p, stderr=float(np.sqrt(p * (1 - p) / n)))
            print(k, res["depths"][str(k)], flush=True)
    ks = list(ks)
    loc = []
    for a, b in zip(ks, ks[1:]):
        pa, pb = res["depths"][str(a)]["fraction"], res["depths"][str(b)]["fraction"]
        sa, sb = res["depths"][str(a)]["stderr"], res["depths"][str(b)]["stderr"]
        r = (pb / pa) ** (1 / (b - a))
        se = r / (b - a) * np.sqrt((sa / pa) ** 2 + (sb / pb) ** 2)
        loc.append(dict(span=f"{a}->{b}", rate=float(r), stderr=float(se)))
    res["local_rates"] = loc
    return res


def _sens(args):
    k, seed, batches = args
    RD._init()
    CT, RT = RD._S["CT"], RD._S["RT"]
    rng = np.random.default_rng(seed)
    cuts = {"2.999": 2.999, "3-1e-6": 3 - 1e-6, "3-1e-8": 3 - 1e-8, "3-1e-10": 3 - 1e-10}
    cnt = {c: 0 for c in cuts}
    n = 0
    for _ in range(batches):
        B = 512
        U = np.broadcast_to(np.eye(3, dtype=complex), (B, 3, 3)).copy()
        idx = rng.integers(8, size=(B, k - 1))
        for j in range(k - 1):
            U = np.einsum('wij,wjk->wik', RT[idx[:, j]], U)
        U = np.einsum('wij,wjk->wik', CT[rng.integers(216, size=B)], U)
        ov = E.overlaps(U)
        for c, th in cuts.items():
            cnt[c] += int((ov > th).sum())
        n += B
    return cnt, n


def threshold_sensitivity(ks=(12, 24, 30), batches=1200, nproc=8):
    out = {}
    with Pool(nproc) as pool:
        for k in ks:
            parts = pool.map(_sens, [(k, 114597000 + 100 * k + i, batches // nproc) for i in range(nproc)])
            cnt = {c: sum(p[0][c] for p in parts) for c in parts[0][0]}
            out[str(k)] = dict(samples=sum(p[1] for p in parts), counts=cnt,
                               false_fraction_at_2999=1 - cnt["3-1e-8"] / max(cnt["2.999"], 1))
            print("sens", k, out[str(k)], flush=True)
    return out


def no_merging(ks=(12, 16, 20, 24), per_k=8_000_000, nproc=8):
    """MECHANISM TEST.  The identity coset lets adjacent T gates merge (T^2, and T^3 = Z1 is Clifford), lowering the
    effective magic of deep words.  The same walk with the identity coset excluded (7 inner choices): if its local rates
    stay flat, T-merging is what slows the decay."""
    out = dict(depths={})
    with Pool(nproc) as pool:
        for k in ks:
            parts = pool.map(_sample, [(k, 114595000 + 1000 * k + i, per_k // nproc, 1) for i in range(nproc)])
            rev, n = sum(p[0] for p in parts), sum(p[1] for p in parts)
            p = rev / n
            out["depths"][str(k)] = dict(samples=n, reversible=rev, fraction=p, stderr=float(np.sqrt(p * (1 - p) / n)))
            print("no-merge", k, out["depths"][str(k)], flush=True)
    ks = list(ks)
    out["local_rates"] = []
    for a, b in zip(ks, ks[1:]):
        pa, pb = out["depths"][str(a)]["fraction"], out["depths"][str(b)]["fraction"]
        sa, sb = out["depths"][str(a)]["stderr"], out["depths"][str(b)]["stderr"]
        r = (pb / pa) ** (1 / (b - a)) if pa > 0 and pb > 0 else float("nan")
        se = r / (b - a) * np.sqrt((sa / pa) ** 2 + (sb / pb) ** 2) if pa > 0 and pb > 0 else float("nan")
        out["local_rates"].append(dict(span=f"{a}->{b}", rate=float(r), stderr=float(se)))
    return out


def main():
    if "--sens" in sys.argv:
        res = json.load(open(OUT))
        res["threshold_sensitivity"] = threshold_sensitivity()
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        return
    if "--nomerge" in sys.argv:
        res = json.load(open(OUT))
        res["no_identity_coset_walk"] = no_merging()
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        print(json.dumps(res["no_identity_coset_walk"]["local_rates"], indent=1))
        return
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps(res["local_rates"], indent=1))


if __name__ == "__main__":
    main()
