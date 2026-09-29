#!/usr/bin/env python3
"""Pass 11171: the temporal CGLMP value never reaches 4 in any finite dimension, and its growth beyond d = 10.

THEOREM (4 is never attained).  One qudit, state rho, rank-1 Lueders measurements A_x then B_y:
P(a,b|x,y) = <a_x|rho|a_x> |<b_y|a_x>|^2.  I_d = 4 needs every term (x,y) to equal 1, i.e. P supported on the 'correct'
relation b = a + s_xy with s_00 = s_01 = s_11 = 0, s_10 = 1.  Let S_x = supp p(.|x).  Then for a in S_0, |a_0> is an
eigenvector of both B-bases with labels a, a: |a_0> = |b0_a> = |b1_a>; for a in S_1, |a_1> = |b0_{a+1}> = |b1_a>.  Since
<a_x|rho|a_x> = 0 forces rho|a_x> = 0, supp rho lies in span{b1_a : a in S_0} and in span{b1_a : a in S_1}, hence in
span{b1_a : a in S_0 cap S_1} -- so S_0 cap S_1 is nonempty, and for a in it b0_a = b1_a = b0_{a+1} up to phase:
two distinct vectors of one orthonormal basis coincide.  Contradiction.  So I_d < 4 for every finite d, although a
classical system with memory (invasive) reaches 4 (Pass 11144).  The memory a quantum qudit can carry between the
two times -- the post-measurement vector |a_x> -- is too small to encode both a and x perfectly.
NUMERICS.  analysis/w33_pass11171_scan_temporal_cglmp_large_d.py: I_d = max over four bases of lambda_max(W) (the optimal
state is the top eigenvector), Riemannian gradient ascent with analytic gradients, 24 restarts per d; frozen in
data/w33_pass11171_temporal_cglmp_large_d.json.  Values are LOWER bounds (local optimisation); d <= 10 are checked
against the committed optima (Passes 11147, 11155, 11160).
RESULTS.  d = 12, 14, 16, 20, 24: 3.7147, 3.7562, 3.7873, 3.8310, 3.8528 (all < 4, monotone).  The d = 4, 5, 8, 9, 10
optima are reproduced (to 2e-4); d = 6, 7 IMPROVE to 3.4552, 3.5209 (committed 3.4531, 3.5168 were local optima --
Pass 11155 corrected).  d (4 - I_d) = 3.39, 3.41, 3.42, 3.42, 3.41, 3.40, 3.38 for d = 8..20: the deficit falls like
~3.4/d (fit for d >= 6: exponent 0.97; a free-limit fit gives 4.03), superseding the d^-0.75 fit to d <= 10 (Pass 11160).
Consistent with sup_d I_d = 4, never attained.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCAN = ROOT / "data" / "w33_pass11171_temporal_cglmp_large_d.json"
OUT = ROOT / "data" / "w33_pass11171_temporal_cglmp_summary.json"
COMMITTED = {int(k): v['temporal'] for k, v in json.loads((ROOT / "data" / "w33_pass11160_temporal_cglmp_high_d.json")
              .read_text())['table'].items() if int(k) >= 3}


def fits(ds, vals):
    ds, vals = np.asarray(ds, float), np.asarray(vals, float)
    # (a) power law to 4:  4 - I = c d^-alpha
    A = np.vstack([np.ones_like(ds), np.log(ds)]).T
    co, res_a, *_ = np.linalg.lstsq(A, np.log(4 - vals), rcond=None)
    pred_a = 4 - np.exp(co[0]) * ds ** co[1]
    # (b) power law to a free limit: I = L - c d^-alpha  (grid over L)
    best = None
    for L in np.linspace(vals.max() + 1e-3, 4.5, 2000):
        cb, *_ = np.linalg.lstsq(A, np.log(L - vals), rcond=None)
        pred = L - np.exp(cb[0]) * ds ** cb[1]
        r = float(((pred - vals) ** 2).sum())
        if best is None or r < best[0]:
            best = (r, float(L), float(-cb[1]))
    return dict(to4_alpha=float(-co[1]), to4_rss=float(((pred_a - vals) ** 2).sum()),
                free_limit=best[1], free_alpha=best[2], free_rss=best[0])


def summarize():
    d = json.loads(SCAN.read_text())
    table = {int(k): v['best'] for k, v in d.items()}
    check = {k: table[k] - COMMITTED[k] for k in table if k in COMMITTED}
    merged = {**COMMITTED, **{k: max(v, COMMITTED.get(k, -9)) for k, v in table.items()}}
    ds = sorted(merged)
    fit_all = fits(ds, [merged[k] for k in ds])
    big = [k for k in ds if k >= 6]
    fit_big = fits(big, [merged[k] for k in big])
    res = dict(pass_id=11171, scan_best={str(k): v for k, v in sorted(table.items())},
               validation_vs_committed={str(k): v for k, v in sorted(check.items())},
               merged={str(k): v for k, v in sorted(merged.items())},
               monotone=all(merged[a] < merged[b] for a, b in zip(ds, ds[1:])),
               below_4=all(v < 4 for v in merged.values()), fit_all=fit_all, fit_d_ge_6=fit_big)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
