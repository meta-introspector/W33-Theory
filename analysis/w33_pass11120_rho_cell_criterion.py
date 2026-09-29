#!/usr/bin/env python3
"""Pass 11120: toward a closed-form rho-tachyon criterion -- a sharp necessary condition from the Wilson-line cells, and
why no single lattice-norm invariant suffices.

For every model (491: the 104 of Pass 11095 + 387 of the Pass 11108 rescan) and every Wilson-line cell (N, kappa)
(N in Z3^2 the winding classes of the two Wilson-line tori, kappa = (pi + V0).W mod 1 in (Z3)^2), the minimal norm
l^2_min(N, kappa) of the beta-projected gauge coset E8^2 + V0 + N.W with l^2 < 2 is computed (the only norms that can
carry a tachyon, since Delta = -1/2 + pR^2/2 with l^2 + pL^2 - pR^2 = 1).  Frozen: data/w33_pass11120_rho_cell_invariants_491.json.

Results:
  * minimal norms over N != 0 fall in five patterns (min, min over single-torus N, min over mixed N):
        (1, 1, 1) 185 | (1, 1, 5/3) 95 | (1, 5/3, 1) 120 (5 of them rho-safe) | (5/3, 5/3, 5/3) 91
    all 5 rho-safe models sit in (1, 5/3, 1);
  * sharper: every rho-safe model has EXACTLY 2 cells of norm 1 (both with mixed N and kappa nonzero in one torus only)
    and 10 cells of norm 5/3 -- the sparsest structure present; 15 tachyonic models share it (model 10968 even has the same
    norm-1 cells as the safe model 8666);
  * the universal rho levels: Delta = -1/6 <=> pR^2 = 2/3 and Delta = -1/18 <=> pR^2 = 8/9 (the WL-torus Narain vectors at
    the SU(3) point with the order-3 offsets).
Reading: a sharp necessary condition (sparsest cell structure) but NOT a sufficient one: the decision needs the torus
offsets o_t(N, kappa) = kappa_t/3 + V0.W_t + W_t.(N.W)/2 and the torus norms at rho -- the finite cell x torus lookup the
explicit enumerator (Pass 11107) performs.  No single gauge-lattice norm invariant separates rho-safe models.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11120_rho_cell_invariants_491.json"
OUT = ROOT / "data" / "w33_pass11120_rho_cell_criterion.json"


def feats(r):
    m = {eval(k): v for k, v in r['minl2'].items() if eval(k) != (0, 0)}
    single = [v for (a, b), v in m.items() if (a == 0) != (b == 0) and v is not None]
    mixed = [v for (a, b), v in m.items() if a and b and v is not None]
    allv = [v for v in m.values() if v is not None]
    return (round(min(allv), 3), round(min(single), 3), round(min(mixed), 3))


def cells(r, val):
    return [(N, k) for N, d in r['full'].items() for k, v in d.items() if abs(v - val) < 1e-6]


def summarize():
    res_all = json.loads(DATA.read_text())
    pat = Counter()
    for r in res_all:
        pat[(str(feats(r)), r['rho'] is None)] += 1
    safe = [r for r in res_all if r['rho'] is None]
    sparse = [r for r in res_all if len(cells(r, 1.0)) == 2 and len(cells(r, 5 / 3)) == 10]
    res = dict(pass_id=11120, models=len(res_all),
               patterns={f"{k[0]}|safe={k[1]}": v for k, v in pat.items()},
               safe_models=[r['model'] for r in safe],
               safe_all_sparse=all(len(cells(r, 1.0)) == 2 and len(cells(r, 5 / 3)) == 10 for r in safe),
               sparse_models=len(sparse), sparse_tachyonic=sum(1 for r in sparse if r['rho'] is not None),
               safe_norm1_cells_mixed_single_kappa=all(all(all(eval(N)) and sum(1 for x in eval(k) if x) == 1
                                                           for N, k in cells(r, 1.0)) for r in safe))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    summarize()
