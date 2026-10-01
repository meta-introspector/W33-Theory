#!/usr/bin/env python3
"""Pass 11229: is every AME(10,3) stabilizer state F9-linear?  An independent, larger census.

Pass 11224 showed that the 71 AME(10,3) graph states of master's census (300 tabu restarts) are all F9-linear up to
local Cliffords and all Glynn's code, and that the F9-type state is unique up to LC and relabelling (|Aut_LC| = 2880,
one orbit of 79 888 260 016 373 760).  Whether a non-F9-linear AME(10,3) stabilizer state exists was left open.

Every stabilizer state is LC-equivalent to a graph state, so a search over weighted qutrit graphs reaches every LC class.
This pass reruns master's tabu search (Pass 11186) with fresh seeds and more restarts and applies Pass 11224's
F9-linearity test and Schur-square test to every AME(10,3) graph found.  Outcome: reported in the certificate.  A
negative search is evidence, not proof; a single non-F9-linear hit would settle the question the other way.
"""
from __future__ import annotations

import json
import sys
import time
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11186_scan_ame10_tabu as T  # noqa: E402
import w33_pass11175_scan_five_qutrit as F  # noqa: E402
import w33_pass11224_ame10_f9_structure as A  # noqa: E402

OUT = ROOT / "data" / "w33_pass11229_ame10_f9_census.json"


def classify(w):
    G = F.to_mat(np.array(w))
    st = A.f9_structures(G)
    tr1 = [s for s in st if all(int(np.trace(A.ORD8[k]) % 3) == 1 for k in s)]
    out = dict(f9_structures=len(st))
    if tr1:
        rows, _ = A.to_f9_code(G, tr1[0])
        k, sq = A.schur_square_dim(rows)
        out.update(f9_dimension=k, schur=sq)
    return out


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    restarts = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    minutes = float(sys.argv[3]) if len(sys.argv) > 3 else 70
    t0 = time.time()
    seeds = [11229_000 + k for k in range(workers)]
    with Pool(workers) as pool:
        res = pool.map(T.worker, [(s, restarts, 4000, minutes) for s in seeds])
    dist, found = Counter(), []
    for d, f in res:
        dist.update(d)
        found += f
    old = {tuple(w) for w in json.loads(T.OUT.read_text())["graphs"]}
    new = [w for w in found if tuple(w) not in old]
    cls = [classify(w) for w in found]
    summary = Counter((c["f9_structures"], c.get("schur")) for c in cls)
    out = dict(pass_id=11229, seeds=seeds, restarts_run=sum(dist.values()), best_distribution={str(k): v for k, v in
               sorted(dist.items())}, ame_found=len(found), new_vs_master_census=len(new),
               classes={f"structures={a},schur={b}": n for (a, b), n in summary.items()},
               all_f9_linear_glynn=all(c["f9_structures"] == 4 and c.get("schur") == 10 for c in cls),
               graphs=found, seconds=round(time.time() - t0, 1))
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "graphs"}, indent=1))


if __name__ == "__main__":
    main()
