"""Pass 11358: the three-qutrit one-gate mechanism -- solution structure of (S) against the verdict.

At n = 2 (Pass 11351) the verdict on the non-collinear cell was decided by the number of symplectic solutions of
(S) M s^-k Q (J M J) = Q, Q z1 = z1 (bad only with exactly 3 solutions, then half of them, exchanged by -I).  Here the
same statistics are collected at n = 3 on uniformly sampled Sp(6,3) classes (Pass 11352's sampler), with:
  * the cell of the magic axis, the number of symplectic solutions for each k, the number of distinct affine pieces in
    (F) and whether their union is affine;
  * the parity partner -M (its verdict, from the same decider);
to see which n = 2 mechanisms persist and where the extra badness (0.1333 vs 1/8) comes from.
"""

from __future__ import annotations

import itertools
import json
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11352_three_qutrit_geometry as G  # noqa: E402

OUT = ROOT / "data" / "w33_pass11358_three_qutrit_mechanism.json"
CAP = 3 ** 10


def nsol(D, M, k):
    s = D.symplectic_solutions(M, k)
    if s is None:
        return 0
    x0, b = s
    if 3 ** len(b) > CAP:
        return None
    n = 0
    N2 = D.N2
    for cs in itertools.product(range(3), repeat=len(b)):
        q = x0.copy()
        for c, v in zip(cs, b):
            q = (q + c * v) % 3
        n += D.is_symplectic(q.reshape(N2, N2))
    return n


def _job(seed):
    D = G._S["D"]
    rng = np.random.default_rng(seed)
    M = G.random_symplectic(rng, G._S["gens"])
    cell = G.cell(M, D.z1, D.wl.Om)
    ns = tuple(nsol(D, M, k) for k in range(3))
    g = D.good_frames(M)
    gm = D.good_frames((-M) % 3)
    if g is None or gm is None:
        return None
    return dict(cell=cell, nsol=ns, bad=bool((~g).any()), bad_minus=bool((~gm).any()),
                frames=int((~g).sum()))


def run(n=3000):
    with Pool(10, initializer=G._init) as pool:
        rows = [r for r in pool.map(_job, [900000 + i for i in range(n)], chunksize=10) if r]
    by_cell_sol = defaultdict(Counter)
    parity = defaultdict(Counter)
    for r in rows:
        key = f"{r['cell']} | nsol={r['nsol']}"
        by_cell_sol[key]["bad" if r["bad"] else "good"] += 1
        parity[r["cell"]][f"{int(r['bad'])}{int(r['bad_minus'])}"] += 1
    res = dict(pass_id=11358, sampled=n, decided=len(rows),
               bad_fraction=sum(r["bad"] for r in rows) / len(rows),
               cell_and_solution_count=dict(sorted({k: dict(v) for k, v in by_cell_sol.items()}.items(),
                                                   key=lambda kv: -sum(kv[1].values()))),
               parity_pairing_by_cell={k: dict(v) for k, v in parity.items()})
    print(json.dumps(res, indent=1), flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
