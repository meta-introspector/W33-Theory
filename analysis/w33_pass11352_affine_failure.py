"""Pass 11352 (supplement): the affine law fails at n = 3 -- locate the classes whose reversible frames are not an
affine subspace, and confirm them independently with the Weyl criterion of Pass 11252.

Re-scan the uniform seeds of Pass 11352 (both batches) restricted to the 'same line, other' cell (the cell is computed
first and is cheap), decide those classes with the linear decider, and for every class whose good-frame set is not an
affine subspace: (i) record its size and test whether it is a union of hyperplanes; (ii) decide 12 frames with the Weyl
criterion (6 inside the good set, 6 outside) and compare.
"""

from __future__ import annotations

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11252_exact_reversibility as R  # noqa: E402
import w33_pass11352_three_qutrit_geometry as G  # noqa: E402

OUT = ROOT / "data" / "w33_pass11352_affine_failure.json"


def is_affine(points):
    keys = {tuple(p) for p in points}
    P = np.array(points)
    return all(tuple((x + y - z) % 3) in keys for x in P[:40] for y in P[:40] for z in P)


def _scan(seed):
    D = G._S["D"]
    rng = np.random.default_rng(seed)
    M = G.random_symplectic(rng, G._S["gens"])
    if G.cell(M, D.z1, D.wl.Om) != "same line, other":
        return None
    g = D.good_frames(M)
    if g is None:
        return None
    ng = int(g.sum())
    if ng in (0, 729, 243, 81, 27, 9, 3, 1):
        labels = D.wl.labels[g]
        if ng == 0 or is_affine(labels):
            return None
    labels = D.wl.labels
    good = labels[g]
    # Weyl confirmation
    r2 = np.random.default_rng(seed + 7)
    inside = [int(x) for x in r2.choice(np.flatnonzero(g), 6, replace=False)]
    outside = [int(x) for x in r2.choice(np.flatnonzero(~g), 6, replace=False)]
    V = D.weil(M)
    agree = 0
    for a in inside + outside:
        U = D.wl.W[a] @ V @ D.T1
        rev = R.decide(U, 3, np.random.default_rng(a))[0] is True
        agree += rev == bool(g[a])
    # union-of-hyperplanes test: count affine hyperplanes contained in the good set
    hyper = 0
    for l in labels[1:]:
        if tuple(l) != tuple(labels[np.flatnonzero(l)[0]] * 0 + l):
            pass
    lv = labels[1:]
    lv = lv[[next(j for j in range(6) if v[j]) for v in lv] and [v[next(j for j in range(6) if v[j])] == 1 for v in lv]]
    dots = (labels @ lv.T) % 3                                       # (729, 364)
    for h in range(lv.shape[0]):
        for c in range(3):
            if g[dots[:, h] == c].all():
                hyper += 1
    union = np.zeros(729, bool)
    for h in range(lv.shape[0]):
        for c in range(3):
            if g[dots[:, h] == c].all():
                union |= dots[:, h] == c
    return dict(seed=seed, good_frames=ng, affine=False, weyl_checks=12, weyl_agree=agree,
                hyperplanes_inside=hyper, good_set_is_union_of_those_hyperplanes=bool((union == g).all()),
                good_size=len(good))


def run():
    seeds = [500000 + i for i in range(8000)] + [600000 + i for i in range(24000)]
    with Pool(10, initializer=G._init) as pool:
        rows = [r for r in pool.map(_scan, seeds, chunksize=50) if r]
    res = dict(pass_id=11352, supplement="affine law at n = 3", scanned_seeds=len(seeds), non_affine_classes=rows)
    print(json.dumps(res, indent=1), flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
