"""Pass 11437: four qutrits cell by cell -- why the one-gate bad-class fraction drops from 437/3276 (n = 3) to about
1/8 (n = 4, Pass 11421).

The bad-class fraction is an exact average over the cells of the magic axis (Pass 11373):
    P_n = sum_cells (share of Sp(2n,3)) x (bad fraction in the cell),
with shares fixed by vector counts among the 3^(2n) - 1 nonzero images M z1:
    fixed / reversed: 1 each;  collinear (omega(z1, Mz1) = 0, Mz1 != +-z1): 3^(2n-1) - 3;  non-collinear: 3^(2n) - 3^(2n-1).
(The collinear vectors split further by M^-1 z1 into 'same line ...' and 'collinear, different lines'.)  Here each
cell is sampled directly at n = 4: M = h g with g a long random word and h a transvection map sending g z1 to a
uniformly random vector of the target cell (Pass 11352's map_to), decided by Pass 11421's sparse decider.  The
same estimator at n = 3 is compared with the exact cell fractions of Pass 11373 (a control on the sampler).
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11350_linear_decider as L  # noqa: E402
import w33_pass11352_three_qutrit_geometry as GEO  # noqa: E402
import w33_pass11421_four_qutrits as F4  # noqa: E402

OUT = ROOT / "data" / "w33_pass11437_four_qutrit_cells.json"
_S = {}


def _init(n):
    L.CAP = 3 ** 11
    _S["n"] = n
    _S["D"] = F4.SparseDecider(n)
    _S["gens"] = GEO.gen_mats(n)


def _random_vector(rng, kind, z, Om, N2):
    while True:
        y = rng.integers(3, size=N2)
        if not y.any() or ((y - z) % 3 == 0).all() or ((y + z) % 3 == 0).all():
            continue
        col = int(z @ Om @ y) % 3 == 0
        if (kind == "collinear") == col:
            return y


def _job(args):
    seed, kind = args
    D, gens = _S["D"], _S["gens"]
    rng = np.random.default_rng(seed)
    z, Om = D.z1, D.wl.Om
    g = GEO.random_symplectic(rng, gens, length=60 * _S["n"])
    y = _random_vector(rng, kind, z, Om, D.N2)
    M = (GEO.map_to((g @ z) % 3, y, Om, rng) @ g) % 3
    cell = GEO.cell(M, z, Om)
    good = D.good_frames(M)
    if good is None:
        return dict(kind=kind, cell=cell, undecided=True)
    return dict(kind=kind, cell=cell, bad=bool((~good).any()))


def sample(n, per_kind, nproc=6, seed0=0):
    jobs = [(seed0 + i, kind) for kind in ("collinear", "non-collinear") for i in range(per_kind)]
    with Pool(nproc, initializer=_init, initargs=(n,)) as pool:
        rows = pool.map(_job, jobs, chunksize=4)
    out = defaultdict(Counter)
    for r in rows:
        out[r["cell"]]["undecided" if r.get("undecided") else ("bad" if r["bad"] else "good")] += 1
    return {k: dict(v) for k, v in out.items()}


def estimate(n, cells):
    """P_n from exact shares and sampled within-cell fractions (collinear cells weighted by their sampled split)"""
    tot = 3 ** (2 * n) - 1
    col = 3 ** (2 * n - 1) - 3
    nonc = 3 ** (2 * n) - 3 ** (2 * n - 1)
    colcells = {k: v for k, v in cells.items() if k != "non-collinear"}
    ncol = sum(v.get("bad", 0) + v.get("good", 0) for v in colcells.values())
    bad_col = sum(v.get("bad", 0) for v in colcells.values())
    vc = cells.get("non-collinear", {})
    nnc = vc.get("bad", 0) + vc.get("good", 0)
    p_col, p_nc = bad_col / ncol, vc.get("bad", 0) / nnc
    se = np.sqrt((col / tot) ** 2 * p_col * (1 - p_col) / ncol + (nonc / tot) ** 2 * p_nc * (1 - p_nc) / nnc)
    P = (1 + col * p_col + nonc * p_nc) / tot                         # + the always-bad fixed cell (share 1/tot)
    return dict(collinear_bad_fraction=p_col, noncollinear_bad_fraction=p_nc, P_estimate=float(P), stderr=float(se),
                shares=dict(collinear=str(Fraction(col, tot)), noncollinear=str(Fraction(nonc, tot))))


def run():
    res = dict(pass_id=11437)
    # control at n = 3 against Pass 11373's exact cells
    c3 = sample(3, 1500, seed0=114370000)
    res["n3_control"] = dict(sampled_cells=c3, estimate=estimate(3, c3),
                             exact=dict(collinear="65/432", noncollinear="10/81", P="437/3276"))
    print(json.dumps(res["n3_control"], indent=1, default=str), flush=True)
    c4 = sample(4, 3000, seed0=114375000)
    res["n4"] = dict(sampled_cells=c4, estimate=estimate(4, c4))
    print(json.dumps(res["n4"], indent=1, default=str), flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
