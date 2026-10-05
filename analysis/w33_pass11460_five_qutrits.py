"""Pass 11460: five qutrits, one magic gate -- does the bad-class fraction tend to 1/8?

P_n = exact cell shares x sampled cell fractions (Pass 11437's estimator), at n = 5 with Pass 11421's sparse decider
(243-dimensional states, 59,049 frames; validated against the table decider at n = 2, 3 in Pass 11421).
Sanity at n = 5: the fixed-axis rule and the law F' (Passes 11420/11457) on targeted fixed-axis classes.
Comparison: P_1 = P_2 = 1/8, P_3 = 437/3276 = 0.13339 (exact), P_4 = 0.128 +- 0.005 (Pass 11437).
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11352_three_qutrit_geometry as GEO  # noqa: E402
import w33_pass11437_four_qutrit_cells as C  # noqa: E402

OUT = ROOT / "data" / "w33_pass11460_five_qutrits.json"


def _init_lowcap():
    C._init(5)
    C.L.CAP = 3 ** 9                      # fixed-axis classes have huge (S) spaces; report undecided instead of stalling


def _fixed_job(seed):
    D, gens = C._S["D"], C._S["gens"]
    rng = np.random.default_rng(seed)
    z, Om = D.z1, D.wl.Om
    M = GEO.random_symplectic(rng, gens, length=300)
    M = (GEO.map_to((M @ z) % 3, z, Om, rng) @ M) % 3
    g = D.good_frames(M)
    if g is None:
        return "undecided"
    bad = ~g
    xne = D.wl.labels[:, 0] != 0
    return "F' equality" if (bad == xne).all() else ("F' containment" if (bad >= xne).all() else "FAILS")


def run(per_kind=4000, fixed=40):
    res = dict(pass_id=11460)
    cells = C.sample(5, per_kind, seed0=114600000)
    res["n5"] = dict(sampled_cells=cells, estimate=C.estimate(5, cells))
    print(json.dumps(res["n5"], indent=1, default=str), flush=True)
    json.dump(res, open(OUT, "w"), indent=1, default=str)              # save the census before the slow check
    with Pool(6, initializer=_init_lowcap) as pool:
        res["n5_fixed_axis_F_prime"] = dict(Counter(pool.map(_fixed_job, range(114608000, 114608000 + fixed))))
    res["n5_fixed_axis_cap"] = "3^9"
    print(res["n5_fixed_axis_F_prime"], flush=True)
    p5, s5 = res["n5"]["estimate"]["P_estimate"], res["n5"]["estimate"]["stderr"]
    res["sequence"] = {"1": 0.125, "2": 0.125, "3": 437 / 3276, "4": [0.1284, 0.0045], "5": [p5, s5]}
    res["n5_z_vs_one_eighth"] = (p5 - 0.125) / s5
    res["n5_z_vs_437_3276"] = (p5 - 437 / 3276) / s5
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
