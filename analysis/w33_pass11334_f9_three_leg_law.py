"""Pass 11334: the F9-linear three-leg law, exhaustively -- magic on three legs of an F9-linear perfect two-qutrit tick
always breaks substrate time reversal.

Pass 11311 sampled it (6000/6000); Codex's audit correctly classed that as sample evidence.  Here it is exhaustive.

REDUCTION (exact).  A dressing exponent e in Z_9 gives T^e = Z^{(e - e mod 3)/3} T^{e mod 3} (T^3 = Z).  Z-type Paulis
on output legs multiply V = W(a) V_M on the left; on input legs V Z^k = W(M z^k-label) V up to phase -- in both cases
only the Pauli frame a changes, and every frame is enumerated.  So all 64 * 81 * (4 * 6^3 * 3) = 13,436,928 three-leg
cases of Pass 11311 reduce to
        U = (T^{s1} (x) T^{s2}) W(a) V_M (T^{s3} (x) T^{s4}),  s in {1, 2} on the three magic legs, 0 on the fourth,
i.e. 64 classes * 81 frames * 4 choices of the plain leg * 2^3 = 165,888 exact decisions (Pass 11252 criterion), with
the canonical Weil unitary V_M.  The one-magic-leg cases (41,472) are decided the same way for comparison.
"""

from __future__ import annotations

import itertools
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402
import w33_pass11252_exact_reversibility as R  # noqa: E402
import w33_pass11311_perfect_edge_law as PE  # noqa: E402
import w33_pass11330_orbit_census as O  # noqa: E402

OUT = ROOT / "data" / "w33_pass11334_f9_three_leg_law.json"
_S = {}


def lt(s1, s2):
    return np.kron(np.linalg.matrix_power(P2.T, s1), np.linalg.matrix_power(P2.T, s2))


def _init(f9):
    R.WEYL[2] = R.Weyl(2)
    wl = R.WEYL[2]
    rng = np.random.default_rng(0)
    _S["V"] = {i: R.weil(wl, M, rng) for i, M in f9}


def _job(args):
    i, a, k = args
    wl = R.WEYL[2]
    V = wl.W[a] @ _S["V"][i]
    rng = np.random.default_rng(a)
    viol = 0
    total = 0
    if k == 3:
        patterns = []
        for plain in range(4):
            for ss in itertools.product((1, 2), repeat=3):
                s = list(ss)
                s.insert(plain, 0)
                patterns.append(s)
    else:
        patterns = []
        for leg in range(4):
            for s0 in (1, 2):
                s = [0, 0, 0, 0]
                s[leg] = s0
                patterns.append(s)
    for s in patterns:
        U = lt(s[0], s[1]) @ V @ lt(s[2], s[3])
        total += 1
        viol += R.decide(U, 2, rng)[0] is False
    return viol, total


def run():
    wl = R.Weyl(2)
    R.WEYL[2] = wl
    Ms, keys = O.all_symplectic(wl)
    reps = np.load(P2.CACHE)
    f9 = []
    for i in range(len(Ms)):
        if PE.is_local(reps[i]) or not PE.perfect(reps[i]):
            continue
        M = Ms[i]
        if (((M @ PE.JB - PE.JB @ M) % 3) == 0).all():
            f9.append((i, M))
    res = dict(pass_id=11334, f9_perfect_classes=len(f9))
    print(res, flush=True)
    out = {}
    with Pool(11, initializer=_init, initargs=(f9,)) as pool:
        for k in (3, 1):
            jobs = [(i, a, k) for i, _ in f9 for a in range(81)]
            r = pool.map(_job, jobs, chunksize=16)
            v = sum(x for x, _ in r)
            n = sum(y for _, y in r)
            out[str(k)] = dict(violating=int(v), cases=int(n), fraction=v / n)
            print(k, out[str(k)], flush=True)
    res["magic_legs"] = out
    res["three_leg_law_exhaustive"] = out["3"]["violating"] == out["3"]["cases"]
    res["covers_unreduced_cases"] = 64 * 81 * 4 * 216 * 3
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
