#!/usr/bin/env python3
"""Pass 11220: how many qubits a random Clifford tick protects -- E_n[c] exactly through n = 6, and the limit.

Pass 11217 proved A(S) = n - c(S) for qubits and gave the closed form c = m1(1)/2 + chi2 m2(1) + m1(x^2+x+1) (exact on
all 320 classes with n <= 5).  Master's Pass 11215 computed the qutrit mean E_n[c] for every n and its limit
0.72553535...  For qubits the cycle-index route needs Hesselink's index; this pass instead sums over GAP's conjugacy
classes, now including Sp(12,2) (477 classes, analysis/gap/w33_pass11220_class_reps_q2_n6.g).

Results (exact fractions in the certificate):
    n = 2  29/40                     = 0.725
    n = 3  5981/8960                 = 0.667522
    n = 4  12945139/19496960         = 0.663957
    n = 5  6791923649683/10212039720960 = 0.665090
    n = 6  180847576047477961/271885345530839040 = 0.665161
Successive differences -0.0575, -0.00356, +0.00113, +0.0000713: the last ratio is ~1/16, so the limit is
0.665166 +- 0.00001 (geometric extrapolation; not a proof).  Qubits protect fewer than qutrits (0.72554).
The n = 6 value uses the closed form (brute-force-verified on all 320 classes with n <= 5).  `--brute` checks it at
n = 6 over all 1 397 760 nondegenerate planes of F_2^12, but the exhaustive search explodes on classes with many
invariant planes and c < 6; it was not completed in this session.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11208_arrow_universality as U  # noqa: E402
import w33_pass11217_qubit_arrow_law as Q  # noqa: E402

OUT = ROOT / "data" / "w33_pass11220_qubit_mean_protection.json"
N6 = ROOT / "data" / "w33_pass11220_gap_class_reps_q2_n6.txt"


def class_file(n):
    return U.GAP_SMALL if n <= 4 else (U.GAP_BIG if n == 5 else N6)


def exact_mean(n):
    cls, head = U.load_classes(class_file(n), 2, n)
    order = head[0]
    assert sum(c["size"] for c in cls) == order
    E = Fraction(0)
    dist = {}
    for c in cls:
        cc = Q.c_formula(c["S"] % 2, n)
        E += Fraction(c["size"] * cc, order)
        dist[cc] = dist.get(cc, 0) + c["size"]
    return dict(n=n, classes=len(cls), group_order=str(order), E_c=str(E), E_c_float=float(E),
                c_distribution={str(k): str(Fraction(v, order)) for k, v in sorted(dist.items())},
                p_arrow_free=str(Fraction(dist.get(n, 0), order)), p_max_arrow=str(Fraction(dist.get(0, 0), order)))


def run(brute=False):
    rows = [exact_mean(n) for n in range(2, 7)]
    vals = [r["E_c_float"] for r in rows]
    d = [b - a for a, b in zip(vals, vals[1:])]
    ratio = d[-1] / d[-2]
    limit = vals[-1] + d[-1] * ratio / (1 - ratio)
    res = dict(pass_id=11220, means=rows, differences=d, last_ratio=ratio, extrapolated_limit=limit,
               qutrit_limit_master_11215=0.72553535079905834944)
    prev = json.loads(OUT.read_text()) if OUT.exists() else {}
    if "n6_brute_force" in prev:
        res["n6_brute_force"] = prev["n6_brute_force"]
    if brute:
        import numpy as np
        cls, _ = U.load_classes(N6, 2, 6)
        J = U.form(6, 2)
        pts = U.points(12, 2)
        idx = {tuple(x): i for i, x in enumerate(pts)}
        W = (pts @ J @ pts.T) % 2
        planes = []
        for a in range(len(pts)):
            for b in np.flatnonzero(W[a, a + 1:]) + a + 1:
                c = idx[tuple((pts[a] + pts[b]) % 2)]
                if a < b < c:
                    planes.append(np.stack([pts[a], pts[b]], 1))
        Pl = np.array(planes)
        bad = sum(Q.c_brute(c["S"] % 2, Pl, 6) != Q.c_formula(c["S"] % 2, 6) for c in cls)
        res["n6_brute_force"] = dict(planes=len(Pl), classes=len(cls), mismatches=int(bad))
    return res


def main():
    res = run(brute="--brute" in sys.argv)
    OUT.write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k != "means"}, indent=1))
    for r in res["means"]:
        print(r["n"], r["classes"], r["E_c"], r["E_c_float"])


if __name__ == "__main__":
    main()
