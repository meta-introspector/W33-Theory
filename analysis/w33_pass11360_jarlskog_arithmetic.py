"""Pass 11360: the arithmetic of the substrate Jarlskog invariant -- exact values of J_6 on short one-qutrit words.

J_6(U) = avg_C |tr(U^dag C U C^dag)|^6 - avg_C |tr(U^dag C U^T C^dag)|^6 (Pass 11355).  For Clifford+T words the traces
lie in Q(zeta_9) up to powers of sqrt 3, so |tr|^2 lies in the real cubic field Q(cos 2 pi/9) (Pass 11253's argument)
and so does J_6.  Here J_6 is recomputed at 60 digits with the 216 Cliffords rebuilt exactly (Pass 11237) for every
violating word with one and two cubic gates (coset-reduced), and each distinct value is identified as
a + b c + e c^2, c = 2 cos(2 pi/9), by PSLQ (accepted only if the relation holds to 10^-50).  Reports the distinct
values ('the Jarlskog spectrum') and whether each is rational.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11213_cubic_t_violation as P1  # noqa: E402
import w33_pass11237_fidelity_levels as FL  # noqa: E402
import w33_pass11312_depth_law as DL  # noqa: E402
import w33_pass11355_substrate_jarlskog as JS  # noqa: E402

OUT = ROOT / "data" / "w33_pass11360_jarlskog_arithmetic.json"


def run():
    mp, H, S, X, Z, T = FL.mp_setup()
    mp.mp.dps = 60
    exact = FL.exact_cliffords(mp, [H, S, X, Z])
    Cf = list(np.array(P1.clifford1()))
    Cl = [exact[P1.key(V)] for V in Cf]
    Rf = DL.coset_reps(Cf)
    Rl = [exact[P1.key(V)] for V in Rf]
    c = 2 * mp.cos(2 * mp.pi / 9)

    def tr(M):
        return M[0, 0] + M[1, 1] + M[2, 2]

    def J6(Ue):
        Ud = Ue.H
        UT = Ue.T
        a = b = mp.mpf(0)
        for C in Cl:
            a += abs(tr(Ud * C * Ue * C.H)) ** 6
            b += abs(tr(Ud * C * UT * C.H)) ** 6
        return (a - b) / 216

    vals = {}
    words = []
    for i in range(216):
        words.append((1, Cf[i] @ P1.T1, Cl[i] * T))
    for r in range(24):
        for i in range(216):
            words.append((2, Rf[r] @ P1.T1 @ Cf[i] @ P1.T1, Rl[r] * T * Cl[i] * T))
    for k, Uf, Ue in words:
        jf = JS.J(Uf)
        if jf < 1e-9:
            continue
        key = (k, round(jf, 9))
        if key in vals:
            continue
        j = J6(Ue)
        rel = mp.pslq([j, 1, c, c ** 2], maxcoeff=10 ** 12, maxsteps=10 ** 6)
        ok = rel is not None and rel[0] != 0 and abs(rel[0] * j + rel[1] + rel[2] * c + rel[3] * c ** 2) < mp.mpf(10) ** -50
        if ok:
            a0, a1, a2 = (-mp.mpf(rel[i]) / rel[0] for i in (1, 2, 3))
            from fractions import Fraction
            abc = [str(Fraction(int(-rel[i]), int(rel[0]))) for i in (1, 2, 3)]
        else:
            abc = None
        vals[key] = dict(gates=k, J6=mp.nstr(j, 30), in_Q_cos_2pi_9=ok, a_b_e=abc,
                         rational=bool(ok and abc[1] == "0" and abc[2] == "0"))
        print(key, vals[key], flush=True)
    res = dict(pass_id=11360, distinct_values=len(vals), values=list(vals.values()),
               all_in_cubic_field=all(v["in_Q_cos_2pi_9"] for v in vals.values()))
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k != "values"}, indent=1))


if __name__ == "__main__":
    main()
