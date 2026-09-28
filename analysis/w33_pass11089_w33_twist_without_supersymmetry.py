#!/usr/bin/env python3
"""Pass 11089: the W(3,3) twist without supersymmetry.

(1) SUSY-independence.  The obstructions of Passes 11024/11090 that close the Z3xZ3 route are statements about
    massless FERMIONS and LEFT-MOVING quantum numbers: a Dirac mass needs a conjugate (SM x hidden) pair, so the chiral
    excess of the hidden-refined index stays massless in any vacuum preserving the hidden gauge group, whether or not
    supersymmetry is broken (by F-terms, gaugino condensation, or anything that preserves the gauge group); the charge
    class c = 3Q + t is fixed by the gauge-lattice shift of each fixed point, so it is the same for any right-moving
    completion.  Breaking supersymmetry does not rescue those vacua.

(2) No non-supersymmetric prime-order W(3,3) orbifold.  A non-SUSY twin of a Z_N orbifold attaches (-1)^F to the
    twist: v -> v' = v + n, n an integer vector of odd sum (a 2 pi rotation, -1 on spinors).  For T6/Z3 the twin is
    v' = (1/3,1/3,1/3), v'^2 = 1/3, and modular invariance needs 3(V^2 - v'^2) in 2Z, i.e. 9(V1^2 + V2^2) = 3 mod 6.
    But for EVERY order-3 shift V of E8 (3V in the E8 lattice), 9V^2 is even: integral class 9V^2 = sum u_i^2 with
    sum u_i even; half-integral class 9V^2 = sum u_i(u_i + 1) + 2.  So no non-supersymmetric T6/Z3 of the E8 x E8
    string exists for any shifts -- in particular none with the W(3,3) (A8 Kac) twist.  Checked by enumeration and by
    orbifolder, which rejects the twin: '3 (V_M^2 - v_M^2) = 2.33 = 0 mod 2 failed'.

(3) Which W(3,3) families have non-SUSY twins.  A twin is non-supersymmetric iff no supercharge q in {+-1/2}^4 (even
    number of minus signs) has q.v' in Z for every generator; it keeps the SAME gauge shifts and Wilson lines iff
    N(v'^2 - v^2) in 2Z for each generator (and the cross term for Z_M x Z_N is unchanged mod 1).  The table below is
    computed from the orbifolder geometry twists.  Adding (-1)^F turns every sign combination +-v1 +-v2 +-v3 by an
    odd integer, so a twin can be non-SUSY only if no sign combination of v is an odd integer.
(4) Z6-I: all 87 W(3,3) Standard Models have non-SUSY twins with the same shifts and Wilson lines -- orbifolder loads
    all 87 with N = 0 (the originals with N = 1).  orbifolder 1.2 cannot build an N = 0 spectrum (its spectrum code
    groups states into supermultiplets and segfaults), so their spectra need a non-supersymmetric spectrum engine: the
    named next step.  Without superpartners there are no dimension-4 or dimension-5 proton-decay operators, the
    obstruction that closed the supersymmetric Z6-I route.
"""
from __future__ import annotations

import itertools
import json
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11089_w33_twist_without_supersymmetry.json"

# twists (components 1-3) of the geometries used by the nine W(3,3) families (orbifolder 1.2 Geometry/*.txt)
FAMILIES = {
    "Z3": (3, [(F(1, 3), F(1, 3), F(-2, 3))]),
    "Z3xZ3": (3, [(F(1, 3), F(0), F(-1, 3)), (F(1, 3), F(-1, 3), F(0))]),
    "Z6-I": (6, [(F(1, 6), F(1, 6), F(-1, 3))]),
    "Z6-II": (6, [(F(1, 6), F(1, 3), F(-1, 2))]),
    "Z12-I": (12, [(F(1, 12), F(-5, 12), F(1, 3))]),
    "Z2xZ6-I": (6, [(F(1, 2), F(0), F(-1, 2)), (F(0), F(1, 6), F(-1, 6))]),
    "Z2xZ6-II": (6, [(F(1, 2), F(0), F(-1, 2)), (F(1, 6), F(1, 6), F(-1, 3))]),
    "Z3xZ6": (6, [(F(1, 3), F(0), F(-1, 3)), (F(0), F(1, 6), F(-1, 6))]),
    "Z6xZ6": (6, [(F(1, 6), F(0), F(-1, 6)), (F(0), F(1, 6), F(-1, 6))]),
}
ORDERS = {"Z3": [3], "Z3xZ3": [3, 3], "Z6-I": [6], "Z6-II": [6], "Z12-I": [12], "Z2xZ6-I": [2, 6],
          "Z2xZ6-II": [2, 6], "Z3xZ6": [3, 6], "Z6xZ6": [6, 6]}


def supercharges():
    for s in itertools.product((F(1, 2), F(-1, 2)), repeat=4):
        if sum(1 for x in s if x < 0) % 2 == 0:
            yield s


def n_susy(twists):
    return sum(1 for q in supercharges()
               if all((q[1] * v[0] + q[2] * v[1] + q[3] * v[2]).denominator == 1 for v in twists))


def norm(v):
    return sum(x * x for x in v)


def twin_table():
    out = {}
    for fam, (_, twists) in FAMILIES.items():
        orders = ORDERS[fam]
        best = None
        for g in range(len(twists)):                      # attach (-1)^F to generator g
            for n in itertools.product(range(-2, 3), repeat=3):
                if sum(n) % 2 == 0:
                    continue
                tw = [tuple(x + (ni if j == g else 0) for x, ni in zip(v, n)) for j, v in enumerate(twists)]
                ns = n_susy(tw)
                if ns:
                    continue
                same_shifts = all((orders[j] * (norm(tw[j]) - norm(twists[j]))) % 2 == 0 for j in range(len(tw)))
                if len(tw) == 2:
                    from math import gcd
                    gg = gcd(*orders)
                    # gcd(M,N) (v1'.v2' - v1.v2) must be EVEN (orbifolder: "= 0 mod 2"); an earlier version tested only
                    # integrality and listed twists that fail it (corrected in Pass 11096)
                    same_shifts = same_shifts and (gg * (sum(a * b for a, b in zip(tw[0], tw[1])) -
                                                          sum(a * b for a, b in zip(twists[0], twists[1])))) % 2 == 0
                cand = dict(generator=g, n=list(n), twist=[[str(x) for x in v] for v in tw], same_shifts=same_shifts)
                if best is None or (same_shifts and not best["same_shifts"]):
                    best = cand
        out[fam] = dict(n_susy_original=n_susy(twists),
                        non_susy_twin_exists=best is not None,
                        twin=best,
                        odd_integer_sign_combination=[any((s1 * v[0] + s2 * v[1] + s3 * v[2]).denominator == 1 and
                                                          (s1 * v[0] + s2 * v[1] + s3 * v[2]) % 2 == 1
                                                          for s1, s2, s3 in itertools.product((1, -1), repeat=3))
                                                      for v in twists])
    return out


def nine_v_squared_even(bound=3):
    vals = set()
    for u in itertools.product(range(-bound, bound + 1), repeat=8):
        if sum(u) % 2 == 0:
            vals.add(sum(x * x for x in u) % 2)
    half = [F(2 * x + 1, 2) for x in range(-bound, bound)]
    for u in itertools.product(half, repeat=8):
        if sum(u) % 2 == 0:
            vals.add(int(sum(x * x for x in u)) % 2)
    return sorted(vals)


def order3_twin_obstruction():
    """for every order-3 generator of Z3, Z3xZ3 and every odd-sum n, 9 (v+n)^2 is odd (while 9V^2 is always even):
    no (-1)^F can sit on an order-3 generator modular-invariantly, whatever the shifts."""
    res = {}
    for fam in ("Z3", "Z3xZ3"):
        par = set()
        for v in FAMILIES[fam][1]:
            for n in itertools.product(range(-3, 4), repeat=3):
                if sum(n) % 2:
                    par.add(int(9 * norm(tuple(x + ni for x, ni in zip(v, n)))) % 2)
        res[fam] = sorted(par)
    return res


def main():
    parity = nine_v_squared_even()
    ob3 = order3_twin_obstruction()
    print("9 v'^2 mod 2 for order-3 twins:", ob3, "   9V^2 mod 2 for order-3 shifts:", parity)
    table = twin_table()
    for k, v in table.items():
        print(k, "N_orig", v["n_susy_original"], "non-SUSY twin", v["non_susy_twin_exists"],
              "same shifts", v["twin"]["same_shifts"] if v["twin"] else None, v["twin"]["twist"] if v["twin"] else "")
    orbifolder = {
        "Z3 W(3,3) model (A8 + (1/3,1/3,0^6)) with twist (1/3,1/3,1/3)": "rejected: 3 (V_M^2 - v_M^2) = 2.33 != 0 mod 2",
        "Z3 W(3,3) model, SUSY control (twist (1/3,1/3,-2/3))": "loads, N = 1, gauge SU(9) x SO(14) x U(1)",
        "Z6-I W(3,3) Standard Models, SUSY originals": "87/87 load, N = 1",
        "Z6-I W(3,3) Standard Models, twins (twist (1/6,1/6,2/3))": "87/87 load, N = 0; spectrum builder segfaults (no N=0 path in orbifolder 1.2)",
        "Z6-II twin by (0,0,1)": "loads with N = 1: supersymmetric again (sign combination 1/6+1/3-1/2 = 0)",
    }
    out = dict(pass_id=11089, nine_V_squared_mod_2_over_order3_shifts=parity, nine_vprime_squared_mod_2_order3_twins=ob3,
               twins=table, orbifolder=orbifolder)
    OUT.write_text(json.dumps(out, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
