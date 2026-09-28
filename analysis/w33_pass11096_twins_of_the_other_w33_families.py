#!/usr/bin/env python3
"""Pass 11096: the non-supersymmetric twins of the other W(3,3) families -- Z12-I, Z2xZ6-I, Z3xZ6, Z6xZ6.

A. Correction to Pass 11089.  For Z_M x Z_N, modular invariance needs gcd(M,N)(v1'.v2' - v1.v2) to be EVEN; Pass 11089
   tested only integrality, and its twin twists for Z2xZ6-I, Z3xZ6 and Z6xZ6 fail (orbifolder rejects them:
   "2 (V_N x V_M - v_N x v_M) = 1 = 0 mod 2 failed").  Correct twins exist for all three (and for Z2xZ6-II); Z3xZ3 still
   has none.  The representatives used: v2 -> (0, 7/6, -1/6) (Z2xZ6-I, Z3xZ6), v1 -> (7/6, 0, -1/6) (Z6xZ6),
   v -> (1/12, -5/12, -2/3) (Z12-I).
B. Exact right-mover tachyon levels of every twin, sector by sector (pure Python).
C. Frozen evidence (patched orbifolder, analysis/orbifolder_n0_drivers): all 333 W(3,3) Standard-Model twins load with
   N = 0 and are anomaly-free up to one Green--Schwarz U(1); tachyons from the mass-level engine at every level of B.
D. Z2xZ6-I: the 8 tachyon-free twins, with the SUSY parent's hypercharge applied to the twin's fermions, never have the
   Standard-Model chiral content.
"""
from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from fractions import Fraction as F
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11089_w33_twist_without_supersymmetry as Q  # noqa: E402
import w33_pass11092_nonsusy_z6i_twins_are_tachyonic as P  # noqa: E402

EVID = ROOT / "data" / "w33_pass11096_family_twin_evidence.json"
CHIRAL = ROOT / "data" / "w33_pass11096_z2xz6i_twin_chiral_content.json"
OUT = ROOT / "data" / "w33_pass11096_twins_of_the_other_w33_families.json"

Z = F(0)
ORIGINAL = {
    "Z12-I": ((12,), [(F(1, 12), F(-5, 12), F(1, 3))]),
    "Z2xZ6-I": ((2, 6), [(F(1, 2), Z, F(-1, 2)), (Z, F(1, 6), F(-1, 6))]),
    "Z2xZ6-II": ((2, 6), [(F(1, 2), Z, F(-1, 2)), (F(1, 6), F(1, 6), F(-1, 3))]),
    "Z3xZ6": ((3, 6), [(F(1, 3), Z, F(-1, 3)), (Z, F(1, 6), F(-1, 6))]),
    "Z6xZ6": ((6, 6), [(F(1, 6), Z, F(-1, 6)), (Z, F(1, 6), F(-1, 6))]),
    "Z3xZ3": ((3, 3), [(F(1, 3), Z, F(-1, 3)), (F(1, 3), F(-1, 3), Z)]),
}
TWIN = {
    "Z12-I": [(F(1, 12), F(-5, 12), F(-2, 3))],
    "Z2xZ6-I": [(F(1, 2), Z, F(-1, 2)), (Z, F(7, 6), F(-1, 6))],
    "Z3xZ6": [(F(1, 3), Z, F(-1, 3)), (Z, F(7, 6), F(-1, 6))],
    "Z6xZ6": [(F(7, 6), Z, F(-1, 6)), (Z, F(1, 6), F(-1, 6))],
}
OLD_1089 = {  # the twists Pass 11089 listed (integrality test only)
    "Z2xZ6-I": [(F(1, 2), Z, F(-1, 2)), (F(-2), F(-11, 6), F(-7, 6))],
    "Z3xZ6": [(F(1, 3), Z, F(-1, 3)), (F(-2), F(-11, 6), F(-7, 6))],
    "Z6xZ6": [(F(-11, 6), F(-2), F(-7, 6)), (Z, F(1, 6), F(-1, 6))],
}


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def conditions(orders, orig, tw):
    """modular-invariance conditions a twin with the SAME shifts and Wilson lines must satisfy"""
    order_ok = all((orders[j] * (dot(tw[j], tw[j]) - dot(orig[j], orig[j]))) % 2 == 0 for j in range(len(tw)))
    mixed = None
    if len(tw) == 2:
        mixed = gcd(*orders) * (dot(tw[0], tw[1]) - dot(orig[0], orig[1]))
    return dict(n_susy=Q.n_susy(tw), order_conditions=order_ok,
                mixed_value=str(mixed) if mixed is not None else None,
                mixed_even=(mixed % 2 == 0) if mixed is not None else None,
                spinorial_order_unchanged=all(sum(o * x for x in v) % 2 == 0 for o, v in zip(orders, tw)))


def count_twins(orders, orig):
    """single-generator (-1)^F insertions n (|n_i| <= 2, odd sum) giving a consistent non-SUSY twin with the same shifts"""
    good = 0
    for g in range(len(orig)):
        for n in itertools.product(range(-2, 3), repeat=3):
            if sum(n) % 2 == 0:
                continue
            tw = [tuple(x + (ni if j == g else 0) for x, ni in zip(v, n)) for j, v in enumerate(orig)]
            c = conditions(orders, orig, tw)
            if c["n_susy"] == 0 and c["order_conditions"] and (c["mixed_even"] in (None, True)):
                good += 1
    return good


def rm_levels(orders, tw):
    tw4 = [(Z,) + tuple(v) for v in tw]
    if len(tw4) == 1:
        tw4.append((Z,) * 4)
        orders = (orders[0], 1)
    out = {}
    for k in range(orders[0]):
        for l in range(orders[1]):
            vk = tuple(k * a + l * b for a, b in zip(tw4[0], tw4[1]))
            om = [P.frac(x) for x in vk[1:]]
            a_R = F(-1, 2) + sum(w * (1 - w) for w in om) / 2
            nums = [w for w in om if w] + [1 - w for w in om if w]
            osc = P.osc_counts(nums, F(1, 2))
            lev = Counter()
            for qs, n2 in P.so8_weights_near(list(vk), 2 * (F(1, 2) - a_R)):
                for N, c in osc.items():
                    m = n2 / 2 + N + a_R
                    if m < 0:
                        lev[str(m)] += c
            if lev:
                out[f"{k},{l}"] = dict(lev)
    return out


def main():
    A = {}
    for fam, (orders, orig) in ORIGINAL.items():
        A[fam] = dict(consistent_single_generator_twins=count_twins(orders, orig))
        if fam in TWIN:
            A[fam]["used_twin"] = conditions(orders, orig, TWIN[fam])
        if fam in OLD_1089:
            A[fam]["pass11089_twin"] = conditions(orders, orig, OLD_1089[fam])
        print(fam, A[fam], flush=True)
    B = {fam: rm_levels(ORIGINAL[fam][0], tw) for fam, tw in TWIN.items()}
    for fam, v in B.items():
        print(fam, 'tachyonic right-mover sectors', v, flush=True)
    ev = json.loads(EVID.read_text())
    chiral = json.loads(CHIRAL.read_text())
    z2 = ev["Z2xZ6-I"]["per_model"]
    tf = [k for k, v in z2.items() if v["total"] == 0]
    sm_like = {"(3,2)_1/6": 3, "(-3,1)_-2/3": 3, "(-3,1)_1/3": 3, "(1,2)_-1/2": 3, "(1,1)_1": 3}
    sm_like_conj = {"(-3,2)_-1/6": 3, "(3,1)_2/3": 3, "(3,1)_-1/3": 3, "(1,2)_1/2": 3, "(1,1)_-1": 3}
    D = {k: dict(label=z2[k]["label"], net_chiral=chiral[k]["net_chiral"],
                 net_fractional_colourless=chiral[k]["net_fractional_colourless"],
                 is_three_sm_generations=chiral[k]["net_chiral"] in (sm_like, sm_like_conj)) for k in tf}
    summary = {fam: dict(models=e["models"], n0=e["nsusy"].get("0", 0), anomaly_free_up_to_one_gs_u1=e["anomaly_free_up_to_one_gs_u1"],
                         tachyonic=e["tachyonic"], tachyon_free=e["tachyon_free"],
                         susy_parent_twisted_states_at_minus_1_6=e["susy_parent_twisted_states_at_minus_1_6"])
               for fam, e in ev.items()}
    summary["Z2xZ6-I_tachyon_free_with_three_sm_generations_under_parent_hypercharge"] = sum(d["is_three_sm_generations"] for d in D.values())
    summary["Z2xZ6-I_all_29_with_three_sm_generations_under_parent_hypercharge"] = sum(
        1 for v in chiral.values() if v["net_chiral"] in (sm_like, sm_like_conj))
    OUT.write_text(json.dumps(dict(pass_id=11096, modular_invariance=A, right_mover_tachyon_levels=B, summary=summary,
                                   z2xz6i_tachyon_free=D), indent=1, sort_keys=True))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
