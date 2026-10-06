"""Pass 11535: what the ten exact reversible fractions P_1..P_10 say about the decay rate -- and what they do not.

Data: Passes 11312, 11422, 11436, 11488, 11500, 11513 (exact rationals, denominators 2^(3k-1) 3^[k odd]).
Tests: (1) integer sequence a_k = P_k 216 8^(k-1) / 18 -- any linear recurrence of order <= 4 with constant coefficients
(homogeneous or with constant term), on the whole sequence or on the two parity subsequences?  (2) two-step local rates
sqrt(P_(k+2)/P_k) on both parity chains; (3) the gap to the sampled constant 0.72 (Pass 11459) and to rho_(8,8) = 0.72276.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11535_decay_exact_data.json"


def load_P():
    d = json.load(open(ROOT / "data" / "w33_pass11513_exact_depth10.json"))
    P = {}
    for k in range(1, 11):
        if str(k) in d["levels"] and "P" in d["levels"][str(k)]:
            P[k] = Fraction(d["levels"][str(k)]["P"])
    d9 = json.load(open(ROOT / "data" / "w33_pass11500_exact_depth9.json"))
    for k in range(1, 10):
        if k not in P and str(k) in d9["levels"]:
            P[k] = Fraction(d9["levels"][str(k)]["P"])
    return [P[k] for k in range(1, 11)]


def recurrence_tests(seq, maxorder=4):
    out = {}
    for o in range(1, maxorder + 1):
        for const in (False, True):
            m = o + (1 if const else 0)
            rows = [list(seq[i:i + o]) + ([1] if const else []) for i in range(len(seq) - o)]
            rhs = [seq[i + o] for i in range(len(seq) - o)]
            if len(rows) <= m:
                out[f"order {o}{' +const' if const else ''}"] = "too few terms"
                continue
            A = sp.Matrix(rows[:m])
            if A.det() == 0:
                out[f"order {o}{' +const' if const else ''}"] = "singular fit"
                continue
            x = A.LUsolve(sp.Matrix(rhs[:m]))
            ok = all(sum(x[j] * rows[i][j] for j in range(m)) == rhs[i] for i in range(len(rows)))
            out[f"order {o}{' +const' if const else ''}"] = "FITS" if ok else "no"
    return out


def run():
    P = load_P()
    a = [int(P[k] * 216 * 8 ** k / 18) for k in range(10)]
    assert all(P[k] * 216 * 8 ** k / 18 == a[k] for k in range(10))
    res = dict(pass_id=11535, P=[str(p) for p in P], integers=a)
    res["recurrences_full"] = recurrence_tests(a)
    res["recurrences_odd_k"] = recurrence_tests(a[0::2], 2)
    res["recurrences_even_k"] = recurrence_tests(a[1::2], 2)
    two = [float(P[i + 2] / P[i]) ** 0.5 for i in range(8)]
    res["two_step_rates_odd_start"] = two[0::2]
    res["two_step_rates_even_start"] = two[1::2]
    res["one_step_ratios"] = [float(P[i + 1] / P[i]) for i in range(9)]
    res["rho88"] = float(max(np.roots([2592, -1332, -585, 140]).real))
    print(json.dumps(res, indent=1), flush=True)
    return res


def main():
    json.dump(run(), open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
