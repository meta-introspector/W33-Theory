"""Regression for Passes 11213-11216 (frozen certificates plus fast recomputations)."""
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def load(name):
    return json.loads((ROOT / "data" / name).read_text())


def test_11213_cubic_t_violation():
    import w33_pass11213_cubic_t_violation as P
    C = P.clifford1()
    assert len(C) == 216
    assert all(P.t_invertible(V, C) for V in C) and P.t_invertible(P.T1, C)
    assert not P.t_invertible(P.T1 @ P.X1, C)
    assert P.t_invertible(P.T1 @ P.S1, C)
    d = load("w33_pass11213_cubic_t_violation.json")
    assert d["shortest"]["depth"] == 1 and d["shortest"]["violating_words"] == 18
    f = [d["by_depth"][str(k)]["fraction"] for k in range(9)]
    assert f[0] == 0 and f[8] > 0.9


def test_11214_rules():
    d = load("w33_pass11214_projective_line_rule.json")
    assert set(d["ame6_harmonic_labellings_per_sample"]) == {120}
    assert sum(d["ame10_labellings_and_rule"].values()) == 71
    for k in d["ame10_labellings_and_rule"]:
        assert k.startswith("2 labellings") and ("(True, False), (False, True)" in k or "(False, True), (True, False)" in k)


def test_11214_field():
    import w33_pass11214_projective_line_rule as R
    K9 = R.Field(3, ext=True)
    i = K9.enc(0, 1)
    assert K9.mul(i, i) == K9.neg(K9.ONE)
    assert all(K9.mul(x, K9.inv(x)) == K9.ONE for x in range(1, 9))


def test_11215_cycle_index():
    import w33_pass11215_protected_qutrit_limit as P
    assert P.steinberg_checks(5) == (True, True)
    E = P.expected_c(5)
    assert E[1] == 1 and E[2] == Fraction(133, 180) and E[3] == Fraction(106927, 147420)
    assert E[4] == Fraction(142080247, 195832728)
    d = load("w33_pass11215_protected_qutrit_limit.json")
    assert d["limit"]["value"].startswith("0.7255353507990583494490873")


def test_11216_circulation():
    d = load("w33_pass11216_clock_circulation.json")
    assert d["3-4"]["current_eigenvalues_over_i_sqrt3"] == [24.0, 36.0]
    assert d["12-13"]["current_eigenvalues_over_i_sqrt3"] == [72.0, 324.0]
    assert d["3-4"]["oriented_3_cycles_per_point"] == 7168
