"""Regression for Passes 11498-11505 (the magic-axis law for all n; the good rules; exact P9; six qutrits extended; the
1/108 identity completed and the time-odd module; the rule-cell shares; the cubed-phase test)."""

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def load(name):
    return json.load(open(ROOT / "data" / name))


def test_11498_theorem_ingredients():
    d = load("w33_pass11498_magic_axis_law_all_n.json")
    for n, r in d["lemmas"].items():
        assert r["t_k_in_span_z1"] and r["r_Q"]["r_Q = 0"] == r["r_Q"]["tested"] > 0
    n2 = d["n2_all_classes"]
    assert n2["classes"] == 1296 and n2["solutions with QJv = -v"] == n2["solutions enumerated"] == 23328
    assert n2["decided: every reversible frame has omega(v,a) = 0"] == n2["decided"] == 1295
    n3 = d["n3_orbit_reps"]
    assert n3["solutions with QJv = -v"] == n3["solutions enumerated"] and n3["decided: every reversible frame has omega(v,a) = 0"] == n3["decided"]
    lin = d["linear_check_whole_solution_space"]
    for n in ("n2", "n3", "n4"):
        assert lin[n]["spaces with QJv = -v throughout"] == lin[n]["affine solution spaces"] > 0
    assert lin["n4"]["classes"] == 549


def test_11498_proof_step2_fast():
    """A = QJ with A = M s^-k A M and Q z1 = z1 gives A v = -v (on a random n = 2 class with M^2 z1 = z1)"""
    import w33_pass11350_linear_decider as L
    import w33_pass11498_magic_axis_law_all_n as T
    D = L.Decider(2)
    M = np.eye(4, dtype=np.int64)                     # identity: M^2 z1 = z1, v = 2 z1
    out = T.linear_check(D, [M])
    assert out["spaces with QJv = -v throughout"] == out["affine solution spaces"]


def test_11499_good_rules():
    d = load("w33_pass11499_good_rules.json")
    c = d["n2"]["counts"]
    assert c["M^2 z1 = -z1"] == 1944 and c["M z1 = -z1"] == 648
    assert c["M^2 z1 = -z1: hypothesis holds"] == c["M^2 z1 = -z1: hypothesis holds and all frames reversible"] == 1890
    assert c["M z1 = -z1: hypothesis holds"] == c["M z1 = -z1: hypothesis holds and all frames reversible"] == 620
    u = d["n2_failures_by_union"]
    assert u["M^2 z1 = -z1: union of k=0 criteria covers every frame"] == u["M^2 z1 = -z1: union frames all reversible (sound)"] == 54
    assert u["M z1 = -z1: union of k=0 criteria covers every frame"] == u["M z1 = -z1: union frames all reversible (sound)"] == 27
    assert u["M z1 = -z1: capped and unresolved"] == 0
    c3 = d["n3"]["counts"]
    for cell in ("M z1 = -z1", "M^2 z1 = -z1"):
        assert c3[f"{cell}: hypothesis holds"] == c3[f"{cell}: hypothesis holds and all frames reversible"]


def test_11500_exact_depth9():
    d = load("w33_pass11500_exact_depth9.json")
    assert d["levels"]["9"]["P"] == "11519231/201326592" and d["levels"]["8"]["matches_known"] is True
    assert all(d["denominator_law"]) and d["counts_over_18"][-1] == "11519231"
    assert d["control_depth8_words"]["agree"] and d["overlap_gap"]["max_prefilter_passing_violator"] == 0.0
    ck = d["dedupe_checks"]["7672320"]
    assert ck["at_1e5"] == ck["at_1e6"] == 5567184 <= ck["at_1e9"]


def test_11500_prefilter_fast():
    import w33_pass11422_depth5_exact as E
    import w33_pass11436_depth_reduction as RD
    import w33_pass11500_exact_depth9 as P9
    E._init()
    RD._init()
    CT, RT = RD._S["CT"], RD._S["RT"]
    U = np.einsum('rij,wjk->rwik', RT, CT).reshape(-1, 3, 3)          # all depth-2 words
    rev = P9.reversible(U)[0]
    assert Fraction(int(rev.sum()), len(U)) == Fraction(19, 32)


def test_11501_six_qutrits():
    d = load("w33_pass11501_six_qutrits_extended.json")
    e = d["estimate"]
    assert d["classes"] >= 2000 and min(d["completed_by_kind"].values()) >= 1000
    assert 0.05 < e["P_estimate"] < 0.25 and e["stderr"] < 0.012


def test_11502_rank_one_and_module():
    d = load("w33_pass11502_time_odd_module.json")
    r = d["rank_one_census"]
    assert r["rank_one"] == 1259712 and r["non_stabiliser"] == 0 and r["uniform_over_144_pairs"]
    assert r["pair_counts_min"] == r["pair_counts_max"] == 8748 and r["singular_values"] == [1.0]
    assert d["even_generator_degrees"] == [1, 3, 4, 5, 6, 7, 8, 9]
    assert d["odd_generator_degrees"] == [6, 7, 8, 9, 10]
    for k, v in d["module_structure"].items():
        for kept, dropped in v["gaps"]:
            assert kept > 1e-8 and dropped < 1e-12


def test_11503_cell_shares():
    d = load("w33_pass11503_rule_cell_shares.json")
    assert d["n2"]["all_match"] and d["n3"]["all_match"] and d["n3"]["group_order"] == 9170703360
    assert d["n3"]["shares"]["M^2 z1 = -z1"] == "3/728" and d["n3"]["shares"]["M^2 z1 = z1, v != 0"] == "1/364"
