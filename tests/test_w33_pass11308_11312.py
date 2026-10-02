"""Regression for Passes 11308-11312 (frozen certificates plus fast exact recomputations)."""
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
FMIN = (1 + 2 * np.cos(2 * np.pi / 9)) / 3


def load(name):
    return json.loads((ROOT / "data" / name).read_text())


def test_11308_shared_level():
    d = load("w33_pass11308_shared_level_and_clifford_bound.json")
    i = d["identities"]
    assert i["trT_over3_equals_Fmin"] and i["T_tensor_T_trace_equals_sumset"] and i["equals_Fmin_squared"]
    assert d["holds_FT_ge_2Fcl2_minus_1"] and d["holds_FT_ge_Fcl2_on_sample"]
    assert d["holds_FT_ge_Fcl_two_qutrit_sample"] is False
    assert d["bounds"]["one_qutrit"]["violating_words"] == 1344 and d["bounds"]["two_qutrit"]["violating_words"] == 159


def test_11308_identity_fast():
    import w33_pass11308_shared_level_and_clifford_bound as M
    i = M.exact_identities()
    assert i["T_tensor_T_trace_equals_sumset"] and abs(i["T3_max_stabilizer_probability"] - FMIN ** 2) < 1e-14


def test_11310_not_u81():
    d = load("w33_pass11310_magic_orientation_cocycle.json")
    assert d["magic_groups"]["T"]["invariants"] == {"order": 81, "center": 3, "derived": 9, "gamma3": 3, "order9": 18}
    assert d["U81_order9"] == {"1": 36, "2": 36} and d["non_isomorphic"]
    assert d["gap"]["magic_sylow_X_S_T_mod_scalars"] == [81, 9] and d["gap"]["U81_s1"] == [81, 7]
    assert d["T_and_T_inverse_generate_same_group_mod_scalars"]


def test_11310_fast():
    import w33_pass11310_magic_orientation_cocycle as M
    G = M.generate([M.X, M.S, M.T_gate(1)])
    assert len(G) == 81 and M.invariants(G)["order9"] == 18 and M.u81_order9(1) == 36


def test_11311_edge_law():
    d = load("w33_pass11311_perfect_edge_law.json")
    assert d["class_sizes"] == {"local": 1152, "perfectF9": 64, "perfect": 13760, "other": 36864}
    c = d["cells"]
    assert all(c[f"{k}:k=0"]["violating"] == 0 for k in ("local", "perfectF9", "perfect", "other"))
    assert c["perfectF9:k=3"]["violating"] == c["perfectF9:k=3"]["samples"] == 6000
    assert not d["perfect_one_magic_leg_never"] and not d["perfect_three_magic_legs_always"]
    assert d["F9_perfect_three_magic_legs_always"]


def test_11309_one_gate_law():
    d = load("w33_pass11309_one_gate_law.json")
    assert d["bad_classes"] == 6480 == 51840 // 8
    assert d["bad_class_frame_counts"] == {"54": 5160, "72": 24, "81": 1296}
    assert d["violating_cliffords_exact"] == 385344 and d["probability_exact"] == "223/2430"
    assert "NOT_AFFINE" not in d["nonviolating_frames_affine_subspace_dims"]
    assert Fraction(5160 * 54 + 24 * 72 + 1296 * 81, 51840 * 81) == Fraction(223, 2430)


def test_11312_depth_law():
    d = load("w33_pass11312_depth_law.json")
    assert d["matches_pass_11266"]
    assert d["exact"]["4"]["probability"] == "1425/2048" and d["exact"]["4"]["words"] == 2985984
    s = d["sampled"]
    assert s["5"]["probability"] < s["6"]["probability"] < s["8"]["probability"] < s["12"]["probability"] < 1


def test_11312_reduction_fast():
    import w33_pass11312_depth_law as D
    D._init()
    assert len(D._S["R"]) == 24
    assert D._block((1, ())) == 18                                   # k = 1: 18 of 216
