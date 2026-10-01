"""Regression for Passes 11265-11269 (frozen certificates plus fast recomputations)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
FMIN = (1 + 2 * np.cos(2 * np.pi / 9)) / 3


def load(name):
    return json.loads((ROOT / "data" / name).read_text())


def test_11265_formula_and_bb():
    d = load("w33_pass11265_exact_fidelity.json")
    assert d["formula_max_deviation_vs_twirl"] < 1e-12
    assert d["bb_agrees_all"]
    t3 = d["three_qutrit"]
    for k in ("(T(x)T)SUM (x) I", "(T(x)T(x)T) SUM12"):
        assert t3[k]["search_exhausted"] and t3[k]["numerical_pruning_tolerance"] == 1e-12
        assert abs(t3[k]["F_T_lower"] - FMIN) < 1e-9 and abs(t3[k]["F_T_upper"] - FMIN) < 1e-9
    open_ = t3["(T(x)T(x)T) SUM23 SUM12"]
    assert not open_["search_exhausted"] and open_["F_T_lower"] < float(open_["F_T_upper"]) < 1


def test_11265_bb_fast():
    import w33_pass11227_two_qutrit_cp_mixing as P2
    import w33_pass11252_exact_reversibility as R
    import w33_pass11265_exact_fidelity as X
    R.WEYL[2] = R.Weyl(2)
    bb = X.best_fidelity(R.WEYL[2], np.kron(P2.T, P2.T) @ P2.SUM)
    assert bb[1] and abs(bb[0] - FMIN) < 1e-9


def test_11267_selector():
    d = load("w33_pass11267_tm1_selector.json")
    assert d["phi_dot_chi_on_vacuum_manifold"]["TM1"] == [-4, 4]
    assert d["phi_dot_chi_on_vacuum_manifold"]["theta13 = 0 column"] == [0]
    assert d["kappa_selects_TM1_both_signs"] and d["kappa_beats_wrong_sign_eps_5e-4"]
    assert d["boson_loops_select_TM1"] and d["fermion_loops_select_theta13_0"]


def test_11269_dictionary():
    d = load("w33_pass11269_g26_qutrit_dictionary.json")
    assert d["jacobian_equals_19440_SIC_x_stab2_exact"]
    assert d["jacobian_symbolic_identity"] == "Q[w]/(w^2+w+1)[x,y,z]"
    assert d["jacobian_nonzero_monomials"] == 42
    rg = d["reflection_group"]
    assert rg["order"] == 1296 and rg["collineation_order"] == 216 and rg["all_reflections_clifford"]
    inv = d["invariants"]
    assert inv["u6_is_sixth_SIC_power_sum"] and inv["u12_is_twelfth_SIC_power_sum"] and inv["u9_vs_SIC_product"] == "1"
    assert d["codex_control_point"]["regular"]


def test_11269_sic_numeric():
    w = np.exp(2j * np.pi / 3)
    sic = [np.array(v) / np.sqrt(2) for v in
           [(1, -w ** k, 0) for k in range(3)] + [(1, 0, -w ** k) for k in range(3)] + [(0, 1, -w ** k) for k in range(3)]]
    G = np.abs(np.array([[np.vdot(a, b) for b in sic] for a in sic])) ** 2
    assert np.allclose(G[~np.eye(9, dtype=bool)], 0.25)


def test_11265_residuals():
    d = load("w33_pass11265_exact_fidelity.json")
    assert d["levels_with_a_quantised_optimal_residual"] == [0.666666667, 0.712386014, 0.844029629]
    assert all(v["is_Fmin2"] and v["sumset_residual"] for v in d["fmin2_words"].values())


def test_11266_one_magic_gate():
    from fractions import Fraction
    d = load("w33_pass11266_one_magic_gate.json")
    assert d["one_qutrit_k1"]["probability"] == "1/12" and d["one_qutrit_k1"]["rule_holds"]
    assert Fraction(d["one_qutrit_k2"]["violating"], d["one_qutrit_k2"]["total"]) == Fraction(13, 32)
    assert Fraction(d["one_qutrit_k3"]["violating"], d["one_qutrit_k3"]["total"]) == Fraction(421, 768)
    assert d["two_qutrit_k1"]["samples"] == 200_000 and 0.05 < d["two_qutrit_k1"]["probability"] < 0.15


def test_11266_rule_fast():
    import w33_pass11266_one_magic_gate as M
    r = M.one_qutrit_rule()
    assert r["violating"] == 18 and r["rule_holds"]


def test_11268_perfect_edges():
    d = load("w33_pass11268_perfect_magic_ticks.json")["perfect_magic"]
    two = d["2"]["violating_by_magic_legs"]
    assert two["1"] == [0, 648] and two["3"] == [2592, 2592] and two["2"] == [216, 1944] and two["4"] == [648, 1296]
    e3 = d["3"]["exhaustive_edges"]
    assert e3["one_magic_leg"] == [0, 8748] and e3["five_magic_legs"] == [139968, 139968]
