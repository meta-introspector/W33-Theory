"""Regression for Passes 11199-11206 (frozen certificates plus cheap recomputations)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def load(name):
    return json.loads((ROOT / "data" / name).read_text())


def test_11199_local_spreads():
    d = load("w33_pass11199_ame10_local_spreads.json")
    c = d["counts"]
    assert c["spreads_1"] == c["regular_True"] == c["block_labellings_0"] == 8520
    assert c["regulus_True|block"] == 2 * 8520 and c["regulus_True|split_A1"] == 3 * 8520
    assert "regulus_True|split_A0" not in c


def test_11199_one_graph_recomputed():
    import w33_pass11199_ame10_local_spreads as S
    r = S.run(1)
    assert r["counts"]["spreads_1"] == 120 and r["counts"]["regular_True"] == 120


def test_11200_family():
    d = load("w33_pass11200_sign_method_ame6_ame8.json")
    assert d["ame6"]["all_isomorphic"] and set(d["ame6"]["aut_orders"]) == {120}
    a8 = d["ame8"]
    assert a8["with_duality"] in ("OPTIMAL", "FEASIBLE") and a8["n_solutions"] == 30
    assert a8["aut_orders_first40"] == {"1344": 30} and a8["four_sets_minus_singleton_form_S3_4_8"]
    assert 30 * 1344 == 40320


def test_11202_reflections():
    d = load("w33_pass11202_reflection_double_six.json")
    assert d["outer_involutions"] == {"12": 36, "24": 540}
    assert d["double_six_each"] and d["distinct_double_sixes"] == 36
    assert d["fixed_splits_are_planes_of_fixed_lines"] and d["fixed_geometry_is_GQ22"] and d["swaps_factors_in_all_15"]


def test_11203_certified():
    d = load("w33_pass11203_certified_bnb.json")
    assert d["complete"] and d["A_up_encloses"] and d["negative_control_rejected"]
    assert d["optimum_neighbourhood_accepted"] and d["target_below_4_over_sqrt15"]
    import w33_pass11203_certified_bnb as B
    assert B.TARGET_LO < 4 / np.sqrt(15)
    # a box at the corner of the simplex is closed by the rigorous bound
    assert B.monotone_up(np.array([[0.9, 0.0, 0.0, 0.0]]), np.array([[0.91, 0.01, 0.01, 0.01]]))[0] < B.TARGET_LO


def test_11205_one_way_clock():
    d = load("w33_pass11205_one_way_clock.json")
    assert all(d["reverse_is_conjugate"].values())
    assert all(d["directed_have_complex_eigenvalues"][k] for k in ("3", "4", "7", "9", "12", "13"))
    assert not d["directed_have_complex_eigenvalues"]["1"]
    ev = [complex(a, b) for a, b in d["spectra"]["3"]]
    assert any(abs(e - complex(-14, 18 * np.sqrt(3))) < 1e-5 for e in ev)
    assert d["commute_3_4"] and not d["commute_3_12"]


def test_11206_t_violation():
    d = load("w33_pass11206_t_violation.json")
    assert d["n2"]["all_T_invertible"] and d["n3"]["all_T_invertible"]
    assert d["n2"]["only_antiunitary_fraction"] == "77/162" and d["n3"]["only_antiunitary_fraction"] == "84013/157464"
    assert d["python_check_n2"]["every_class_inverted_by_an_anti_symplectic_map"]
    assert all(o % 3 == 0 for o in d["n2"]["non_real_orders"] + d["n3"]["non_real_orders"])


def test_11201_chiral_gauss_sum():
    d = load("w33_pass11201_twist_chiral_invariant.json")
    assert d["exact_matches_contraction"] and d["same_orientation_patterns_real"] and d["imag_equals_Q_sqrt3_over_243"]
    assert d["all_members_imag_in_units_i_sqrt3_over_243"] == {"3": {"-1.0": 256}, "4": {"1.0": 256}}
    assert d["pair_6912"]["flips"] and not d["symmetrised_flips"]
    assert abs(abs(d["pair_6912"]["value_12"][1]) - np.sqrt(3) / 729) < 1e-10


def test_11201_evaluator_one_case():
    import pytest
    import w33_pass11182_paper_ticks_mereology as PT
    import w33_pass11189_oriented_holonomy as H
    import w33_pass11201_twist_chiral_invariant as C
    if not H.CACHE.exists():
        pytest.skip("orbit labels (~25 min orbit computation) not cached")
    Bs = PT.load_facs()
    lab = np.load(H.CACHE)
    v = C.gauss_inv(C.lagrangian(Bs[np.flatnonzero(lab == 4)[0]] % 3), C.CHIRAL)
    assert abs(v - 1j * np.sqrt(3) / 243) < 1e-12


def test_11204_second_law():
    d = load("w33_pass11204_second_law_counting.json")
    assert d["validation_n2"] and d["validation_n3"]
    e = d["estimates"]
    assert e["4"]["fraction"] < 0.01 and e["5"]["arrow_free"] == 0
    f = [1.0, 133 / 360, 55241 / 884520, e["4"]["fraction"]]
    assert all(f[i + 1] / f[i] < f[i] / f[i - 1] for i in range(1, 3))   # accelerating decay


def test_11204_decomposition_known_positive():
    import w33_pass11204_second_law_counting as S
    J = S.form(2)
    V = S.all_vectors(2)
    I4 = np.eye(4, dtype=np.int64)
    assert S.decomposes(S.invariant_planes(I4, V, J), 2, J)
