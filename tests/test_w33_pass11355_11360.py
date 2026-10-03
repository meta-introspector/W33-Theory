"""Regression for Passes 11355-11360 (frozen certificates plus fast exact recomputations)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def load(name):
    return json.loads((ROOT / "data" / name).read_text())


def test_11355_complete_one_qutrit():
    d = load("w33_pass11355_substrate_jarlskog.json")
    assert d["identity_J6_equals_half_norm_sq_max_dev"] < 1e-10 and d["J2_J4_max_abs_on_random_words"] < 1e-10
    o = d["one_qutrit"]
    assert [o[k]["J6_positive"] for k in ("1", "2", "3", "4")] == [18, 2106, 68202, 2077650]
    assert all(o[k]["counts_equal"] and o[k]["per_word_mismatches"] == 0 for k in o)
    assert abs(o["1"]["min_positive_J6"] - 9 / 8) < 1e-12 and d["J6_complete_on_all_checked"]


def test_11355_reversal_transpose_fast():
    import w33_pass11213_cubic_t_violation as P1
    import w33_pass11252_exact_reversibility as R
    import w33_pass11355_substrate_jarlskog as J
    R.WEYL[1] = R.Weyl(1)
    rng = np.random.default_rng(1)
    for _ in range(40):
        C1, C2 = J.CF[rng.integers(216)], J.CF[rng.integers(216)]
        U = C1 @ P1.T1 @ C2 @ P1.T1
        assert (J.J(U) > 1e-9) == (R.decide(U, 1, rng)[0] is False)
        assert abs(J.J(U, 1)) < 1e-10 and abs(J.J(U, 2)) < 1e-10


def test_11356_two_qutrits():
    d = load("w33_pass11356_jarlskog_two_qutrits.json")
    assert d["checked"] == 240 and d["mismatches"] == 0 and d["violators"] == 173
    assert d["named"]["(T(x)T)SUM"]["J6"] > 1e-3 and abs(d["named"]["(I(x)T)SUM"]["J6"]) < 1e-10
    assert abs(d["named"]["(I(x)T)SUM(I(x)T^2)"]["J6"] - 363 / 160) < 1e-9


def test_11357_degree_by_dimension():
    d = load("w33_pass11357_jarlskog_degree_by_dimension.json")["dims"]
    assert d["2"]["first_not_identically_zero_degree"] == 10 and d["2"]["first_complete_degree"] == 10
    assert d["3"]["first_not_identically_zero_degree"] == 6 and d["3"]["first_complete_degree"] == 6
    assert d["5"]["first_not_identically_zero_degree"] == 4 and d["5"]["first_complete_degree"] == 6
    assert [d[k]["clifford_order_mod_phase"] for k in ("2", "3", "5")] == [24, 216, 3000]


def test_11358_mechanism():
    d = load("w33_pass11358_three_qutrit_mechanism.json")
    t = d["cell_and_solution_count"]
    for k, v in t.items():
        if "nsol=(np.int64(1), 0, 0)" in k:
            assert "bad" not in v
        if "nsol=(0, 0, 0)" in k:
            assert "good" not in v
    assert d["parity_pairing_by_cell"]["non-collinear"].get("11", 0) == 66


def test_11359_rates():
    d = load("w33_pass11359_higher_degree_rates.json")["reps"]
    assert abs(d["(4,4)"]["rho"] - 3 / 8) < 1e-9 and abs(d["(5,5)"]["rho"] - 11 / 24) < 1e-9
    assert abs(d["(9,0)"]["rho"] - 1 / 2) < 1e-9
    assert all(v["basis_rank"] == v["clifford_invariants_character"] for v in d.values())


def test_11360_arithmetic():
    d = load("w33_pass11360_jarlskog_arithmetic.json")
    assert d["distinct_values"] == 9 and d["all_in_cubic_field"]
    rats = sorted(v["a_b_e"][0] for v in d["values"] if v["rational"])
    assert set(rats) == {"9/8", "15/8"}
