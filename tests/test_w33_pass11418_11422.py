"""Regression for Passes 11418-11422 (J8 completeness tests; anatomy of J6's spurious set; the magic-axis rules;
four qutrits; exact depth 5)."""

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def load(name):
    return json.load(open(ROOT / "data" / name))


def test_11418_j8():
    d = load("w33_pass11418_j8_completeness.json")
    pr = d["pseudo_reflections"]
    assert pr["J6"]["zeros_off_real_type"] > 0 and pr["J8"]["zeros_off_real_type"] == 0
    assert pr["J6"]["second_singular_value_relative"] < 1e-9 < pr["J8"]["second_singular_value_relative"]
    s = d["seeded_from_J6_spurious"]
    assert s["j6_spurious_starts"] > 100 and s["certified_j8_spurious"] == 0
    r = d["ratio_test"]
    for delta in ("0.1", "0.3"):
        assert abs(r[f"J6, dist>={delta}"]["J_there"]) < 1e-11                  # J6: genuine zeros far from reversible
    assert r["J8, dist>=0.3"]["J_there"] > 1e-7                                # J8: min over dist >= 0.3 positive
    a = d["identity_asymptotics"]["J8_over_lam12_Delta6"]
    assert all(row[0] < row[1] < row[2] < 252 and row[2] > 235 for row in a)    # -> 252 (smallest beta is rounding-limited)
    assert all(225 < row[3] < 252 for row in a)


def test_11419_anatomy():
    d = load("w33_pass11419_spurious_family.json")
    assert d["factorisation"]["singular_values_relative"][1] < 1e-10
    assert d["factorisation"]["max_rel_dev_from_abs_1_minus_e_ib_pow12"] < 1e-7
    assert d["projector_identity"]["max_abs_dev_of_J6U_over_abs_lam12_minus_J6P_relative_to_max_J6P"] < 1e-5
    odd = {int(k): v["tau_odd"] for k, v in d["odd_ray_invariants"].items()}
    assert odd == {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 1, 7: 2}
    m = {int(k): v for k, v in d["moment_differences_max"].items()}
    assert all(m[k] < 1e-12 for k in range(1, 6)) and m[6] > 1e-6
    p = d["populations"]
    assert p["near_degenerate"] == p["snapped_to_pseudo_reflection_component"]["reached_g_zero"] > 0
    assert d["g_zero_set"]["zeros_off_real_type"] > 50 and all(h["rank"] == 1 for h in d["g_zero_set"]["hessian_at_zeros"])
    assert d["S3_control"]["max_abs_imag_S3"] < 1e-12


def test_11419_identity_fast():
    import w33_pass11357_jarlskog_degree_by_dimension as P7
    import w33_pass11419_spurious_family as S
    Cl = P7.clifford_group(3)
    psi = np.array([1.0, 0.3 + 0.7j, -0.4 + 0.2j])
    psi /= np.linalg.norm(psi)
    lam = np.exp(2.2j) - 1
    U = np.eye(3) + lam * np.outer(psi, psi.conj())
    assert abs(P7.J(U, 3, Cl) - abs(lam) ** 12 * S.ray_witness(psi, Cl)) < 1e-9


def test_11420_rules():
    d = load("w33_pass11420_magic_axis_rules.json")
    lc = d["lemma_checks_n2"]
    assert lc["L1_checked_classes"] == 1296 and lc["L2_frames_with_U3_clifford"] == 4 * 648
    f2 = d["conjecture_F_n2"]
    assert f2["Mz1 = z1"] == {"bad == {omega(a, z1+Mz1) != 0}": 624, "bad contains it": 24}
    assert f2["same line, M^2 z1 = z1"] == {"bad == {omega(a, z1+Mz1) != 0}": 648}
    f3 = d["conjecture_F_n3"]
    assert not any("FAILS" in k for part in f3.values() for k in part)
    assert sum(f3["slow_classes_by_frame_orbits"].values()) > 0
    k = d["L3_kernel_condition"]
    assert k["n2_classes"]["dim ker(M+I) = 1"] == 414 and sum(k["n2_classes"].values()) == 648


def test_11421_four_qutrits():
    d = load("w33_pass11421_four_qutrits.json")
    assert all(v["classes_checked"] == v["identical_good_frame_sets"] for v in d["validation"])
    t = d["targeted_rule_cells"]
    assert t["fixed"].get("good", 0) == 0 and t["fixed"]["bad"] > 50
    assert t["reversed"].get("bad", 0) == 0 and t["reversed"]["good"] > 50
    assert t["line, M^2z1=z1"].get("good", 0) == 0 and t["line, M^2z1=-z1"].get("bad", 0) == 0
    u = d["uniform"]
    assert u["decided"] == 3000 and abs(u["bad_class_fraction"] - 0.125) < 3 * u["stderr"]


def test_11422_depth5():
    d = load("w33_pass11422_depth5_exact.json")
    e = d["exact"]
    assert e["4"]["violating"] == 2077650 and e["5"]["words"] == 71663616
    assert all(v["min_overlap_reversible"] > 2.9999 and v["max_overlap_violating"] < 2.999 for v in e.values())
    assert Fraction(e["5"]["reversible_fraction"]) == Fraction(e["5"]["reversible"], 71663616)
