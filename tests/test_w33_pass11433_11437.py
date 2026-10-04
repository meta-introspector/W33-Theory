"""Regression for Passes 11433-11437 (grading lemma and n = 4 exhaustive fixed axis; explicit h6; J8 positivity on
J6's blind spot; depth reduction and exact depths 6-7; four-qutrit cells)."""

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def load(name):
    return json.load(open(ROOT / "data" / name))


def test_11433_grading_and_n4():
    d = load("w33_pass11433_magic_axis_law.json")
    g = d["grading_n2"]
    assert g["Pi commutes with V_M T1"] == g["classes"] == 1296
    assert g["frames with c != 0"] == g["U graded by Pi"] == g["spec(U) omega-invariant"] == 69984
    e = d["fixed_axis_n4_exhaustive"]
    assert e["classes"] == 549 and e["group_order"] == 20056328248320
    assert "FAILS" not in e["verdict_classes"]


def test_11434_h6():
    d = load("w33_pass11434_h6_explicit.json")
    assert d["extended_clifford_order"] == 432 and d["monomial_orbits"] == 74
    assert d["all_nonzero_proportional_to_h6"] and len(d["nonzero_odd_orbit_sums"]) == 8
    assert all(sorted(n["mubs"]) == [0, 0, 0, 1, 1, 2] for n in d["nonzero_odd_orbit_sums"])
    assert abs(d["explicit_h6"]["Delta6_over_h6_squared"] - 349920) < 1e-3
    g = d["generic_component_eigenvectors"]["fraction_min_Delta6_below_1e_8"]
    assert g["spurious"] > 3 * g["haar"]


def test_11434_h6_fast():
    import w33_pass11357_jarlskog_degree_by_dimension as P7
    import w33_pass11419_spurious_family as S
    import w33_pass11434_h6_explicit as H
    Cl = P7.clifford_group(3)
    psi = np.array([0.4, 0.7 - 0.2j, -0.3 + 0.5j])
    assert abs(S.ray_witness(psi, Cl) - 349920 * H.h6(psi, Cl) ** 2) < 1e-12


def test_11435_j8_positive():
    d = load("w33_pass11435_j8_positivity.json")
    p = d["pseudo_reflection_component"]
    assert p["condition_holds_everywhere"] and p["min_Delta8_over_Delta7"] > 2 and p["max_law_rel_err"] < 1e-4
    assert d["generic_component"]["min_J8"] > 1e-6
    assert all(s["mpmath_positive"] for s in d["small_distance"])


def test_11436_reduction():
    d = load("w33_pass11436_depth_reduction.json")
    assert d["K_order"] == 27 and d["coset_reps"] == 8
    e = d["exact"]
    assert e["5"]["reversible_fraction"] == "11003/49152" and e["6"]["reversible_fraction"] == "20327/131072"
    assert d["denominator_law_2^(3k-1)3^[k odd]"]
    assert all(v["min_overlap_reversible"] > 2.9999 and v["max_overlap_violating"] < 2.999 for v in e.values())


def test_11437_cells():
    d = load("w33_pass11437_four_qutrit_cells.json")
    c3 = d["n3_control"]["estimate"]
    assert abs(c3["P_estimate"] - 437 / 3276) < 3 * c3["stderr"]
    e4 = d["n4"]["estimate"]
    assert abs(e4["P_estimate"] - 0.125) < 3 * e4["stderr"]
    s = d["n4"]["sampled_cells"]["same line, other"]
    assert s["bad"] / (s["bad"] + s["good"]) < Fraction(77, 342) - 0.05
