"""Regression for Pass 11105: Delta(54) = H27 : <-I> is the flavour group of the Z3 rule; one family torus in 104/104."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11105_delta54_flavour_join as P  # noqa: E402


def test_rule_symmetry_is_delta54():
    a = P.part_a()
    assert a["rule_symmetry_group_order"] == 54 and a["rule_group_equals_Delta54"]
    assert a["H27_order"] == 27 and a["heisenberg_ZX_eq_P_XZ"] and a["centre_is_point_group"]
    assert sorted(a["conjugation_image_on_F3sq"].values()) == [27, 27]          # image {+I, -I}
    assert a["clifford648_index"] == 12 and not a["quadratic_phase_is_rule_symmetry"]


def test_geometry_w33_and_schur():
    b = P.part_b()
    assert b == dict(AG33_lines=117, AG23_hesse_lines=12, theta_theta2_untwisted_needs_same_point=True)
    assert P.w33_from_two_qutrit_pauli() == dict(points=40, degree=[12], lam=[2], mu=[4])
    assert P.part_d() == dict(stabiliser_order=18, commutant_dim_on_charm_up=1)


def test_one_family_torus_in_every_model():
    c = json.loads(P.OUT.read_text())["C_spectra"]
    assert c["count_distribution"] == {"1": 104}
    assert sum(c["torus_distribution"].values()) == 104
