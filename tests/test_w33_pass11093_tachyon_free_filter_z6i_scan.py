"""Regression for Pass 11093: in the Z6-I W(3,3) scan, tachyon-free twins and Standard-Model parents never coincide."""
import json
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11093_tachyon_free_filter_z6i_scan as P  # noqa: E402

C = json.loads(P.OUT.read_text())
S = C["summary"]


def test_47_of_58_base_shifts_are_tachyonic_without_wilson_lines():
    assert S["base_shifts"] == 58 and S["tachyonic_at_wilson_line_free_fixed_point"] == 47
    bases = {b["label"]: b for b in json.loads(P.BASES.read_text())["models"]}
    for lab in ("Z6I_00", "Z6I_14", "Z6I_27", "Z6I_56"):          # recompute a sample exactly
        assert P.base_fixed_point([F(x) for x in bases[lab]["V"]]) == C["base_fixed_points"][lab]


def test_tachyon_freedom_and_the_sm_exclude_each_other():
    t = S["scan_totals"]
    assert (t["tries"], t["tachyon_free"], t["parent_sm"]) == (66000, 9576, 705)
    assert S["sm_parents_among_tachyon_free"] == 0 and t["twin_sm"] == 0
    assert S["expected_overlap_if_independent"] > 70 and S["poisson_probability_of_zero"] < 1e-30


def test_sm_parents_tachyonic_through_exotics_only():
    assert S["sm_parents_checked_for_excited_states"] == 705 and S["sm_parents_without_excited_theta_state"] == 0
    assert set(S["sm_parent_excited_theta_labels"]) <= {"v", "w", "x", "bv", "bw", "bx"}
