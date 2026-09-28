"""Regression for Pass 11096: twins of Z12-I, Z2xZ6-I, Z3xZ6, Z6xZ6 -- no tachyon-free Standard Model."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11096_twins_of_the_other_w33_families as P  # noqa: E402

C = json.loads(P.OUT.read_text())


def test_pass11089_mixed_condition_correction():
    A = C["modular_invariance"]
    for fam in ("Z2xZ6-I", "Z3xZ6", "Z6xZ6"):
        assert A[fam]["pass11089_twin"]["mixed_even"] is False and A[fam]["used_twin"]["mixed_even"] is True
        assert A[fam]["used_twin"]["n_susy"] == 0 and A[fam]["consistent_single_generator_twins"] > 0
    assert A["Z3xZ3"]["consistent_single_generator_twins"] == 0


def test_right_mover_levels_recomputed():
    assert P.rm_levels(P.ORIGINAL["Z12-I"][0], P.TWIN["Z12-I"]) == C["right_mover_tachyon_levels"]["Z12-I"]


def test_every_family_twin_consistent_and_no_tachyon_free_sm():
    S = C["summary"]
    for fam, n in (("Z12-I", 289), ("Z2xZ6-I", 29), ("Z3xZ6", 5), ("Z6xZ6", 10)):
        assert S[fam]["models"] == S[fam]["n0"] == S[fam]["anomaly_free_up_to_one_gs_u1"] == n
        assert S[fam]["susy_parent_twisted_states_at_minus_1_6"] == 0
    assert S["Z12-I"]["tachyonic"] == 289 and S["Z3xZ6"]["tachyonic"] == 5 and S["Z6xZ6"]["tachyonic"] == 10
    assert S["Z2xZ6-I"]["tachyon_free"] == 8
    assert S["Z2xZ6-I_tachyon_free_with_three_sm_generations_under_parent_hypercharge"] == 0
