"""Regression for Pass 11095: the non-SUSY orbifolder cross-check and the tachyon-free SO(16)xSO(16) A8 Standard Models."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11095_so16xso16_a8_standard_models as P  # noqa: E402

C = json.loads(P.OUT.read_text())
S = C["summary"]


def test_independent_code_agrees_on_every_sample():
    assert S["crosscheck_models"] == 19 and S["crosscheck_identical"] == 19
    assert "(10,1,1)" in S["tachyon_control"]


def test_witten_shift_enumeration():
    W = S["witten_shifts"]
    assert W["half_norm4_vectors"] == 2160 and W["admissible_V0"] == 793508
    assert sum(W["gauge_classes"].values()) == 793508 and len(W["gauge_classes"]) == 6
    assert all("A8" not in g for g in W["gauge_classes"])            # SU(9) never survives the Witten shift


def test_every_sm_like_model_verified_by_both_engines():
    n = S["sm_like_models"]
    assert n >= 80
    assert S["their_sm_and_tachyon_free"] == n and S["our_tachyonic_fields_at_minus_half"] == 0
    assert S["verified_three_generations"] == n
    assert S["verified_no_net_fractional_colourless"] == n and S["verified_with_higgs_scalar"] == n
    for cls in S["scan"].values():
        assert cls["tachyon_free"] == cls["ok"]                        # every draw tachyon-free


def test_the_open_problem_is_recorded():
    # every model carries vector-like fractionally charged fermions: masses for them are NOT shown
    assert S["models_without_any_fractional_fermion"] == 0 and S["vectorlike_fractional_fermion_states_min"] == 88 and S["vectorlike_fractional_fermion_states_max"] == 230
