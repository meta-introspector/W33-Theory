"""Regression for Pass 11092: the 87 non-supersymmetric Z6-I W(3,3) twins are consistent three-generation models,
and every one is tachyonic in its theta sector."""
import json
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11092_nonsusy_z6i_twins_are_tachyonic as P  # noqa: E402

C = json.loads(P.OUT.read_text())
MODELS = {f"{m['index']:02d}": m for m in json.loads(P.MODELS.read_text())["models"]}


def test_every_twin_is_tachyonic_with_a_checkable_witness():
    s = C["summary"]
    assert s["models"] == s["tachyonic_twins"] == s["modular_invariant_twins"] == 87
    for key, r in C["models"].items():
        m = MODELS[key]
        V, W5 = [F(x) for x in m["V"]], [F(x) for x in m["W5"]]
        assert P.modular_invariant(V, P.V_TWIN)
        assert r["tachyonic_left_states"] > 0 and P.check_witness(V, W5, r["witness"])


def test_independent_engine_reproduces_orbifolder_massless_count():
    assert C["summary"]["massless_theta_count_matches_orbifolder"] == 87
    for key in ("02", "17", "85"):
        m = MODELS[key]
        ts = P.theta_sector([F(x) for x in m["V"]], [F(x) for x in m["W5"]])
        assert ts == C["models"][key]["theta_sector"]


def test_right_movers_tachyon_only_in_theta_sectors_and_only_without_susy():
    rm = P.right_mover_table(P.V_TWIN)
    assert {k: v["tachyonic"] for k, v in rm.items() if v["tachyonic"]} == {"1": {"-1/6 scalar": 1}, "5": {"-1/6 scalar": 1}}
    assert all(not v["tachyonic"] for v in P.right_mover_table(P.V_SUSY).values())


def test_font_hernandez_tachyon_free_control():
    fh = P.font_hernandez_control()
    assert fh["modular_invariant"] and fh["tachyon_free"]
    assert fh["roots_kept"] == {"first E8": 42, "second E8": 112}   # SO10 x SU2 (x U1^2) and SO16


def test_twins_are_consistent_three_generation_models():
    s = C["summary"]
    assert s["twin_generations"] == {'{"d": 3, "e": 3, "l": 3, "q": 3, "u": 3}': 87}
    assert s["theta_sector_fermion_chirality_flipped"] == 87
    assert sum(s["twin_higgs_doublet_scalars"].values()) == 87
    assert s["gauge_group_equal_to_parent"] == 87


def test_tachyon_iff_parent_has_massless_oscillator_states():
    assert C["summary"]["tachyonic_iff_parent_has_massless_oscillator_states"] == 3 * 87
