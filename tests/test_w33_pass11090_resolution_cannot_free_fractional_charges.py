"""Regression for Pass 11090: dropping selection rules (resolution) never frees all fractional charges."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
C = json.loads((ROOT / "data" / "w33_pass11090_resolution_cannot_free_fractional_charges.json").read_text())


def test_gauge_only_bound():
    s = C["summary"]
    assert s["condensate"]["vacua"] == 55 and s["singlet"]["vacua"] == 28
    assert s["condensate"]["min_light_fractional"]["gauge only"] >= 18
    assert s["singlet"]["min_light_fractional"]["gauge only"] >= 15
    for v in C["vacua"]:
        lf = v["light_fractional"]
        assert lf["all"] >= lf["no SG"] >= lf["no SG, PG"] >= lf["gauge only"] > 0   # monotone, never zero


def test_stuck_states_live_at_unresolved_fixed_points():
    s = C["summary"]
    for tag in ("condensate", "singlet"):
        assert set(s[tag]["isolated"]) == {"elsewhere"}


def test_charge_class_is_a_space_group_character():
    ch = C["charge_character"]
    assert len(ch) == 13
    assert sum(v["fixed_points"] for v in ch.values()) == 946
    assert sum(v["fractional_fixed_points"] for v in ch.values()) == 516
    for v in ch.values():
        assert v["constant_on_every_fixed_point"] and len(v["space_group_characters"]) == 1
        assert v["sm_singlets_at_fractional_fixed_points"] == 0 and v["untwisted_class"] == 0
