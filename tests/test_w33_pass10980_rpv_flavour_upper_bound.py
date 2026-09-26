"""Regression for Pass 10980: flavour structure does not hide the regenerated RPV in the Z6-I D-flat vacua."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
C = json.loads((ROOT / "data" / "w33_pass10980_rpv_flavour_upper_bound.json").read_text())


def test_every_vacuum_fails_independently_of_proton_decay():
    s = C["summary"]
    assert s["models"] == 23
    assert s["no_light_Hu"] + s["five_light_d"] == 23
    assert s["no_light_Hu"] == 7 and s["five_light_d"] == 16


def test_no_flavour_suppression_below_the_bound():
    s = C["summary"]
    assert s["below_1e-10"] == 0
    assert s["min_P"] > 1e-4
