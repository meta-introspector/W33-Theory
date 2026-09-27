"""Regression for Pass 11076: the flat-or-massive pattern is an unbroken R-symmetry, not a law."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
C = json.loads((ROOT / "data" / "w33_pass11076_flat_or_massive_is_an_unbroken_r_symmetry.json").read_text())


def test_table():
    T = C["table"]
    assert T["Z3xZ3 minimal"] == {"zero (unbroken symmetry) | massless (unbroken symmetry)": 36}
    assert T["Z3xZ3 extended"] == {"nonzero | massive": 28}
    assert T["Z12-I Z2^R"]["zero (unbroken symmetry) | massive"] == 2          # W forbidden, exotics massive
    assert sum(v for k, v in T["Z12-I weakest"].items() if k.startswith("nonzero | massless")) == 11  # W != 0, massless
    assert sum(T["Z12-I Z2^R"].values()) == sum(T["Z12-I weakest"].values()) == 14


def test_counterexample_is_the_1063_candidate():
    v = C["vacua"]["Z12-I Z2^R"]["Z12-I|Z12I_1063__SM_20260926_303"]
    assert v["W"] == "zero (unbroken symmetry)" and v["exotic"] == "massive"
