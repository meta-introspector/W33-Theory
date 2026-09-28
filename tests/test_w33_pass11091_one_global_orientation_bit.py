"""Regression for Pass 11091: exactly one global orientation bit (the multiplier character)."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
C = json.loads((ROOT / "data" / "w33_pass11091_one_global_orientation_bit.json").read_text())


def test_orders_and_uniqueness():
    assert C["orders"] == {"Sp": 51840, "GSp": 103680, "PSp": 25920, "PGSp": 51840}
    assert C["Sp_perfect"] and C["multiplier_is_homomorphism"] and C["scalars_have_multiplier_1"]


def test_codec_bits():
    assert C["deck_multiplier"] == {"1": 1, "T": 1, "S": -1, "TS": -1}
    assert C["chi_equals_mu"] and not C["sigma_is_global_character"]
