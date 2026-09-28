"""Regression for Pass 11101: a localized light Higgs gives a single heavy top in the 31 twisted-up models."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11101_localized_higgs_heavy_top as P  # noqa: E402

C = json.loads(P.OUT.read_text())


def test_equilateral_fixed_points():
    g = P.su3_fixed_point_distances()
    assert g["all_fixed_by_theta"] and len(set(g["distances"].values())) == 1


def test_single_heavy_top_and_survivors():
    assert C["classes"]["up"] == {"degenerate": 73, "single_heavy": 31}
    assert C["single_heavy_up_equals_twisted_up"] and C["charm_up_degenerate_in_every_single_heavy_model"]
    assert C["best_models"] == ["2", "10", "13", "14", "15", "35", "53", "57", "69", "77", "78", "102"]
    o = json.loads((ROOT / "data" / "w33_pass11097_orders_33.json").read_text())
    assert all(o[m]["T_star"] == 3 for m in C["best_models"])
