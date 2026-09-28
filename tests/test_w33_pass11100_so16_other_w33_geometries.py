"""Regression for Pass 11100: SO(16)xSO(16) with the W(3,3) shifts -- only T6/Z3 gives tachyon-free SM-like models."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11100_so16_other_w33_geometries as P  # noqa: E402

C = json.loads(P.OUT.read_text())
T = C["table"]


def test_only_prime_z3_works():
    assert C["only_z3_gives_tachyon_free_sm"]
    assert T["T6/Z3 (Pass 11095)"]["sm_and_tachyon_free"] == 153 and T["T6/Z3 (Pass 11095)"]["tachyon_free_fraction"] == 1.0
    for k in ("T6/(Z3xZ3)", "T6/Z6-I", "T6/Z12-I"):
        assert T[k]["sm_and_tachyon_free"] == 0 and T[k]["tachyon_free_fraction"] < 0.4


def test_scan_sizes_and_bases():
    ev = json.loads(P.EVID.read_text())
    assert ev["Z6-I"]["bases"] == 58 and ev["Z12-I"]["bases"] == 117 and ev["Z3xZ3"]["distinct"]["bases"] == 2
    assert all(v.get("V0") for v in ev["Z6-I"]["per_base"].values())
