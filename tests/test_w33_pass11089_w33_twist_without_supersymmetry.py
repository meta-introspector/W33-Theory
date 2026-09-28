"""Regression for Pass 11089: no non-SUSY prime-order W(3,3) orbifold; which families have twins."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
C = json.loads((ROOT / "data" / "w33_pass11089_w33_twist_without_supersymmetry.json").read_text())


def test_order_three_parity_theorem():
    assert C["nine_V_squared_mod_2_over_order3_shifts"] == [0]          # 9V^2 always even
    assert C["nine_vprime_squared_mod_2_order3_twins"] == {"Z3": [1], "Z3xZ3": [1]}   # 9v'^2 always odd


def test_twin_table():
    T = C["twins"]
    assert all(v["n_susy_original"] > 0 for v in T.values())
    assert not T["Z6-II"]["non_susy_twin_exists"]
    for fam in ("Z3", "Z3xZ3"):
        assert T[fam]["non_susy_twin_exists"] and not T[fam]["twin"]["same_shifts"]
    for fam in ("Z6-I", "Z12-I", "Z2xZ6-I", "Z2xZ6-II", "Z3xZ6", "Z6xZ6"):
        assert T[fam]["non_susy_twin_exists"] and T[fam]["twin"]["same_shifts"]


def test_orbifolder_evidence_recorded():
    o = C["orbifolder"]
    assert "rejected" in o["Z3 W(3,3) model (A8 + (1/3,1/3,0^6)) with twist (1/3,1/3,1/3)"]
    assert o["Z6-I W(3,3) Standard Models, twins (twist (1/6,1/6,2/3))"].startswith("87/87 load, N = 0")
