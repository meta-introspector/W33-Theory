"""Regression for Pass 11088: baryonic and mixed hidden condensates in Z3xZ3 never clean the spectrum."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
C = json.loads((ROOT / "data" / "w33_pass11088_baryonic_hidden_condensates.json").read_text())
M = C["models"]
C3, C4, C1 = ("z3z3|Z3Z3_0001_c3__SM_20260934_739", "z3z3|Z3Z3_0001_c4__SM_20260935_3729",
              "z3z3|Z3Z3_0001_c1__SM_20260932_2822")


def test_every_viable_baryonic_vacuum_keeps_massless_fractional_charges():
    for name, r in M.items():
        for c in r["candidates"]:
            assert c["parity_composites_even"] == "viable"
            assert c["light_fractional"] >= c["light_fractional_protected_by_symmetry"] > 0
    assert M[C3]["min_light_fractional_viable"] == 148 and M[C4]["min_light_fractional_viable"] == 142


def test_census_counts():
    assert (M[C3]["census"]["fi_cancelling_rays"], M[C3]["census"]["realizable_curing_rays"], M[C3]["with_baryonic"]) == (3069, 78, 39)
    assert (M[C4]["census"]["fi_cancelling_rays"], M[C4]["census"]["realizable_curing_rays"], M[C4]["with_baryonic"]) == (4417, 39, 24)
    assert (M[C1]["census"]["fi_cancelling_rays"], M[C1]["census"]["realizable_curing_rays"], M[C1]["with_baryonic"]) == (68, 4, 0)


def test_su4_is_fully_broken_and_only_antibaryons_condense():
    for name in (C3, C4):
        for c in M[name]["candidates"]:
            assert 0 in c["full_breaking_slots"] or 1 in c["full_breaking_slots"]
            assert "antibaryon" in c["kinds"] and "baryon" not in c["kinds"]
    assert M[C1]["kinds"] == {"su2pair": 277}   # SU(2): every invariant is a pair, nothing baryonic exists
