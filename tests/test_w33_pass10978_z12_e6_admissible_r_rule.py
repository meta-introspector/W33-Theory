"""Regression for Pass 10978: E6 centralizer, Z12-I verdicts, the 1063 flat-or-light trade, the 10967 audit."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
C = json.loads((ROOT / "data" / "w33_pass10978_z12_e6_admissible_r_rule.json").read_text())


def test_e6_centralizer():
    A = C["A"]
    assert A["W_order"] == 51840 and A["coxeter_order"] == 12
    assert A["centralizer_order"] == 12 and A["centralizer_is_cyclic_c"]
    assert A["minus_one_in_W"] is False
    assert A["eigenvalues_minus_c6"] == [-1.0, -1.0, 1.0, 1.0, 1.0, 1.0]
    assert A["eigenphases_c"] == [1, 4, 5, 7, 8, 11]


def test_all_fourteen_decided():
    models = C["B"]["models"]
    assert len(models) == 14
    verdicts = [m["verdict"] for m in models.values()]
    assert sum(v.startswith("dead") for v in verdicts) == 13
    assert sum(v.startswith("candidate") for v in verdicts) == 1
    assert models["Z12-I|Z12I_1063__SM_20260926_303"]["verdict"].startswith("candidate")


def test_candidate_is_flat_or_light():
    c = C["C"]
    assert tuple(c["higgs_lattice"]) == (3, 4) and tuple(c["exotic_d_lattice"]) == (4, 4)
    assert c["W_S_possible"] is False
    assert c["linear_terms"][0] == [4, "n_36"]
    assert all(tuple(v["higgs_lattice"]) == (4, 4) for v in c["single_absorptions"].values())


def test_audit_of_pass10967():
    for mode in ("orbifolder", "prime_planes", "gamma_corrected"):
        assert C["D"][mode] == {"order3": 0, "order4": 0}
