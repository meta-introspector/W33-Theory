"""Regression for Pass 11217: the arrow law A = n - c for qubits (constructions recomputed; class census frozen)."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_constructions_recomputed():
    import w33_pass11217_qubit_arrow_law as Q
    assert Q.check_V(3)["all_split"] and Q.check_V(4)["all_split"]
    for k in (4, 5, 6, 7):
        r = Q.check_W(k)
        assert r["construction_split"] and r["Q_not_in_P"]
        assert r["Pperp_contains_unit"] if k % 2 == 0 else (not r["Pperp_contains_unit"] and r["Pperp_not_in_t2"])
    assert Q.check_hyperbolic([1, 0, 1, 1], 1)["split"]
    assert Q.check_selfdual([1, 1, 1], 2)["split"] and Q.check_selfdual([1, 1, 1, 1, 1], 1)["split"]
    _, S, O = Q.W_piece(2)
    assert Q.exception("W(2)", S, O)["split"]
    assert Q.verifier_control(trials=40)["accepted"] == 0


def test_closed_form_n2_recomputed():
    import w33_pass11208_arrow_universality as U
    import w33_pass11217_qubit_arrow_law as Q
    cls, _ = U.load_classes(U.GAP_SMALL, 2, 2)
    Pl = U.all_planes(2, 2)
    for c in cls:
        S = c["S"] % 2
        A, cc, _ = U.arrow_and_c(S, Pl, 2, 2)
        assert A == 2 - cc and Q.c_formula(S, 2) == cc


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11217_qubit_arrow_law.json").read_text())
    assert d["all_constructions_split"] and d["W_containments_fail"]
    assert d["verifier_control"]["accepted"] == 0
    assert [r["forms"] for r in d["V2k"]] == [4, 8, 16, 32, 64]
    cl = d["classes"]
    assert [cl[n]["classes"] for n in ("2", "3", "4", "5")] == [11, 30, 81, 198]
    assert all(cl[n]["c_formula_exact"] for n in ("2", "3", "4", "5"))
    for n in ("2", "3", "4"):
        assert cl[n]["law"] and cl[n]["criterion_on_every_c0_class"] and cl[n]["split_built_on_every_c0_class"]
