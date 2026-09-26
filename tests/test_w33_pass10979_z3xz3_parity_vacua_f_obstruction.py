"""Regression for Pass 10979: T6/Z3 and Z3xZ3 with the W(3,3) twist -- scans, parity census, the 28 extended
vacua, and the exact F-flatness obstruction (exponent vectors independent through degree 11)."""
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
C = json.loads((ROOT / "data" / "w33_pass10979_z3xz3_parity_vacua_f_obstruction.json").read_text())
SPEC = importlib.util.spec_from_file_location("p10979", ROOT / "analysis" / "w33_pass10979_z3xz3_parity_vacua_f_obstruction.py")
P = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(P)


def test_scans():
    s = C["A"]["scan"]
    assert (s["T6/Z3"]["tries"], s["T6/Z3"]["standard_models"]) == (64000, 0)
    assert (s["Z3xZ3"]["tries"], s["Z3xZ3"]["standard_models"]) == (40000, 13)
    assert C["A"]["distinct_spectra"] == 12


def test_census_all_rules_established():
    assert C["B"]["summary"] == {"models": 13, "parity_exists": 4, "closed_farkas": 1, "fi_rays_realizable": 3}


def test_minimal_supports_flat_or_massive():
    viable = 0
    for name, r in C["C"].items():
        viable += r["summary"].get("viable", 0)
        if r["summary"].get("viable"):
            assert r["exotic_d_massless_all"] and r["W_S_zero_all"] and r["max_linear_order"] <= 4
    assert viable == 36


def test_extended_vacua_and_f_obstruction():
    E = C["E"]["summary"]
    assert E["extended_vacua"] == 28 and E["dflat_full"] == 28
    assert E["higgs_massless_all_orders"] == 22 and E["higgs_massive"] == 6
    assert E["min_first_dependent_degree"] == 12 and E["independent_through_14"] == 12
    for vs in C["E"]["vacua"].values():
        for v in vs:
            if "dflat_full" in v:
                assert v["first_dependent_degree"] in (None, 12)
                assert min(int(k) for k in v["monomials_by_degree"]) in (3, 4)


def test_obstruction_lemma_control():
    # non-vacuous control: x^3 and x^6 have proportional exponent vectors -> dependent at degree 6
    assert P.f_obstruction([(3, (3, 0)), (6, (6, 0))]) == 6
    assert P.f_obstruction([(3, (2, 1)), (4, (1, 3))]) is None


def test_recompute_one_vacuum_exactly():
    blob = P.load()
    name = "z3z3|Z3Z3_0001_c4__SM_20260935_3729"
    v = next(x for x in C["E"]["vacua"][name] if x.get("first_dependent_degree") == 12)
    b = name.split("|")[1]
    Z = P.P8.Z6IIR(name, {name: blob["ledger"][name]}, {b: blob["disc"][b]}, {})
    mons = P.Tools(Z).monomials(v["support"], 12)
    assert {str(d) for d, _ in mons} == set(v["monomials_by_degree"]) & {str(d) for d in range(2, 13)}
    assert P.f_obstruction([m for m in mons if m[0] <= 11]) is None
    assert P.f_obstruction(mons) == 12
