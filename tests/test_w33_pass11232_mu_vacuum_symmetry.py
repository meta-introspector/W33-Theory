"""Regression for Pass 11232: vacuum U(1)s that forbid mu but keep the top Yukawa, over the 215-model ledger."""
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_span_and_matching():
    import w33_pass11232_mu_vacuum_symmetry as M
    V = M.Span([[Fraction(1), Fraction(2), Fraction(0)], [Fraction(0), Fraction(1), Fraction(1)]])
    assert [Fraction(1), Fraction(3), Fraction(1)] in V and [Fraction(0), Fraction(0), Fraction(1)] not in V
    assert M.structural_rank(range(2), range(2), lambda a, b: a == 0) == 1


def test_one_model_recomputed():
    import w33_pass11232_mu_vacuum_symmetry as M
    import w33_pass10960_matter_even_dflat_closure as P
    L, _ = P.load_ledger()
    name = "Z6-II|Z6II_01__SM_20260917_1388"
    r = M.analyse(name, L[name])
    frozen = json.loads((ROOT / "data" / "w33_pass11232_mu_vacuum_symmetry.json").read_text())["models"][name]
    for k in ("fi_rays", "vacuum_spans", "mu_protected_vacua", "mu_protected_clean"):
        assert r[k] == frozen[k]


def test_frozen_summary():
    d = json.loads((ROOT / "data" / "w33_pass11232_mu_vacuum_symmetry.json").read_text())
    s = d["summary"]
    assert s["models"] == 215
    assert set(s["models_with_clean_protected_vacuum"]) <= set(s["models_with_protected_and_generator"])
    assert set(s["models_with_protected_and_generator"]) <= set(s["models_with_mu_protected_vacuum"])


def test_clean_vacuum_charges():
    d = json.loads((ROOT / "data" / "w33_pass11232_mu_vacuum_symmetry.json").read_text())
    c = d["clean_vacua"]["Z6-II|Z6II_06__SM_20260917_1204"]
    assert c["unbroken_dim"] == 2
    q = {k: Fraction(v) / Fraction(17, 4) for k, v in c["U1prime_charges"].items()}
    assert q["bl_1"] == -7 and q["q_1"] == 6 and q["bu_1"] == 1          # top Yukawa q_1 H_u u^c neutral
    assert all(q["bl_1"] + q[f"l_{j}"] != 0 for j in range(1, 11))        # every mu_1j charged
    assert all(q[s] == 0 for s in c["support"])
