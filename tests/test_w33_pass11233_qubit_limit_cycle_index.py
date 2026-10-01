"""Regression for Pass 11233: the char-2 cycle index reproduces every exact E_n[c] and fixes the qubit limit."""
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_exact_part_recomputed():
    import w33_pass11233_qubit_limit_cycle_index as P
    r = P.exact_part(N=5, NU=24)
    assert all(r["reconstructs_pass11220"].values())
    assert all(c["steinberg"] and c["E_2_to_d_formula_ok"] for c in r["unipotent_checks"])
    assert r["unitary_steinberg"]
    assert r["unitary_limit"].startswith("0.2726118656378")
    # q-binomial identity: P(1) = prod 1/(1 - 2^{-1-2i}) equals the Steinberg sum
    assert abs(sum(float(P.w(m)) for m in range(60)) - float(r["P1"])) < 1e-12


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11233_qubit_limit_cycle_index.json").read_text())
    assert Fraction(d["F_m_exact"]["6"]["value"]) < Fraction(d["F_m_exact"]["5"]["value"])
    lo, hi = (float(x) for x in d["limit_bracket"])
    assert lo < 0.6651608 < hi
    assert d["E_n_not_monotone"]
    assert abs(float(d["limit_if_F_m_flat"]) - 0.66516078) < 1e-7
