"""Regression for Passes 11227-11228 (frozen certificates plus fast one-qutrit recomputations)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
QUANTUM = (1 + 2 * np.cos(2 * np.pi / 9)) / 3


def load(name):
    return json.loads((ROOT / "data" / name).read_text())


def test_11227_two_qutrit_verdicts():
    d = load("w33_pass11227_two_qutrit_cp_mixing.json")
    c = d["circuits"]
    assert d["symplectic_classes"] == 51840
    assert c["(T(x)T) SUM"] == "T-VIOLATING" and c["(T(x)T) SUM (T(x)T) SUM"] == "T-VIOLATING"
    for k in ("T(x)T (control: diagonal cubic only)", "SUM (Clifford control)", "(I(x)T) SUM", "(T(x)I) SUM",
              "SUM (I(x)T) SUM^dag", "SUM (I(x)T) SUM", "(I(x)T) SUM (I(x)T)"):
        assert c[k] == "reversible"
    assert all(abs(v["fidelity"] - QUANTUM) < 1e-9 for v in d["fidelities"].values())


def test_11228_witnesses():
    d = load("w33_pass11228_t_witness.json")
    assert d["J3_max_abs_on_240_words"] < 1e-9 and all(d["J3_symmetries"].values())
    mw = d["moment_witness_vs_exact"]
    assert "reversible|moment!=0" not in mw and mw["violating|moment!=0"] > 0
    assert abs(d["T_X_fidelity"] - QUANTUM) < 1e-12 and d["minimal_violator_fidelities"] == {"0.844029628746": 18}


def test_11228_fidelity_recomputed():
    import w33_pass11213_cubic_t_violation as P
    import w33_pass11228_t_witness as W
    C = P.clifford1()
    assert abs(W.fidelity(P.T1 @ P.X1, C) - QUANTUM) < 1e-12
    assert abs(W.fidelity(P.T1, C) - 1) < 1e-12
    assert abs(W.J3(P.T1 @ P.X1)) < 1e-9
