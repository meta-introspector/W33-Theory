"""Regression for Passes 11234-11235 (frozen certificates plus a fast recomputation on Sp(4,2))."""
import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
QUANTUM = (1 + 2 * np.cos(2 * np.pi / 9)) / 3


def load(name):
    return json.loads((ROOT / "data" / name).read_text())


def test_11234_frozen():
    d = load("w33_pass11234_qubit_bruteforce_check.json")
    assert d["n2"]["E_c"] == "29/40" and d["n2"]["arrow_law_violations"] == 0 and d["n2"]["arrow_law_checked"] == 720
    assert d["n3"]["E_c"] == "5981/8960" and d["n3"]["order"] == 1451520 and d["n3"]["arrow_law_violations"] == 0
    assert "2" not in d["n3"]["c_distribution"]


def test_11234_sp42_recomputed():
    import w33_pass11234_qubit_bruteforce_check as Q
    G = Q.enumerate_sp(2)
    assert len(G) == 720
    dist = Counter(Q.protected(M, 2) for M in G)
    assert Fraction(sum(c * k for c, k in dist.items()), 720) == Fraction(29, 40)


def test_11235_spectrum():
    d = load("w33_pass11235_two_qutrit_t_spectrum.json")
    assert abs(d["control_TT_SUM"] - QUANTUM) < 1e-12 and abs(d["control_SUM"] - 1) < 1e-12
    by = d["by_depth"]
    for depth in ("1", "2"):
        assert max(float(k) for k in by[depth]["fidelities"]) <= QUANTUM + 1e-9
    assert d["max_violator_fidelity"] > QUANTUM + 0.05 and d["milder_than_quantum"] >= 1
    assert [by[str(k)]["violating"] for k in range(1, 5)] == sorted(by[str(k)]["violating"] for k in range(1, 5))
