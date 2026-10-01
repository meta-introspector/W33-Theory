"""Regression for Pass 11220: exact qubit means E_n[c] (n <= 4 recomputed; n = 5, 6 frozen)."""
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_recomputed_small():
    import w33_pass11220_qubit_mean_protection as M
    assert M.exact_mean(2)["E_c"] == "29/40"
    assert M.exact_mean(3)["E_c"] == "5981/8960"


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11220_qubit_mean_protection.json").read_text())
    E = {r["n"]: Fraction(r["E_c"]) for r in d["means"]}
    assert E[6] == Fraction(180847576047477961, 271885345530839040)
    assert E[5] == Fraction(6791923649683, 10212039720960)
    assert abs(d["extrapolated_limit"] - 0.665166) < 1e-5
    six = next(r for r in d["means"] if r["n"] == 6)
    assert "5" not in six["c_distribution"]
