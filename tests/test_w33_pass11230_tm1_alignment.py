"""Regression for Pass 11230: TM1 sum rule and the S4 vacuum-alignment scan."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_sum_rule_and_involutions():
    import w33_pass11230_tm1_alignment as M
    assert M.tm1_exact_check()
    s12, cd = M.cos_delta((2 / 3, 1 / 6, 1 / 6), 0.02195, 0.43)
    assert abs(s12 - 0.3184) < 1e-3 and -0.32 < cd < -0.30
    inv = [g for g, w in M.tm1_involution() if np.allclose(w, [2 / 3, 1 / 6, 1 / 6])]
    assert len(inv) == 3
    g = np.array([[1, 0, 0], [0, 0, -1], [0, -1, 0]])
    assert any(np.allclose(g, h) for h in inv)


def test_small_scan():
    import w33_pass11230_tm1_alignment as M
    rng = np.random.default_rng(1)
    r = M.single_scan(True, (2, 4), 8, rng)
    assert all(k.startswith("(1,0,0)") or k.startswith("(1,1,1)") for k in r["global_minimum_types"])


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11230_tm1_alignment.json").read_text())
    assert sum(sum(p["outcomes"].values()) for p in d["two_flavons"]) == 1500
    assert all(p["outcomes"].get("TM1", 0) == 0 for p in d["two_flavons"])
    ren = [s for s in d["single_flavon"] if 6 not in s["degrees"]]
    assert all("(1,1,0)" not in k for s in ren for k in s["global_minimum_types"])
    cliff = next(s for s in d["sum_rules"] if s["column"].startswith("Clifford (2/3)"))
    assert cliff["s23_allowed"][0] > 0.58
