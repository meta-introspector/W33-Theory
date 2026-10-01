"""Regression for Pass 11231: circuit entanglement versus spacetime min-cut."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_recomputed_small():
    import w33_pass11231_dynamical_area_law as D
    import w33_pass11193_optimal_perfect_gate_compiler as C
    p = C.fixed_perfect_gate() % 3
    r = D.run_tick("p", D.dressed_perfect(p, 50), 2, N=8, depth=6, seed=4, init="bell")
    assert r["rt_exact"] <= r["pairs"] and r["fraction_exact"] > 0.8
    z = D.run_tick("p", [p], 2, N=8, depth=6, init="z")
    assert max(z["half_chain_entropy"]) == 0


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11231_dynamical_area_law.json").read_text())
    best = {(r["init"], r["tick"]): r["rt_exact"] for r in d["runs"]}
    assert best[("bell", "random perfect gates (dressed p)")] == 548
    assert all(r["rt_exact"] < r["pairs"] for r in d["runs"])
