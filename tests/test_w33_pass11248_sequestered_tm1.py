"""Regression for Pass 11248: the sign of eps (phi.chi)^2 selects TM1 vs theta13 = 0."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_selection():
    import w33_pass11248_sequestered_tm1 as M
    import w33_pass11230_tm1_alignment as A
    for eps, want in ((-0.02, "TM1"), (0.02, "theta13 = 0 column")):
        r = M.global_min(eps)
        f, c, af, ac = M.nearest_symmetric(r.x)
        assert A.mixing_outcome(A.G3P, A.G3P, f, c) == want


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11248_sequestered_tm1.json").read_text())
    assert d["eps_negative_selects_TM1"] and d["eps_positive_selects_theta13_0"]
    assert max(d["eps0_misalignment_deg"]) < 0.1
