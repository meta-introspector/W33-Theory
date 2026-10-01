"""Regression for Pass 11221: Choi-state mode entanglement grows at 2 x (coupling rank) nats per unit squeezing."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_recomputed():
    import w33_pass11221_gaussian_choi_entropy as C
    rng = np.random.default_rng(7)
    for e, h, l in [(1, 1, 0), (0, 0, 1), (1, 0, 1)]:
        S, planes = C.build(e, h, l, rng)
        ranks, sl, ok = C.check(S, planes)
        assert ok and abs(sum(sl) / 2 - 2 * l) < 0.05
        ranks, sl, ok = C.check(S, C.random_split(len(S) // 2, rng))
        assert ok and all(abs(s - 4) < 0.05 for s in sl)


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11221_gaussian_choi_entropy.json").read_text())
    assert d["slope_equals_twice_rank"] and d["optimal_half_total_equals_A"] and len(d["cases"]) == 18
