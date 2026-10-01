"""Regression for Pass 11219: the arrow law for Gaussian dynamics (fully recomputed, seconds)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_recomputed():
    import w33_pass11219_gaussian_arrow as G
    rng = np.random.default_rng(5)
    for e, h, l in [(1, 1, 0), (0, 0, 1), (1, 0, 2)]:
        r = G.build_generic(e, h, l, rng)
        assert r["split"] and r["total_rank"] == 2 * l == r["n"] - r["c"]
    ns = G.nonsemisimple_cases(rng)
    assert all(v["found"] and v["split"] and v["total_rank"] == 2 for v in ns.values())
    op = G.krein_family(-1, [0.05, 0.099, 0.101, 0.3])
    assert [r["A"] for r in op] == [0, 0, 2, 2]
    assert all(r["A"] == 0 for r in G.krein_family(+1, [0.05, 0.3, 0.5]))
    cp = G.collision_points()
    assert cp["opposite signature, eps = |w1-w2|/2"]["A"] == 2 and cp["equal signature, w1 = w2"]["A"] == 0


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11219_gaussian_arrow.json").read_text())
    assert d["generic_law"] and d["nonsemisimple_law"] and len(d["generic"]) == 33
    assert d["krein_arrow_switches_on_at_collision"] and d["krein_same_never"] and d["growth_matches_sqrt"]
