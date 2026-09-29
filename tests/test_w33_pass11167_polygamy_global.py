"""Regression for Pass 11167: 4/sqrt15 exact on its symmetry class; rigorous bound 25/18; dual see-saw never exceeds 4/sqrt15."""
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11167_polygamy_global as P  # noqa: E402


def test_symmetry_class_and_bound():
    ex = P.symmetry_class_exact()
    assert abs(ex['max_value'] - 4 / math.sqrt(15)) < 1e-12 and ex['boundary_max'] <= 1 + 1e-12
    assert ['2/15', '2/15'] in [c[:2] for c in ex['critical_points']]
    B, spec = P.global_bound(n=300)
    assert abs(B - 25 / 18) < 1e-9
    import numpy as np
    for s, t in [(0.1, 0.05), (0.2, 0.1), (0.03, 0.3)]:
        al = 1 - 2 * s - 2 * t
        nab, nac = P.negs(P.family_state(s, t))
        assert abs(nab - (math.sqrt(t * t + 4 * al * s) - t + s)) < 1e-12
        assert abs(nac - (math.sqrt(s * s + 4 * al * t) - s + t)) < 1e-12
    rng = np.random.default_rng(1)
    vals = [P.seesaw(rng) for _ in range(40)]
    assert max(vals) < 4 / math.sqrt(15) + 1e-9
