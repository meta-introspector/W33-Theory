"""Regression for Pass 11178: rank-20 geometry of three-qutrit factorisations; 15 signature classes; perfect = 3456."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11178_three_qutrit_orbitals as P  # noqa: E402


def test_orbitals():
    n, g = P.gap_classes()
    assert n == 20 and len(g) == 15 and sum(sum(v) for v in g.values()) == 110565
    merged = sorted(sorted(v) for v in g.values() if len(v) > 1)
    assert merged == [[256, 256, 2304, 6912, 6912], [6912, 6912]]
    perfect = [v for k, v in g.items() if all(x in (2, 3) for x in k)]
    assert perfect == [[3456]]
    cnt, n_s = P.sample_classes(n=50000, seed=5)
    assert set(cnt) <= set(g)
