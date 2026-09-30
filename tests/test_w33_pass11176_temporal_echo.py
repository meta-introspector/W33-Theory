"""Regression for Pass 11176: all-or-nothing temporal entanglement through Clifford ticks; amnesia and echo."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11176_temporal_echo as P  # noqa: E402


def test_echo():
    K = P.C.clock()
    assert [round(P.tneg(np.eye(27), q), 9) for q in range(3)] == [1.0, 1.0, 1.0]
    assert [round(P.tneg(K, q), 9) for q in range(3)] == [0.0, 1.0, 0.0]
    N = np.ones((3, 3), int)
    T = P.kick_unitary(N) @ K @ P.kick_unitary(N)
    assert P.order(T) == 9
    ech = P.echo(T, 9)
    assert all(abs(x) < 1e-9 for row in ech[:8] for x in row) and all(abs(x - 1) < 1e-9 for x in ech[8])
    d = json.loads((ROOT / "data" / "w33_pass11176_temporal_echo.json").read_text())
    assert d['all_or_nothing_agrees'] and d['all_or_nothing_checked'] == 120
    assert abs(d['cglmp_through_clock'][0] - 3.06062) < 1e-4 and abs(d['cglmp_through_clock'][1] - 3.16281) < 1e-4
    assert abs(d['cglmp_through_perfect_tick']) < 1e-6
