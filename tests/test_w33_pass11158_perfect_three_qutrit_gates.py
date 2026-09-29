"""Regression for Pass 11158: column/row determinant laws; perfect three-qutrit gates exist; columns (1,1,-1)."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11158_perfect_three_qutrit_gates as P  # noqa: E402


def test_three_qutrit():
    r = P.summarize()
    assert r['column_law'] == 1.0 and r['row_law'] == 1.0 and r['perfect_equals_single_blocks']
    assert 0.029 < r['perfect_fraction'] < 0.033
    assert list(r['column_patterns']) == ['[(1, 1, 2), (1, 1, 2), (1, 1, 2)]']
    S = np.array(r['example']) % 3
    J = np.zeros((6, 6), int)
    for k in range(3):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, -1
    assert np.array_equal((S.T @ J @ S) % 3, J % 3)
    dets = [[int(round(S[2*i, 2*j] * S[2*i+1, 2*j+1] - S[2*i, 2*j+1] * S[2*i+1, 2*j])) % 3 for j in range(3)] for i in range(3)]
    assert all(sum(dets[i][j] for i in range(3)) % 3 == 1 for j in range(3)) and all(v != 0 for row in dets for v in row)
