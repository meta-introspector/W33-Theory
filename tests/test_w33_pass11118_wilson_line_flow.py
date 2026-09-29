"""Regression for Pass 11118: model-57 Wilson-line gradient flow and its separatrix."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11118_wilson_line_flow as P  # noqa: E402


def test_flow():
    r = P.summarize()
    assert r['antisymmetric_part_sign'] == {'negative': 15, 'positive': 0} and not r['minimum_inside']
    diag = [f for f in r['flows'] if f['start'][0] == f['start'][1]]
    assert len(diag) == 3 and all(f['first_boundary'].startswith('torus1') for f in diag)
    for top, seq in r['separatrix'].items():
        assert seq[0][1].startswith('torus1') and seq[-1][1].startswith('torus2')
