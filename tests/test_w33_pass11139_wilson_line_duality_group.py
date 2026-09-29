"""Regression for Pass 11139: the Wilson-line duality group is Gamma_0(3); S and Gamma^0(3) broken, Fricke approximate."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11139_wilson_line_duality_group as P  # noqa: E402


def test_group():
    r = P.summarize()
    assert set(r['exact']) == {'T+1', 'T+3', 'G0(3):T/(3T+1)', 'ell3:(-T+1)/(-3T+2)'}
    assert set(r['broken']) == {'S', 'G^0(3):T/(T+1)'} and r['approximate'] == ['Fricke3']
    for m, (k3, k4) in r['deviations']['G0(3):T/(3T+1)'].items():
        assert k4 < k3 and k4 < 1e-8
