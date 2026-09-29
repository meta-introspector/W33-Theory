"""Regression for Pass 11135: FI numbers recorded; the mechanism is flagged as not applying."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11135_fi_scale_no_mechanism as P  # noqa: E402


def test_fi():
    r = P.summarize()
    assert r['n_models'] == 21 and r['one_anomalous_u1'] and r['trq_equals_12_norm2']
    assert r['qTA_abs'] == [3, 4] and not r['fi_mechanism_applies']
    assert 0.14 < r['susy_formula_T_over_MP'][0] < r['susy_formula_T_over_MP'][1] < 0.21
