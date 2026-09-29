"""Regression for Pass 11142: no two-doublet light Higgs splits charm from up."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11142_higgs_mixing_up_sector as P  # noqa: E402


def test_mixing():
    r = P.summarize()
    assert r['mixtures_tested'] == 576 and r['hierarchical_hits'] == 0
    assert r['exponent_patterns'] == ['[0.0, 0.0, 0.0]', '[0.0, 1.0, 1.0]']
