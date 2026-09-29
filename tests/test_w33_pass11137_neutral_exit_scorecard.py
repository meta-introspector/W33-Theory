"""Regression for Pass 11137: the scorecard is consistent with every certificate it reads."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11137_neutral_exit_scorecard as P  # noqa: E402


def test_scorecard():
    r = P.summarize()
    assert r['all_consistent'] and r['works']['bottom_tau_same_order_models'] == 14
