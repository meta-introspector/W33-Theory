"""Regression for Pass 11119: neutral tori across 491 models; 21 survivors; 8 diagonal-neutral; 6 neutral to Im T = 0.6."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11119_neutral_tori_across_491 as P  # noqa: E402


def test_neutral_tori():
    r = P.summarize()
    assert r['models'] == 491 and r['neutral_only'] == 38 and r['separated'] == 22 and r['all_twisted']
    assert len(r['passing_gauntlet']) == 21 and len(r['passing_neutral_only']) == 8
    assert len(r['diagonal_neutral_only']) == 8 and len(r['diagonal_charged_first']) == 2
    assert len(r['neutral_down_to_0p6_on_axis']) == 6
