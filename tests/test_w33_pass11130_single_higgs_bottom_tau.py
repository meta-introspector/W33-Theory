"""Regression for Pass 11130: the top Higgs's conjugates are the heavy doublets; the condensate lowers m_b/m_t to eps^2."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11130_single_higgs_bottom_tau as P  # noqa: E402


def test_single_higgs():
    r = P.summarize()
    assert r['tree_up_is_678'] and r['conjugates_are_heavy_doublets'] and r['light_generations_degenerate']
    assert len(r['bottom_lowered_3_to_2']) == 4
    assert abs(r['eps_needed']['b_n2'] - 0.116) < 1e-3 and abs(r['eps_needed']['tau_n2'] - 0.151) < 1e-3
