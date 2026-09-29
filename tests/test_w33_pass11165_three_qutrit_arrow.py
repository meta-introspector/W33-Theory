"""Regression for Pass 11165: a perfect three-qutrit tick secret-shares every input's past among all three outputs."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11165_three_qutrit_arrow as P  # noqa: E402


def test_arrow():
    r = P.summarize()
    assert r['tetracode3_symplectic'] and r['n_perfect'] == 12 and r['n_controls'] == 12
    assert r['perfect_max_info_in_proper_subsets'] == 0 and r['perfect_full_info'] == [2]
    assert r['controls_with_leak'] == 12
