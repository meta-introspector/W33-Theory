"""Regression for Pass 11113: the rho-safe / quark-sector table over 491 models is not significant."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11113_rho_tachyon_correlation as P  # noqa: E402


def test_table_and_significance():
    r = P.summarize()
    assert r['models'] == 491
    assert r['table'] == {'twisted|rho_safe=True': 0, 'twisted|rho_safe=False': 152,
                          'untwisted|rho_safe=True': 5, 'untwisted|rho_safe=False': 334}
    assert 0.1 < r['fisher_one_sided'] < 0.2 and not r['significant_at_5pct']
