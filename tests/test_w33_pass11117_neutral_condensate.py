"""Regression for Pass 11117: neutral condensate breaks one U(1), preserves the SM, no cubic term."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11117_neutral_condensate as P  # noqa: E402


def test_breaking_pattern():
    r = P.summarize()
    for m in ('57', '53'):
        v = r[m]
        assert v['broken_u1_rank'] == 1 and v['hypercharge_zero'] and v['nonabelian_singlet']
        assert v['charged_under_anomalous'] and v['opposite_pair'] and not v['cubic_allowed']
