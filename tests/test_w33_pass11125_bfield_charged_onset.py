"""Regression for Pass 11125: charged tachyons need B and small radius; the neutral exit comes first."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11125_bfield_charged_onset as P  # noqa: E402


def test_bfield_map():
    r = P.summarize()
    assert r['same_charged_region_all_six'] and r['neutral_only_at_zero_B'] and r['charged_onset_at_most_1p0']
    assert len(r['lifted_by_B_at_1p4']) == 4 and r['lifted_are_golden']
    assert all(o == {'0.0': None, '0.1': 0.8, '0.2': 0.9, '0.3': 1.0, '0.4': 1.0, '0.5': 1.0}
               for o in r['highest_charged_y'].values())
