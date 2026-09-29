"""Regression for Pass 11141: the Gamma_0(3) elliptic point is a charged tachyon; the cusp at 0 is tachyonic."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11141_no_symmetric_stabilization as P  # noqa: E402


def test_no_symmetric_stabilization():
    r = P.summarize()
    assert r['group'] == 'Gamma_0(3)' and r['elliptic_is_charged_disk_centre']
    assert abs(r['elliptic_pR2'] - 2 / 3) < 1e-12 and abs(r['elliptic_Delta'] + 1 / 6) < 1e-12
    assert r['cusp0_diagonal_tachyonic'] == 15 and len(r['cusp0_diagonal_tachyon_free_at_0p2']) == 6
    assert 0.2 < r['golden_lower_root'] < 0.221 and r['fricke_only_approximate']
