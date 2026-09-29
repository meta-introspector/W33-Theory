"""Regression for Pass 11124: the neutral condensate unlocks lepton/down/nu-Dirac Yukawas, never up or mu."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11124_condensate_unlocked_couplings as P  # noqa: E402


def test_unlocked():
    r = P.summarize()
    assert r['qT_outside_singlet_span'] and r['control']
    assert r['needs_T_coeffs'] == ['-1', '1'] and r['needs_T_targets'] == 102
    assert r['lepton_all_six'] and r['up_and_mu_untouched'] and len(r['down_models']) == 4
    assert r['total_new'] == 1863 and r['new_with_twisted_sm_field'] == 1854 and r['orders'] == [2, 3, 4, 5, 6]
    assert r['separate_winding_total'] == 0 and r['lowered'] == 0
