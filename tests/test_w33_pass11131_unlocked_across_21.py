"""Regression for Pass 11131: never up or mu across the 21; down 14, lepton 17, nu 20; reproduces Pass 11124."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11131_unlocked_across_21 as P  # noqa: E402


def test_across_21():
    r = P.summarize()
    assert r['n_models'] == 21 and r['hidden_check_removes_none'] and r['control_all'] and r['reproduces_11124']
    assert r['models_with_unlocked'] == {'up': 0, 'down': 14, 'lepton': 17, 'nu_dirac': 20, 'mu_HuHd': 0}
    assert r['nothing_unlocked'] == ['A8SM_20260982_47048_TF']
