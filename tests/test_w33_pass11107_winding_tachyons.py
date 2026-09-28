"""Regression for Pass 11107: explicit winding tachyons of the SO(16)xSO(16) A8 survivors."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11107_winding_tachyons as P  # noqa: E402


def test_enumerator_controls():
    r = P.tachyons('2', [1j, 1j], P.RHO)
    assert abs(r[0]['Delta'] + (3 - 3 ** 0.5) / 6) < 1e-6          # -(3 - sqrt3)/6 at Im T_WL = 1
    assert P.tachyons('2', [2j, 2j], P.RHO) == []                    # tachyon-free at Im T_WL = 2
    assert P.tachyons('2', [5j, 5j], 5j) == []                       # large radius: 10D SO16^2 has no tachyon


def test_radii_and_charges():
    s = P.summarize()
    assert sorted(v for k, v in s['critical_radius_ImT_WL'].items() if k != '77') == [1.512] * 7 + [1.73] * 4
    assert s['models_with_neutral_tachyon'] == ['53', '57'] and len(s['models_all_fractional']) == 10
