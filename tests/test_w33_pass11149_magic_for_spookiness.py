"""Regression for Pass 11149: Pauli-only optimum = root of x^3 - 6x - 2; magic in state or measurements; no threshold."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11149_magic_for_spookiness as P  # noqa: E402


def test_magic():
    r = P.summarize()
    assert abs(r['pauli_measurements_best_state_value'] - r['cubic_x3_6x_2_root']) < 1e-9
    assert abs(r['trig_form'] - r['cubic_x3_6x_2_root']) < 1e-12
    assert r['pauli_best_state_mana'] > 0.3 and abs(r['omega_mana']) < 1e-9
    assert abs(r['omega_with_cglmp_bases_value'] - 2.8729340512) < 1e-8 and r['pauli_basis_vector_mana'] == [0.0]
    curve = r['min_mana_curve']
    ratios = [m / (v - 2) for m, v in curve.values()]
    assert all(0.55 < x < 0.72 for x in ratios)
