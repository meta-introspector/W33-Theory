"""Regression for Pass 11126: tachyons carry no family-torus quantum numbers."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11126_family_modulus_under_condensation as P  # noqa: E402


def test_decoupling():
    r = P.summarize()
    assert r['family_torus_quantum_numbers_zero'] and r['one_complex_pair_each']
    assert all(len(v) == 1 for v in r['tachyon_torus'].values())
