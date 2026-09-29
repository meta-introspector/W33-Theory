"""Regression for Pass 11157: B^3 - 3B + 1 = (B + 1)(I + U + U^2), U = XZ^2 (x) XZ^2."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11157_ninth_roots_explained as P  # noqa: E402


def test_identity():
    r = P.summarize()
    assert r['identity_deviation'] < 1e-30 and r['commutator'] < 1e-30
    assert r['weyl_terms'] == 8 and r['all_moduli_one_over_sqrt3']
    assert len(r['pauli_symmetries']) == 3
