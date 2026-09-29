"""Regression for Pass 11134: the zero-radius tachyon is one complex winding species in all 21."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11134_dual_tachyon_species as P  # noqa: E402


def test_species():
    r = P.summarize()
    assert r['n_models'] == 21 and r['one_complex_species'] and r['three_W_in_lattice']
    assert r['l1sq_distribution'] == {'0.222222': 13, '0.555556': 5, '0.444444': 2, '0.777778': 1}
