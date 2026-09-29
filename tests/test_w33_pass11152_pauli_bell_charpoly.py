"""Regression for Pass 11152: Pauli-only CGLMP Bell operator charpoly = (x^3-6x-2)(x^3-3x+1)^2."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11152_pauli_bell_charpoly as P  # noqa: E402


def test_charpoly():
    r = P.summarize()
    assert r['charpoly_coeffs'] == [1, 0, -12, 0, 45, -6, -57, 18, 6, -2]
    assert r['deviation'] < 1e-40 and r['max_imag'] < 1e-40
    assert abs(r['top_eigenvalue'] - r['cubic_root']) < 1e-12 and abs(r['trig'] - r['cubic_root']) < 1e-12
