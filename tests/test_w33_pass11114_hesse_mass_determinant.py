"""Regression for Pass 11114: Hesse-pencil mass determinant, singular at rho; sigma2/sigma1 >= |d/s|."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11114_hesse_mass_determinant as P  # noqa: E402


def test_hesse_structure_at_rho():
    lam = P.hesse_lambda(P.Y.RHO)
    assert abs(lam - P.W ** 2) < 1e-12 and abs(lam ** 3 - 1) < 1e-12
    assert P.determinant_identity(P.Y.RHO) < 1e-12
    al = P.alignments(P.Y.RHO)
    assert al['(1,0,0)'] == [1.0, 0.5, 0.5] and al['(1,w^0,w^1)'][2] < 1e-9


def test_alignment_bound():
    for T in (P.Y.RHO, 1j):
        m, ds = P.min_ratio(T, starts=15)
        assert abs(m - ds) < 1e-6
