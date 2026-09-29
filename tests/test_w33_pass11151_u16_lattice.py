"""Regression for Pass 11151: the zero-radius vector lattice is in the genus of Gamma_v^(0) with root system A15."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11151_u16_lattice as P  # noqa: E402


def test_lattice():
    r = P.summarize()
    assert r['n_models'] == 21 and r['identified_in_all'] and r['discriminant_form_hyperbolic_checked']
