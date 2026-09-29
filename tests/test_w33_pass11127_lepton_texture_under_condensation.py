"""Regression for Pass 11127: light-Higgs textures unchanged by the condensate; m_mu = m_e reflection symmetry."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11127_lepton_texture_under_condensation as P  # noqa: E402


def test_textures():
    r = P.summarize()
    assert r['light_exponents_unchanged_0_1_1'] and r['light_lepton_reflection_symmetric']
    assert r['exponent_changes_only_heavy_doublets']
    assert r['new_entries_total'] == {'down': 324, 'lepton': 405}
