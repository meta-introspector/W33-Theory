"""Regression for Pass 11145: neutrino anarchy refuted; Dirac texture is up-like at tree level."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11145_neutrino_anarchy as P  # noqa: E402


def test_neutrino():
    r = P.summarize()
    assert r['neutrino_dirac_tree_up_like'] and r['condensate_leaves_neutrino_texture'] and r['charged_leptons_anarchic']
    assert r['hypothesis_neutrino_anarchy'] == 'refuted'
