"""Regression for Pass 11241: full-symmetry mu scan."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_lattice():
    import w33_pass11241_discrete_mu_symmetry as D
    L = D.Lattice([[2, 0], [0, 3]])
    assert [4, 3] in L and [1, 0] not in L


def test_frozen():
    s = json.loads((ROOT / "data" / "w33_pass11241_discrete_mu_symmetry.json").read_text())["summary"]
    assert s["z6i_dflat"] == 23 and s["z6i_exhaustive"] and s["z6i_models_mu_protected"] == []
    assert s["models_clean"] == [] and len(s["models_mu_protected"]) == 104
