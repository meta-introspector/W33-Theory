"""Regression for Pass 11192: arrow and mereology of the anti-symplectic half of W(E6) (frozen; recomputation ~1 min)."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11192_time_reversal_mereology.json").read_text())
    assert d["outer_coset_projective_elements"] == 25920
    assert d["arrow_distribution"] == {"0": 9180, "2": 16704, "4": 36}
    assert d["never_local"] == 10944 and 10944 * 45 == 19 * 25920
    inv = sorted((c["size"], c["A"], c["local"]) for c in d["involution_classes"])
    assert inv == [(36, 4, 15), (540, 0, 7)]
