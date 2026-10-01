"""Regression for Pass 11225: residual-symmetry mixing from every 3-dim irrep of every subgroup of W(E6) and Sp(4,3)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_s4_image_recomputed():
    import w33_pass11225_flavour_from_symmetry as F
    reps = F.parse()
    assert len(reps) == 606
    s4 = next(r for r in reps if r["meta"]["image_struct"] == "S4")
    a = F.analyse_image(F.unitarise(s4["mats"]))
    cols = {tuple(c["column"]) for c in a["columns"] if c["viable"]}
    assert (0.166667, 0.166667, 0.666667) in cols and (0.333333, 0.333333, 0.333333) in cols
    assert not a["full_viable"] and not a["cabibbo_any"]


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11225_flavour_from_symmetry.json").read_text())
    assert d["representations"] == 606 and d["distinct_images"] == 34
    assert not d["any_full_viable"] and not d["any_cabibbo"] and not d["legacy_numerology_realised"]
    cols = {tuple(c) for c in d["viable_columns"]}
    assert cols == {(0.077985, 0.275451, 0.646564), (0.095492, 0.25, 0.654508), (0.166667, 0.166667, 0.666667),
                    (0.276393, 0.361803, 0.361803), (0.333333, 0.333333, 0.333333)}
    clif = np.array([(1 - np.cos(2 * np.pi / 9)) / 3, (1 - np.cos(4 * np.pi / 9)) / 3, (1 - np.cos(8 * np.pi / 9)) / 3])
    assert np.allclose(sorted(clif), [0.077985, 0.275451, 0.646564], atol=1e-6)
    tm1 = next(p for p in d["column_predictions"] if p["column"] == [0.166667, 0.166667, 0.666667])
    assert 0.316 < tm1["sin2_th12"][0] < tm1["sin2_th12"][1] < 0.321
