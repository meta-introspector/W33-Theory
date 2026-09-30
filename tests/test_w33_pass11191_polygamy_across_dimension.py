"""Regression for Pass 11191: family maxima of the normalised N_AB + N_AC for D = 2..7."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_closed_forms():
    import w33_pass11191_polygamy_across_dimension as P
    a2 = (4 - np.sqrt(2)) / 14
    assert abs(P.f_family([np.sqrt(1 - 2 * a2), np.sqrt(a2), np.sqrt(a2)], 2) - (8 * np.sqrt(2) - 4) / 7) < 1e-12
    assert abs(P.f_family([np.sqrt(7 / 15), np.sqrt(2 / 15), np.sqrt(2 / 15)], 3) - 4 / np.sqrt(15)) < 1e-12


def test_frozen_decreasing():
    d = json.loads((ROOT / "data" / "w33_pass11191_polygamy_across_dimension.json").read_text())
    v = [d["family_max"][str(D)]["max_sum"] for D in range(2, 8)]
    assert all(x > y for x, y in zip(v, v[1:])) and all(x > 1 for x in v)
    assert d["full_space_D2_matches_family"]
