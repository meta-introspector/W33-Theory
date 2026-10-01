"""Regression for Pass 11249: AME(10,3) states with a twisted cyclic symmetry are all F9-linear."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_frozen_and_spotcheck():
    import w33_pass11249_ame10_cyclic_search as A
    d = json.loads((ROOT / "data" / "w33_pass11249_ame10_cyclic_search.json").read_text())
    s = {x["twist"]: x for x in d["searches"]}
    assert (s["I"]["lagrangians"], s["I"]["ame_states"], s["I"]["f9_linear"]) == (1600, 48, 48)
    assert s["-I"]["ame_states"] == 0 and s["J"]["ame_states"] == 0
    assert all(len(x["non_f9"]) == 0 for x in d["searches"])
    B = np.array(s["I"]["examples"][0]["basis"])
    assert A.isotropic(B) and A.min_weight(B) == 6
    assert A.f9_linear(B)[1] is True
    S = A.sigma(A.TWISTS["I"])
    assert A.key((B @ S.T) % 3) == A.key(B)
