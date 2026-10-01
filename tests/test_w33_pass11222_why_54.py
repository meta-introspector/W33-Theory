"""Regression for Pass 11222: chi as a signed orbit count of the relation's symmetry group."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11222_why_54.json").read_text())
    assert d["chi_is_54Q_on_256"] and d["divisible_by_27_on_256"]
    o = d["orbitals"]
    for v in o.values():
        if v["size"] == 256:
            assert v["stabilizer_order"] == 324 and v["orbit_sizes"] == [27, 54, 81, 162]
            s = sum(int(k) * n for k, n in v["unbalanced_by_size"].items())
            assert s == v["chi"] and abs(v["chi"]) == 54
        else:
            assert v["stabilizer_order"] == 12 and abs(v["chi"]) == 18


def test_recompute_one():
    import sys
    import numpy as np
    sys.path.insert(0, str(ROOT / "analysis"))
    import w33_pass11222_why_54 as W
    reps = json.loads((ROOT / "data" / "w33_pass11212_one_chirality.json").read_text())["representatives"]
    k = next(k for k, v in reps.items() if v["size"] == 6912 and v["tau_image"] != int(k))
    a = W.analyse(np.array(reps[k]["B"], np.int64))
    assert a["stabilizer_order"] == 12 and abs(a["chi"]) == 18
