"""Regression for Pass 11244: char-2 unipotent identities and the conditional qubit limit."""
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_identities_small():
    import w33_pass11244_hesselink_qubit_limit as H
    v = H.verify(maxm=4)
    assert v["identity1_all"] and v["identity2_all"] and v["identity1_covers_all_symplectic_partitions"]
    W, Fw = H.F_series(8)
    assert float(Fw[7] / W[7]) - 0.2645184747 < 1e-9


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11244_hesselink_qubit_limit.json").read_text())
    assert d["identity1_checked"] == 92 and d["identity1_all"] and d["identity2_all"]
    assert d["matches_pass11233_F_m"] and d["steinberg_all_m"]
    lo, hi = (float(x) for x in d["limit_conditional_bracket"])
    assert lo < float(d["limit_conditional"]) < hi and hi - lo < 1e-6
    assert abs(d["mc_F7_z"]) < 3
