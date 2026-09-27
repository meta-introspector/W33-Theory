import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11047_phase_weld_exact_54_right_inverse.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11047_phase_weld_exact_54_right_inverse.json").read_text())
def test_replay(): assert P.payload()==C
def test_rank_and_inverse():
    for rows in C["certificates"].values():
        for x in rows:
            assert x["quotient_rank"]==54 and x["right_inverse_verified"]
def test_stable_pivots():
    assert all(C["pivot_stability"].values())
