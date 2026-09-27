import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11035_cubic_cocycle_weyl_w33_descent.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11035_cubic_cocycle_weyl_w33_descent.json").read_text())
def test_replay(): assert P.payload()==C
def test_w33():
    x=C["past_future_doubling"]
    assert (x["projective_points"],x["lagrangian_lines"],x["transverse_lagrangians"])==(40,40,27)
def test_tetracode(): assert C["tetracode_indexing"]["nonzero_words"]==8
