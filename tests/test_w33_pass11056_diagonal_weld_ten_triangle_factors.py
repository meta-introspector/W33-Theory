import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11056_diagonal_weld_ten_triangle_factors.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11056_diagonal_weld_ten_triangle_factors.json").read_text())
def test_replay(): assert P.payload()==C
def test_factors(): assert len(C["triangle_factors"])==10 and all(x["triangles"]==27 and x["exact_rank"]==54 for x in C["triangle_factors"])
def test_cayley(): assert C["checks"]["support_is_K81_Cayley"]
