import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11061_hesse_slope_factor_alphabet.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11061_hesse_slope_factor_alphabet.json").read_text())
def test_replay(): assert P.payload()==C
def test_alphabet(): assert C["decomposition"]["outer_factors"]==8 and C["decomposition"]["central_factors"]==2
def test_graph(): assert C["commutation_graph"]["symplectic_zero_edges"]==21 and C["commutation_graph"]["symplectic_nonzero_pairs"]==24
