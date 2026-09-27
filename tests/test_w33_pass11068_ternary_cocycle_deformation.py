import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11068_ternary_cocycle_deformation.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11068_ternary_cocycle_deformation.json").read_text())
def test_replay(): assert P.payload()==C
def test_jump(): assert C["invariants"]["0"]["center_order"]==9 and C["invariants"]["1"]["center_order"]==3
def test_orientation(): assert C["orientation_isomorphism"]["verified_pairs"]==6561
