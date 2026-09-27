import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11055_cayley_finite_54_compiler.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11055_cayley_finite_54_compiler.json").read_text())
def test_replay(): assert P.payload()==C
def test_finite(): assert C["compiler_consequence"]["finite_gate_covers_all_54_retyped_directions"]
def test_identity(): assert C["theorem"]["rank_for_diagonal_weld"]==66 and C["theorem"]["retyped_quotient_rank_for_diagonal_weld"]==54
