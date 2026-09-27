import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11057_diagonal_weld_orientation_conjugacy.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11057_diagonal_weld_orientation_conjugacy.json").read_text())
def test_replay(): assert P.payload()==C
def test_generator(): assert C["checks"]["generator_inverse_conjugacy"]
def test_gate(): assert C["checks"]["finite_cayley_inverse_conjugacy"]
