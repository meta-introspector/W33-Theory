import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11079_bell_chart_history_boundary_design.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11079_bell_chart_history_boundary_design.json").read_text())
def test_replay(): assert P.payload()==C
def test_designs(): assert C["boundary_design"]["pair_intersection"]==4 and C["history_design"]["pair_intersection"]==18
def test_simplex(): assert C["simplex"]["dimension"]==39 and C["checks"]["centered_antipodes"]
