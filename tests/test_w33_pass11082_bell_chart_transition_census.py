import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11082_bell_chart_transition_census.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11082_bell_chart_transition_census.json").read_text())
def test_replay(): assert P.payload()==C
def test_profile(): assert C["transition_partition"]=={"boundary_B_boundary_C":4,"boundary_B_finite_C":9,"finite_B_boundary_C":9,"finite_B_finite_C":18,"sum":40}
