import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11083_local_e8_flat_holonomy_frontier.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11083_local_e8_flat_holonomy_frontier.json").read_text())
def test_replay(): assert P.payload()==C
def test_frontier(): assert C["topological_reduction"]["remaining_cycle_edges"]==81 and C["E8_overlap_constraints"]["A2_normalizer_order"]==622080
