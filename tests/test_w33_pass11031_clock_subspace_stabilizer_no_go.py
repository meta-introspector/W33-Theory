import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11031_clock_subspace_stabilizer_no_go.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11031_clock_subspace_stabilizer_no_go.json").read_text())
def test_replay(): assert P.payload()==C
def test_trivial_stabilizer():
    assert C["stabilizer"]["order"]==1
    assert C["stabilizer"]["nonidentity_joint_span_rank_histogram"]=={"6":47}
def test_full_closure(): assert C["checks"]["full_signed_orbit_span_rank24"]
