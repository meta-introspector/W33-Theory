import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11034_independent_flavour_symmetry_search.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11034_independent_flavour_symmetry_search.json").read_text())
def test_replay(): assert P.payload()==C
def test_triplet_no_go():
    assert C["irreducible_triplet_firewall"]["assignments"]==16
    assert C["irreducible_triplet_firewall"]["RPV_free"]==0
def test_det_loophole():
    assert C["survivor_structure"]["count"]==4
    assert C["survivor_structure"]["independent_of_clock_central_parity"] is True
