import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11041_phase_weld_pairwise12_dressing.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11041_phase_weld_pairwise12_dressing.json").read_text())
def test_replay(): assert P.payload()==C
def test_projection():
    r=C["certificates"]["center_plus_external"]["103"]
    assert r["projection_ranks"]=={"C0":1,"C1":6,"C2":12,"C3":8}
def test_transverse():
    r=C["certificates"]["center_minus_external"]["109"]
    assert r["pure_sector_intersection_dimensions"]=={"C0":0,"C1":0,"C2":0,"C3":0}
