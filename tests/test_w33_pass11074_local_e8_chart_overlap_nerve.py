import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11074_local_e8_chart_overlap_nerve.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11074_local_e8_chart_overlap_nerve.json").read_text())
def test_replay(): assert P.payload()==C
def test_nerve(): assert C["nerve"]["f_vector"]==[40,240,160,40] and C["nerve"]["chart_overlap_degree"]==12
def test_h1(): assert C["homology"]["betti"]==[1,81,0,0]
