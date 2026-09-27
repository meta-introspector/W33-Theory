import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11064_five_full54_slope_pair_routes.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11064_five_full54_slope_pair_routes.json").read_text())
def test_replay(): assert P.payload()==C
def test_routes(): assert C["quotient_rank_pattern"]==[54]*5
def test_raw(): assert C["raw_rank_pattern"]==[66,54,54,72,72]
