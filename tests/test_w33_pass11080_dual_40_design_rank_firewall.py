import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11080_dual_40_design_rank_firewall.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11080_dual_40_design_rank_firewall.json").read_text())
def test_replay(): assert P.payload()==C
def test_rank_firewall(): assert C["point_side"]["rank_F3"]==11 and C["line_side"]["rank_F3"]==15
def test_noniso(): assert C["firewall"]["designs_isomorphic"] is False
