import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11060_dark15_triple_graph.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11060_dark15_triple_graph.json").read_text())
def test_replay(): assert P.payload()==C
def test_transverse(): assert C["checks"]["no_pure_sector_dark_vectors"]
def test_projection(): assert C["sector_geometry"]["projection_rank_to_S1"]==C["sector_geometry"]["projection_rank_to_S2"]==C["sector_geometry"]["projection_rank_to_L"]==15
