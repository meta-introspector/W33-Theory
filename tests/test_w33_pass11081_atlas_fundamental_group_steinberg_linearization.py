import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11081_atlas_fundamental_group_steinberg_linearization.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11081_atlas_fundamental_group_steinberg_linearization.json").read_text())
def test_replay(): assert P.payload()==C
def test_free(): assert C["graph_certificate"]["free_rank_value"]==81
def test_linearization(): assert C["checks"]["mod3_is_Steinberg81"] and C["checks"]["U81_restriction_regular"]
