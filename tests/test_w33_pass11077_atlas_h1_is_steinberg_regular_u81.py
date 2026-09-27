import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11077_atlas_h1_is_steinberg_regular_u81.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11077_atlas_h1_is_steinberg_regular_u81.json").read_text())
def test_replay(): assert P.payload()==C
def test_module(): assert C["checks"]["Levi_H1_is_Steinberg81"] and C["checks"]["native_F3_free_rank1"]
