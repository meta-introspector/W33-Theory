import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11076_atlas_nerve_steinberg_chain_map.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11076_atlas_nerve_steinberg_chain_map.json").read_text())
def test_replay(): assert P.payload()==C
def test_h1(): assert C["homology"]["isomorphism"] and C["homology"]["induced_map_rank"]==81
def test_chain(): assert C["checks"]["chain_boundary_identity"] and C["checks"]["triangle_boundaries_killed"]
