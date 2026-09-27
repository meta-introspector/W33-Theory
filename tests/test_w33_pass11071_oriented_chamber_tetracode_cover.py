import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11071_oriented_chamber_tetracode_cover.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11071_oriented_chamber_tetracode_cover.json").read_text())
def test_replay(): assert P.payload()==C
def test_cover(): assert C["homogeneous_cover"]["sheets"]==2 and C["homogeneous_cover"]["oriented_chambers"]==320
def test_tetracode(): assert C["line_tetracode"]["nonzero_words"]==8 and all(len(x["tetracode_words"])==2 for x in C["line_tetracode"]["fibres"])
