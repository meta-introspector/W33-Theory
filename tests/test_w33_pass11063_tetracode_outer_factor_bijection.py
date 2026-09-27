import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11063_tetracode_outer_factor_bijection.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11063_tetracode_outer_factor_bijection.json").read_text())
def test_replay(): assert P.payload()==C
def test_bijection(): assert len(C["factor_word_bijection"])==4 and C["tetracode"]["nonzero_words"]==8
def test_central(): assert C["central_factors_outside_tetracode8"]==[8,9]
