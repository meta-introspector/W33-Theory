import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11062_single_foliation_finite54_transducer.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11062_single_foliation_finite54_transducer.json").read_text())
def test_replay(): assert P.payload()==C
def test_plus(): assert [x["factor"] for x in C["perfect_single_factors"]["plus"]]==[7,8]
def test_minus(): assert [x["factor"] for x in C["perfect_single_factors"]["minus"]]==[6,9]
