import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11039_bell_relative_tomography_rank.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11039_bell_relative_tomography_rank.json").read_text())
def test_replay(): assert P.payload()==C
def test_rank():
    assert C["incidence"]["rank_Q"]==21
    assert C["incidence"]["right_kernel_dimension"]==6
def test_sector_growth(): assert C["sector_growth"]["successive_rank_increments"]==[9,6,4,2]
