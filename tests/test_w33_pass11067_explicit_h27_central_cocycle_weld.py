import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11067_explicit_h27_central_cocycle_weld.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11067_explicit_h27_central_cocycle_weld.json").read_text())
def test_replay(): assert P.payload()==C
def test_nontrivial(): assert C["cohomology_certificate"]["coefficient_rank_mod3"]==25 and C["cohomology_certificate"]["augmented_rank_mod3"]==26
def test_bridge(): assert C["same_carrier_bridge"]["underlying_set"]=="H27 x F3, 81 labels"
