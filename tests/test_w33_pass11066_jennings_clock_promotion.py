import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11066_jennings_clock_promotion.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11066_jennings_clock_promotion.json").read_text())
def test_replay(): assert P.payload()==C
def test_degrees(): assert C["K81"]["Jennings_degrees"]==[1,1,1,2] and C["U81"]["Jennings_degrees"]==[1,1,2,3]
def test_depth(): assert C["K81"]["last_nonzero_augmentation_power"]==10 and C["U81"]["last_nonzero_augmentation_power"]==14
