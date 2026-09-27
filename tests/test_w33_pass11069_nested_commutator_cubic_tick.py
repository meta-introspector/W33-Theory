import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11069_nested_commutator_cubic_tick.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11069_nested_commutator_cubic_tick.json").read_text())
def test_replay(): assert P.payload()==C
def test_split(): assert C["scalar_replays"]["0"]["nested_commutator"]==[0,0,0,0]
def test_twist(): assert C["scalar_replays"]["1"]["nested_commutator"]==[0,0,0,2] and C["scalar_replays"]["2"]["nested_commutator"]==[0,0,0,1]
