import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11065_split_nonsplit_order81_bridge.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11065_split_nonsplit_order81_bridge.json").read_text())
def test_replay(): assert P.payload()==C
def test_classes(): assert C["K81"]["nilpotency_class"]==2 and C["U81"]["nilpotency_class"]==3
def test_quotient(): assert C["U81"]["quotient_by_center"]["is_H27"]
