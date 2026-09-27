import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11073_outer_similitude_tick_reversal.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11073_outer_similitude_tick_reversal.json").read_text())
def test_replay(): assert P.payload()==C
def test_outer(): assert C["checks"]["S_is_minus_symplectic_similitude"] and C["checks"]["highest_root_reversed"]
def test_firewall(): assert C["orientation_firewall"]["consequence"]=="tetracode/chamber sheet sign and cubic-tick chirality are distinct Z2 data"
