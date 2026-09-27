import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11072_internal_borel_deck_involution.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11072_internal_borel_deck_involution.json").read_text())
def test_replay(): assert P.payload()==C
def test_tick_firewall(): assert C["checks"]["tetracode_sign_swapped"] and C["checks"]["highest_root_fixed"]
def test_cocycle(): assert C["checks"]["kappa_invariant"]
