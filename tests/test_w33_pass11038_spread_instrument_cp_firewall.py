import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11038_spread_instrument_cp_firewall.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11038_spread_instrument_cp_firewall.json").read_text())
def test_replay(): assert P.payload()==C
def test_spreads():
    assert C["spread_census"]["spreads"]==36
    assert C["spread_census"]["Bell_containing_spreads"]==9
def test_cp_firewall():
    assert C["complete_positivity"]["CP_projective_points"]==4
    assert C["complete_positivity"]["non_CP_projective_points"]==36
