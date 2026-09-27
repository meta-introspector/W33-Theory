import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11030_g12_reflection_reconstruction_firewall.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11030_g12_reflection_reconstruction_firewall.json").read_text())
def test_replay(): assert P.payload()==C
def test_reflection_group():
    assert C["reflection_realization"]["class_size"]==12
    assert C["reflection_realization"]["generated_subgroup_order"]==48
def test_molien(): assert C["molien"]["basic_invariant_degrees"]==[6,8]
