import importlib.util
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/"data"/"w33_20261001_maslov_borel_factor54.json"

def load_module():
    spec=importlib.util.spec_from_file_location(
        "factor54",ROOT/"analysis"/"w33_20261001_maslov_borel_factor54.py")
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod

def test_frozen_certificate():
    d=json.loads(CERT.read_text())
    assert d["status"].startswith("PASS_MASLOV_54")
    assert d["effective_stabilizer_isomorphism"]["bijective"]
    assert d["effective_stabilizer_isomorphism"]["gap_smallgroup_projective"]==[162,10]
    assert d["relations"]["6"]["U_imbalance"]=={"27":1,"81":-1}
    assert d["relations"]["7"]["U_imbalance"]=={"27":-1,"81":1}
    assert [d["relations"][k]["chi_from_U_orbits"] for k in ("6","7")]==[-54,54]

def test_target_borel_orders():
    mod=load_module()
    U,B=mod.target_borel()
    assert len(U)==81 and len(B)==162
