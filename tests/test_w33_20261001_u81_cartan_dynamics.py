import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_20261001_u81_54d_cocycle_transducer as T
import w33_20261001_u81_antiunitary_orientation as A

def load(n):
    return json.loads((ROOT/"data"/n).read_text())

def test_transducer_frozen_and_replay():
    d=load("w33_20261001_u81_54d_cocycle_transducer.json")
    assert len(d["rows"])==8 and all(d["checks"].values())
    for r in d["rows"]:
        assert r["triangles"]==27 and r["vertices_per_triangle"]==3
        assert all((p["raw_rank"],p["quotient_rank"],p["S2_rank"],p["L_rank"])==(54,54,27,27)
                   for p in r["prime_replays"])
    assert T.payload()==d

def test_antiunitary_orientation_replay():
    d=load("w33_20261001_u81_antiunitary_orientation.json")
    assert all(d["checks"].values())
    assert d["u81_extension"]["element_pair_checks"]==13122
    assert A.payload()==d

def test_cartan_metric_and_local_vacuum_frozen():
    d=load("w33_20261001_cartan_kinetic_wall_vacuum.json")
    g=d["cartan_grade_pair_metric"]
    assert g["matrix"]==[[590,20,0],[20,980,0],[0,0,300]]
    assert g["principal_minors"]==[590,577800,173340000]
    assert g["positive_definite"]
    v=d["zero_fit_wall_vacuum"]
    assert v["strict_local_projective_maximum"]
    assert v["cp_zero_tolerance"]
    assert all(d["checks"].values())
