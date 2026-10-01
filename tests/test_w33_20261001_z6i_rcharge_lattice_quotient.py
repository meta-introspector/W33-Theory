import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_20261001_z6i_rcharge_lattice_quotient as Q

def test_frozen_z6i_rcharge_audit():
    d=json.loads((ROOT/"data/w33_20261001_z6i_rcharge_lattice_quotient.json").read_text())
    a=d["audit"]
    assert (a["models_total"],a["models_dflat"])==(87,23)
    assert a["vacua"]==a["mu_protected_vacua"]==6_695_116
    assert len(a["mu_protected_models"])==23
    assert a["clean_vacua"]==a["clean_and_F_flat"]==0
    assert a["unique_support_lattices"]==1464
    assert all(d["checks"].values())

def test_fast_engine_replays_one_model_and_slow_choices():
    r=Q.model_exact("Z6-I|Z6I_06__SM_20260917_2",validate=3)
    assert r["vacua"]==r["mu_protected"]==82_012
    assert r["clean"]==r["clean_and_F_flat"]==0
    assert r["unique_support_lattices"]==30
    assert r["validation_choices"]==3
