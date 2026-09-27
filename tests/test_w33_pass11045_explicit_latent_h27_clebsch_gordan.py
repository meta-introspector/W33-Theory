import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11045_explicit_latent_h27_clebsch_gordan.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11045_explicit_latent_h27_clebsch_gordan.json").read_text())
def test_replay(): assert P.payload()==C
def test_intertwiner():
    assert all(x["CG_rank"]==27 and x["CG_det_mod_p"] for x in C["split_prime_certificates"])
def test_commutant():
    assert all(x["latent_execution_commute"] for x in C["split_prime_certificates"])
