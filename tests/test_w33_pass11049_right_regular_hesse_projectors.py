import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11049_right_regular_hesse_projectors.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11049_right_regular_hesse_projectors.json").read_text())
def test_replay(): assert P.payload()==C
def test_cross_geometry():
    for c in C["split_prime_certificates"]:
        assert len(c["cross_pairs"])==54
        assert all(r["intersection_dimension"]==1 for r in c["cross_pairs"])
def test_spectrum():
    assert all(r["PQP_spectrum_on_imP"]=={"0":2,"1":1,"1/3":6} for c in C["split_prime_certificates"] for r in c["cross_pairs"])
