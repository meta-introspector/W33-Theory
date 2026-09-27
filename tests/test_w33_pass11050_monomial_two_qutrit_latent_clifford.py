import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11050_monomial_two_qutrit_latent_clifford.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11050_monomial_two_qutrit_latent_clifford.json").read_text())
def test_replay(): assert P.payload()==C
def test_clifford(): assert all(c["clifford_conjugation"]["verified"] for c in C["split_prime_certificates"])
def test_hardware(): assert C["gate_compiler"]["generic_U9_required"] is False
