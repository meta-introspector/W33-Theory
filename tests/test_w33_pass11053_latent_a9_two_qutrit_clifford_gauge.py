import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11053_latent_a9_two_qutrit_clifford_gauge.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11053_latent_a9_two_qutrit_clifford_gauge.json").read_text())
def test_replay(): assert P.payload()==C
def test_exact_gates():
    assert C["exact_gates"]["X_A"].startswith("SUM")
    assert C["exact_gates"]["resource_class"]=="two-qutrit Clifford"
def test_split_primes(): assert all(all(x["SUM_Clifford_conjugation_laws"].values()) for x in C["split_prime_certificates"])
