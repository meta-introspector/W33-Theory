import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11051_tensor_induction_fourier_compiler.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11051_tensor_induction_fourier_compiler.json").read_text())
def test_replay(): assert P.payload()==C
def test_blocks():
    assert C["compiler"]["balanced_blocks"]==9
    assert C["compiler"]["matrix_nonzeros"]==81
def test_fourier(): assert all(c["every_component_dephases_to_F3"] for c in C["split_prime_theta_certificates"])
