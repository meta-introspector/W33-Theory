import importlib.util,json,math
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11042_instrument_dependent_process_memory.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11042_instrument_dependent_process_memory.json").read_text())
def test_replay(): assert P.payload()==C
def test_causal_combs():
    assert C["checks"]["classical_causal_trace"]
    assert C["checks"]["coherent_causal_trace"]
def test_instrument_dependence():
    z=C["instrument_dependence"]["Z_measure_and_reprepare"]["conditional_mutual_information_bits"]
    x=C["instrument_dependence"]["Fourier_measure_and_reprepare"]["conditional_mutual_information_bits"]
    assert z==0.0 and abs(x-math.log2(3))<1e-12
