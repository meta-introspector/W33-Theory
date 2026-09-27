import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11037_history_state_process_memory_firewall.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11037_history_state_process_memory_firewall.json").read_text())
def test_replay(): assert P.payload()==C
def test_history_spectrum(): assert C["relational_history"]["clock_reduced_spectrum"]==["19/27","4/27","4/27"]
def test_memory_witness():
    assert C["causal_break_firewall"]["quantum_store_retrieve_negativity"]>0.999
    assert C["causal_break_firewall"]["classical_measure_prepare_negativity"]==0.0
def test_threshold(): assert C["noise_thresholds"]["entanglement_survives_iff"]=="eta > 1/4"
