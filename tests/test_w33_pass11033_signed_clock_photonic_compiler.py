import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11033_signed_clock_photonic_compiler.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11033_signed_clock_photonic_compiler.json").read_text())
def test_replay(): assert P.payload()==C
def test_network():
    assert C["network"]["arbitrary_routing_swap_count"]==12
    assert C["network"]["nearest_neighbor_adjacent_swap_count"]==148
    assert C["network"]["nearest_neighbor_min_parallel_depth"]==15
def test_noise_margin(): assert C["interferometer"]["guaranteed_separation_epsilon_deg"]>30
