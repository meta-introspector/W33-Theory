import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11052_four_hesse_k81_gauges_and_photonic_lowering.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11052_four_hesse_k81_gauges_and_photonic_lowering.json").read_text())
def test_replay(): assert P.payload()==C
def test_gauges():
    assert len(C["four_compiler_gauges"])==4
    assert C["inflated_Hesse_geometry"]["different_parallel_class_intersection_dimension"]==9
def test_photonics():
    x=C["full_K81_fourier_lowering"]
    assert x["total_F3_operations"]==54
    assert x["nine_tritter_inventory"]["total_resource_waves"]==6
