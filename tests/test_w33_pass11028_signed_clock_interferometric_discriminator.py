import importlib.util, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p11028",ROOT/"analysis"/"w33_pass11028_signed_clock_interferometric_discriminator.py")
P=importlib.util.module_from_spec(S); S.loader.exec_module(P)
C=json.loads((ROOT/"data"/"w33_pass11028_signed_clock_interferometric_discriminator.json").read_text())

def test_replay():
    assert P.payload()==C

def test_visibility_discriminator():
    x=C["coarse_fibre_probe"]
    assert x["signed_overlaps"]==["1/3","-1/3","1/3","-1/3"]
    assert x["signed_model_visibility"]=="1/3"
    assert x["projective_identity_model_visibility"]=="1"

def test_fringe_extrema():
    x=C["equal_arm_fringe"]
    assert (x["signed_Pmax"],x["signed_Pmin"])==("2/3","1/3")
    assert (x["projective_Pmax"],x["projective_Pmin"])==("1","0")

def test_network_and_qutrit_crosscheck():
    x=C["central_operation"]
    assert (x["disjoint_mode_swaps"],x["negative_swap_pairs"])==(12,6)
    assert C["checks"]["independent_qutrit_trace_ratio_one_third"]
