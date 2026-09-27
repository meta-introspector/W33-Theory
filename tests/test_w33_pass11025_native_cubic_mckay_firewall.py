import importlib.util, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p11025",ROOT/"analysis"/"w33_pass11025_native_cubic_mckay_firewall.py")
P=importlib.util.module_from_spec(S); S.loader.exec_module(P)
C=json.loads((ROOT/"data"/"w33_pass11025_native_cubic_mckay_firewall.json").read_text())

def test_replay():
    assert P.payload()==C

def test_schur_cover_correction():
    x=C["cover_and_tensor_correction"]
    assert x["group"]=="GL2(3) = 2^+S4"
    assert x["nonidentity_involutions"]==13
    assert x["binary_octahedral_identification"] is False
    assert x["classical_SL2_McKay_applicable"] is False

def test_corrected_tensor_quiver_and_surviving_saturation():
    x=C["cover_and_tensor_correction"]
    assert x["spin2a_directed_edges"]==14
    assert x["underlying_undirected_edges"]==11
    assert x["conjugate_quiver_relation"]=="B = A^T"
    assert x["clock_saturation_survives"] is True

def test_invariant_cubic_dimensions():
    x=C["invariant_cubic_space"]
    assert (x["total_dimension"],x["Sym3_Vplus"],x["Vplus_tensor_Sym2_Vminus"])==(71,23,48)
    assert x["Sym3_Vminus"]==0 and x["Vminus_tensor_Sym2_Vplus"]==0

def test_native_parity_pattern():
    x=C["native_E6_restriction"]
    assert x["noncentral_triads"]==32
    assert x["central_eigenbasis_distinct_monomials"]==64
    assert x["monomials_by_minus_count"]=={"0":16,"2":48}
