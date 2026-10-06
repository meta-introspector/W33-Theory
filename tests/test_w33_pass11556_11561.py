import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def load(name):
    return json.loads((ROOT/"data"/name).read_text())

def test_pass11556_fission():
    r=load("PART_W33_PASS11556_PSP_HAMMING_DIRECTED_FISSION_SCHEME.json")
    assert r["status"]=="PASS_PSP_ORBITAL_FISSION_OF_H33"
    assert r["scheme"]["valencies"]==[1,6,12,4,4]
    assert r["schurian_orbitals"]["order"]==648
    assert r["hamming_fusion"]["outer_coset_size"]==648

def test_pass11557_event_dirac():
    r=load("PART_W33_PASS11557_PIN_EQUIVARIANT_EVENT_DIRAC.json")
    assert r["status"]=="PASS_PIN_EQUIVARIANT_FIRST_ORDER_EVENT_DIRAC"
    assert r["vector_differences"]["cubic_factorization"]=="D = B1 B2 B3"
    assert r["pin_cover"]["lift_order"]==48
    assert r["pin_cover"]["projective_order"]==24
    assert r["vector_differences"]["D_rank"]==8
    assert r["covariance"]["maximum_Q_covariance_residual"]<1e-9

def test_pass11558_albert_embedding():
    r=load("PART_W33_PASS11558_ALBERT_PIN_EVENT_SPINOR_EMBEDDING.json")
    assert r["status"]=="PASS_EVENT_PIN_CARRIER_EMBEDS_FOURFOLD_IN_ALBERT16"
    assert r["albert_clifford"]["generated_Cl4_algebra_dimension"]==16
    assert r["pin_action"]["lift_order"]==48
    assert r["pin_action"]["order8_eighth_root_multiplicities"]=={"1":4,"3":4,"5":4,"7":4}
    assert r["clock_no_go_reconciliation"]["Pass10959_single16_old_phase_lift_no_go"] is True

def test_pass11559_fixed_dimension_refinement():
    r=load("PART_W33_PASS11559_FIXED_DIMENSION_GAUGE_FRAME_REFINEMENT.json")
    assert r["status"]=="PASS_FIXED_DIMENSION_COVARIANT_DIRAC_PRINCIPAL_SYMBOL"
    assert r["vector_difference_expansion"]["first_moments"]==[[8,0,0],[0,8,0],[0,0,8]]
    assert r["frame_principal_symbol"]["exact_symbol_square_check"] is True
    assert abs(r["tower"]["rows"][-1]["N_squared_gap"]-2*3.141592653589793**2)<1e-3

def test_pass11560_coinvariant_torsion():
    r=load("PART_W33_PASS11560_TRANSLATION_COINVARIANT_ARROW_TORSION.json")
    assert r["status"]=="PASS_ARROW_FROM_TRANSLATION_COINVARIANT_3_TORSION"
    assert r["coinvariant_chain_complex"]["d2"]=="3 I4"
    assert r["coinvariant_chain_complex"]["smith_invariants_d2"]==[3,3,3,3]
    assert r["residual_WD3"]["fixed_dimension_in_Hom_H1_F3"]==1

def test_pass11561_native_e6_hesse():
    r=load("PART_W33_PASS11561_NATIVE_E6_HESSE_POLYNOMIAL_LAW.json")
    assert r["status"]=="PASS_NATIVE_E6_HESSE_POLYNOMIAL_LAW_AND_STABILIZER"
    assert r["native_carrier"]["native_cubic_contraction"]=="d_native(T_frame z,T_frame z,T_frame z)=6 z1 z2 z3"
    assert r["char3_Hesse"]["coefficientwise_mod3_verified"] is True
    assert r["char3_Hesse"]["N_third_cross_effect_on_e1_e2_e3"]==1
    assert r["char3_Hesse"]["T_cubed_third_cross_effect_on_e1_e2_e3"]==0
    assert r["frame_normalizer"]["stabilizer_of_span_N_T3_in_WD3"]==6
