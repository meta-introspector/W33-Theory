import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def load(n,s): return json.loads((ROOT/"data"/f"PART_W33_PASS{n}_{s}.json").read_text())

def test_11562_cl3_no_go():
    r=load(11562,"CANONICAL_CL3_INSIDE_ALBERT_CL9")
    assert r["status"]=="PASS_CL3_EMBEDDING_CLASSIFIED_BUT_NOT_CANONIC"
    assert r["orbit"]["spin9_plane_stabilizer_dimension"]==18
    assert r["orbit"]["centralizer_of_event_spin3_dimension"]==15
    assert r["clifford"]["generated_Cl3_algebra_dimension"]==8

def test_11563_connection():
    r=load(11563,"DISCRETE_SPIN_CONNECTION")
    assert r["status"]=="PASS_UNIQUE_DISCRETE_TORSION_FREE_METRIC_CONNECTION"
    assert r["exact_determinant"]=="det A(E) = -2 (det E)^3"
    assert r["exact_sample"]["det_frame"]==2

def test_11564_lhs():
    r=load(11564,"LHS_ARROW_ORIGIN")
    assert r["status"]=="PASS_LHS_NO_GROUP_H1_ARROW_AND_RESIDUAL_FIXED_LINE"
    assert r["five_term"]["conclusion_H1_affine_PSp_trivial_coefficients_dimension"]==0
    assert r["direction_module"]["fixed_dimension"]==1
    assert r["direction_module"]["H1_W_M_dimension"]==0

def test_11565_f4():
    r=load(11565,"E6_HESSE_STABILIZER_F4")
    assert r["status"]=="PASS_HESSE_TRACE_CUBE_STABILIZER_IS_F4"
    assert r["stabilizer"]["trace_line_in_E6_dimension"]==52
    assert r["stabilizer"]["hesse_two_plane_stabilizer_dimension"]==52
    assert r["executable_albert"]["rank_of_trace_plus_traceless_trace_images"]==27

def test_11566_axis_wilson():
    r=load(11566,"AXIS_DIRAC_HAMMING_WILSON")
    assert r["status"]=="PASS_AXIS_DIRAC_AND_HAMMING_WILSON_REGULATOR"
    assert r["finite_N3"]["exact_arrow_factorization"]=="D = delta_1 delta_2 delta_3"
    assert r["finite_N3"]["exact_square"]=="Q_axis^2 = I4 tensor (6I-A1)"
    assert r["continuum"]["wilson_zero_set"]=="origin only for r != 0 and zero bare mass"
    assert "12 nodal lines" in r["correction_to_null_edge_refinement"]["global_pathology"]

def test_11567_heat_trace():
    r=load(11567,"WILSON_SPECTRAL_ACTION")
    assert r["status"]=="PASS_WILSON_SPECTRAL_ACTION_SINGLE_SPECIES_LIMIT"
    assert abs(r["N81_t1"]["unregulated_ratio"]-8)<0.05
    assert abs(r["N81_t1"]["wilson_ratio"]-1)<0.03

def test_11568_winding():
    r=load(11568,"WILSON_DIRAC_WINDING")
    assert r["status"]=="PASS_INTEGER_WILSON_DIRAC_WINDING_PHASES"
    assert [x["winding"] for x in r["phases"]]==[0,1,-2,1,0]
    for row in r["direct_integral_samples"]:
        assert abs(row["numerical_winding"]-row["corner_formula"])<2e-6

def test_11569_gauge_no_go():
    r=load(11569,"EXCEPTIONAL_GAUGE_SUBGROUP_SEARCH")
    assert r["status"]=="PASS_CURRENT_SELECTORS_STOP_AT_SPIN6_NOT_SM"
    assert r["exact_chain"][-1]["dimension"]==15
    assert r["standard_model_target"]["dimension"]==12
