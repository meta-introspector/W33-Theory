import json
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

CASES=[
 ("11551","w33_pass11551_hamming_arrow_albert_diagonal_cubic.py",
  "PART_W33_PASS11551_HAMMING_ARROW_ALBERT_DIAGONAL_CUBIC.json",
  "PASS_HAMMING_ARROW_IS_ALBERT_DIAGONAL_CUBIC"),
 ("11552","w33_pass11552_hamming_refinement_continuum_firewall.py",
  "PART_W33_PASS11552_HAMMING_REFINEMENT_CONTINUUM_FIREWALL.json",
  "PASS_NAIVE_FULL_NULL_TOWER_IS_EXPANDER_LIKE"),
 ("11553","w33_pass11553_history_arrow_equivariant_cohomology.py",
  "PART_W33_PASS11553_HISTORY_ARROW_EQUIVARIANT_COHOMOLOGY.json",
  "PASS_UNIQUE_PSP_FIXED_MOD3_H1_ARROW_CLASS"),
 ("11554","w33_pass11554_tetracode_orientation_double_cover.py",
  "PART_W33_PASS11554_TETRACODE_ORIENTATION_DOUBLE_COVER.json",
  "PASS_NON_SPLIT_GL23_TO_S4_ORIENTATION_EXTENSION"),
 ("11555","w33_pass11555_temporal_local_dynamics_dirac.py",
  "PART_W33_PASS11555_TEMPORAL_LOCAL_DYNAMICS_DIRAC.json",
  "PASS_EXACT_TEMPORAL_LOCAL_OPERATOR_ALGEBRA"),
]

def test_pass11551_11555_five_attack_packet():
    cert={}
    for key,script,data,status in CASES:
        subprocess.run([sys.executable,str(ROOT/"analysis"/script)],cwd=ROOT,check=True,timeout=90)
        cert[key]=json.loads((ROOT/"data"/data).read_text())
        assert cert[key]["status"]==status

    assert cert["11551"]["clock_albert_frame"]["all_ternary_vectors_checked"]==27
    assert cert["11551"]["clock_albert_frame"]["integer_before_mod3"] is True

    assert cert["11552"]["exact_full_null_relation"]["consequence"].endswith("normalized absolute spectral gap -> 1")
    assert cert["11552"]["nearest_null_distance3"]["weight1_normalized_gap"]=="9/(2n)"

    assert cert["11553"]["complex"]["H1_dimension_F3"]==46
    assert cert["11553"]["cochain_side"]["PSp_fixed_H1_dimension"]==1
    assert cert["11553"]["cochain_side"]["arrow_spans_fixed_line"] is True
    assert cert["11553"]["chain_side"]["homology_warning"].startswith("the arrow is homologically trivial")

    assert cert["11554"]["extension"]["splits"] is False
    assert cert["11554"]["extension"]["F2_augmented_rank"]==cert["11554"]["extension"]["F2_coboundary_rank"]+1

    assert cert["11555"]["operator_D"]["rank"]==8
    assert cert["11555"]["operator_D"]["exact_square_identity"]=="3 D^2 = L(L-6I)(L-12I)"
    assert cert["11555"]["local_kernel_classification"]["translation_and_PSp_invariant_nearest_null_kernel_dimension"]==2


def test_pass11551_11555_correlated_strengthening():
    script=ROOT/"analysis"/"w33_pass11551_11555_hamming_arrow_five_attacks.py"
    data=ROOT/"data"/"PART_W33_PASS11551_11555_HAMMING_ARROW_FIVE_ATTACKS.json"
    subprocess.run([sys.executable,str(script)],cwd=ROOT,check=True,timeout=90)
    c=json.loads(data.read_text())
    assert c["status"]=="PASS_FIVE_HAMMING_ARROW_ATTACKS"
    assert c["passes"]["11553"]["mod3"]["PSp_fixed_H1_dimension"]==1
    assert c["passes"]["11553"]["prime_control"]["2"]["PSp_fixed_H1_dimension"]==0
    assert c["passes"]["11553"]["prime_control"]["5"]["PSp_fixed_H1_dimension"]==0
    assert c["passes"]["11553"]["prime_control"]["7"]["PSp_fixed_H1_dimension"]==0
    assert c["passes"]["11554"]["tetracode_cover"]["non_split"] is True
    assert c["passes"]["11554"]["history_cover"]["split"].startswith("O(3,3)=")
    assert c["passes"]["11554"]["fiber_product"]["is_archived_tomotope"] is False
    assert c["passes"]["11555"]["operator"]["exact_polynomial_identity"]=="3 D^2 = L(L-6I)(L-12I)"
