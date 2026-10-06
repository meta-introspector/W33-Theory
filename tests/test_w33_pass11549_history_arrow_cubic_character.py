import json
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"analysis"/"w33_pass11549_history_arrow_cubic_character.py"
DATA=ROOT/"data"/"PART_W33_PASS11549_HISTORY_ARROW_CUBIC_CHARACTER.json"


def test_pass11549_history_arrow_cubic_character():
    subprocess.run([sys.executable,str(SCRIPT)],cwd=ROOT,check=True,timeout=60)
    c=json.loads(DATA.read_text())
    assert c["status"]=="PASS_CLOSED_FORM_HISTORY_ARROW_CUBIC_CHARACTER"
    assert c["hamming_arrow"]["exact_old_edge_matches"]==108
    assert c["hamming_arrow"]["antisymmetry"]=="chi(-d)=-chi(d)"
    assert c["triangle_foliation"]["temporal_triangles"]==36
    assert c["triangle_foliation"]["each_edge_in_exactly_one_triangle"] is True
    assert c["triangle_foliation"]["factorization"]=="36 = 4 * 9"
    assert c["triangle_foliation"]["chi_oriented_triangle_boundaries_equal_chi_edge_cycle"] is True
    assert c["symmetry"]["repo_PSp_affine_order"]==648
    assert c["symmetry"]["PSp_action_on_chi"].startswith("preserves")
    assert c["symmetry"]["outer_action_on_chi"]=="negates"
