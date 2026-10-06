import json
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"analysis"/"w33_pass11540_history_vo33_orientation_square.py"
DATA=ROOT/"data"/"PART_W33_PASS11540_HISTORY_VO33_ORIENTATION_SQUARE.json"


def test_pass11540_history_vo33_orientation_square_corrected():
    subprocess.run([sys.executable,str(SCRIPT)],cwd=ROOT,check=True,timeout=60)
    c=json.loads(DATA.read_text())
    assert c["status"]=="PASS_EXACT_VO33_WITH_CHARACTER_ASSIGNMENT_CORRECTED_BY_11548"
    assert c["correction_notice"]["corrected_by_pass"]==11548
    assert c["history_graph"]["standard_name"]=="parabolic affine orthogonal polar graph VO(3,3)"
    assert c["history_graph"]["vertices"]==27
    assert c["history_graph"]["degree"]==8
    assert c["history_graph"]["edges"]==108
    assert c["orthogonal_group"]["O_3_3_order"]==48
    assert c["orthogonal_group"]["SO_3_3_order"]==24
    assert c["orthogonal_group"]["Omega_3_3_order"]==12
    assert c["affine_orthogonal_ladder"]["AffO_order"]==1296
    assert c["affine_orthogonal_ladder"]["AffSO_order"]==648
    assert c["affine_orthogonal_ladder"]["AffOmega_order"]==324
    assert c["repo_stabilizer_ladder"]["PSp_order"]==648
    assert "ker(sigma)" in c["repo_stabilizer_ladder"]["PSp_linear"]
    assert "not the repository PSp" in c["affine_orthogonal_ladder"]["warning"]
    assert c["orientation_square_corrected"]["determinant"].startswith("sigma*pi")
