import json
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"analysis"/"w33_pass11548_history_orientation_square_correction.py"
DATA=ROOT/"data"/"PART_W33_PASS11548_HISTORY_ORIENTATION_SQUARE_CORRECTION.json"


def test_pass11548_history_orientation_square_correction():
    subprocess.run([sys.executable,str(SCRIPT)],cwd=ROOT,check=True,timeout=60)
    c=json.loads(DATA.read_text())
    assert c["status"]=="PASS_OBJECTWISE_CORRECTION_OF_PASS11540_CHARACTER_ASSIGNMENT"
    assert c["supersedes"]["verdict"]=="WITHDRAWN"
    assert c["corrected_subgroups"]["repo_PSp_linear"]["definition"].startswith("ker(sigma)")
    assert c["corrected_subgroups"]["repo_PSp_linear"]["order"]==24
    assert c["corrected_subgroups"]["repo_PSp_linear"]["not_equal_SO"] is True
    assert c["corrected_subgroups"]["SO_3_3"]["order"]==24
    assert c["corrected_subgroups"]["common_A4_Omega"]["order"]==12
    assert c["corrected_orientation_dictionary"]["global_history_outer_bit"]["character"].startswith("sigma")
    assert c["corrected_orientation_dictionary"]["inner_line_orientation_bit"]["character"].startswith("pi")
    assert c["oriented_null_tetrahedra"]["PSp_orbits"]==[4,4]
    assert c["oriented_null_tetrahedra"]["PGSp_O_orbit"]==8
