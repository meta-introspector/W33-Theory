import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "analysis" / "w33_pass11540_history_vo33_orientation_square.py"
DATA = ROOT / "data" / "PART_W33_PASS11540_HISTORY_VO33_ORIENTATION_SQUARE.json"


def test_pass11540_history_vo33_orientation_square():
    subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, timeout=60)
    cert = json.loads(DATA.read_text())

    assert cert["status"] == "PASS_EXACT_VO33_ORTHOGONAL_ORIENTATION_SQUARE"
    assert cert["history_graph"]["standard_name"] == "parabolic affine orthogonal polar graph VO(3,3)"
    assert cert["history_graph"]["vertices"] == 27
    assert cert["history_graph"]["degree"] == 8
    assert cert["history_graph"]["edges"] == 108

    assert cert["orthogonal_group"]["O_3_3_order"] == 48
    assert cert["orthogonal_group"]["SO_3_3_order"] == 24
    assert cert["orthogonal_group"]["Omega_3_3_order"] == 12
    assert cert["orthogonal_group"]["Omega_equals_derived_SO"] is True

    assert cert["affine_group_ladder"]["AffO_order"] == 1296
    assert cert["affine_group_ladder"]["AffSO_order"] == 648
    assert cert["affine_group_ladder"]["AffOmega_order"] == 324

    assert cert["orientation_square"]["character_pair_census"] == {
        "++": 12,
        "+-": 12,
        "-+": 12,
        "--": 12,
    }
    assert cert["orientation_square"]["combined_quotient"] == "AffO/AffOmega ~= C2 x C2"
