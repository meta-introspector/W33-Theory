import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def test_exhaustive_restored_r_colored_obstruction():
    d=json.loads((ROOT/"data/w33_20261001_z6i_colored_escape_exhaustive.json").read_text())
    s=d["summary"]
    assert s["models"]==23
    assert s["vacua"]==6_695_116
    assert s["protected_Hu_cases"]==6_695_116
    assert s["colored_escapes"]==0
    assert s["models_with_escape"]==[]
    assert s["unique_support_lattices"]==1464
    assert all(d["checks"].values())
