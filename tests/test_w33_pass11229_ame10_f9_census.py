"""Regression for Pass 11229: independent AME(10,3) census -- every state found is F9-linear (Glynn)."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11229_ame10_f9_census.json").read_text())
    assert d["ame_found"] == 138 and d["new_vs_master_census"] == 138
    assert d["classes"] == {"structures=4,schur=10": 138} and d["all_f9_linear_glynn"]
