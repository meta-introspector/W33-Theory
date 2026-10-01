"""Regression for Pass 11246: no discrete escape from the light-colored-state obstruction."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_frozen():
    s = json.loads((ROOT / "data" / "w33_pass11246_discrete_anomaly_escapes.json").read_text())["summary"]
    assert s["protected"] == 16326 and s["escapes"] == 0 and s["z6i_protected"] == 0
