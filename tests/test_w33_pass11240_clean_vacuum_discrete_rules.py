"""Regression for Pass 11240: the clean U(1)' vacuum of Pass 11232 fails once integer and discrete rules are used."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_recomputed():
    import w33_pass11240_clean_vacuum_discrete_rules as C
    r = C.run()
    t = r["by_ruleset"]
    assert r["choices"] == 4
    assert all(v["protected"] == 4 and v["clean"] == 0 for v in t.values())
    assert t["U(1) lattice only"]["light"] == {"v": 4}
    assert set(t["+ space group + R"]["light"]) == {"q", "d", "v", "extra Higgs"}


def test_frozen_matches():
    d = json.loads((ROOT / "data" / "w33_pass11240_clean_vacuum_discrete_rules.json").read_text())
    assert d["model"] == "Z6-II|Z6II_06__SM_20260917_1204"
    assert all(v["clean"] == 0 for v in d["by_ruleset"].values())
