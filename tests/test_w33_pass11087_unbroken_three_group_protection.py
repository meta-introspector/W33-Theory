"""Regression for Pass 11087: the protecting unbroken 3-group, block by block."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
C = json.loads((ROOT / "data" / "w33_pass11087_unbroken_three_group_protection.json").read_text())


def test_group_orders():
    s = C["summary"]
    assert s["condensate"]["vacua"] == 55 and s["singlet"]["vacua"] == 28
    assert s["condensate"]["orders"] == {"729": 48, "81": 6, "27": 1}
    assert s["singlet"]["orders"] == {"81": 28} and s["singlet"]["unbroken_u1s"] == {"2": 28}
    for v in C["vacua"]:
        d = v["order_with_block_dropped"]
        assert d["PG"]["order"] == v["order"] and d["centre"]["order"] == v["order"]   # no independent part
        if v["tag"] == "condensate" and "c1__" not in v["model"]:                      # the 54 SU(4)-condensed
            assert d["SG"]["order"] < v["order"] and d["R"]["order"] < v["order"]


def test_protection_is_interlocking():
    p = C["summary"]["condensate"]["charge_third_pairs"]
    assert p["forbidden"] == 81876
    necessary = sum(v for k, v in p.items() if k.startswith("necessary"))
    assert p["redundant (no single block necessary)"] + necessary == p["forbidden"]
    sg = sum(v for k, v in p.items() if k.startswith("necessary") and "SG" in k)
    u1 = sum(v for k, v in p.items() if k.startswith("necessary") and "U1" in k)
    assert sg / p["forbidden"] < 0.2 < u1 / p["forbidden"]


def test_singlet_vacua_have_a_z_prime_coupled_to_matter():
    assert C["singlet_vacua_extra_u1"] == {"extra U(1)' coupling to SM matter": 28}
