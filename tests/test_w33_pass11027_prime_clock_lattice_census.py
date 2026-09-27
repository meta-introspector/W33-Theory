import importlib.util, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p11027",ROOT/"analysis"/"w33_pass11027_prime_clock_lattice_census.py")
P=importlib.util.module_from_spec(S); S.loader.exec_module(P)
C=json.loads((ROOT/"data"/"w33_pass11027_prime_clock_lattice_census.json").read_text())

def test_fast_payload_replays():
    assert P.payload()==C

def test_exact_minima_through_13():
    rows=C["exact_census"]
    assert [r["q"] for r in rows]==[3,5,7,11,13]
    assert [r["minimum_glue_norm"] for r in rows]==[2,4,6,10,12]
    assert [r["minimum_codewords"] for r in rows]==[8,72,240,552,756]

def test_q3_is_unique_root_creation_rung():
    rows=C["general_root_firewall"]["rows"]
    assert rows[0]["can_create_norm2_roots_by_bound"] is True
    assert not any(r["can_create_norm2_roots_by_bound"] for r in rows[1:])

def test_root_shells():
    rs=C["root_shells"]
    assert rs[0]["root_count"]==240
    assert all(r["root_count"]==r["q"]*r["rank"] for r in rs[1:])

def test_q13_deep_denominator_frozen():
    assert C["q13_reproduction"]["nonzero_codewords_checked"]==62_748_516
