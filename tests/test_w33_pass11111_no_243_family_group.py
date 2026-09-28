"""Regression for Pass 11111: three generations cannot carry 3^(1+4), alone or glued to colour."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11111_no_243_family_group as P  # noqa: E402


def test_extraspecial_degrees():
    a = P.part_a()
    assert a['k=1']['order'] == 27 and a['k=2']['order'] == 243
    assert a['k=1']['defining_irreducible'] and a['k=2']['defining_irreducible']
    assert a['three_generations_allow_k'] == [1]


def test_no_gauge_centre_gluing():
    assert P.gluing_solutions() == []
    assert len(P.gluing_solutions({k: (0,) + v[1:] for k, v in P.QUARKS.items()})) == 6    # control: the SM Z6 centre
    import json
    q = json.loads(P.OUT.read_text())['B_quark_sectors']
    assert len(q) == 12 and all(v['Q'] == [2] and v['ubar'] == [2] and 2 in v['dbar'] for v in q.values())
