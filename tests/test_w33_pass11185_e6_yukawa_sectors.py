"""Regression for Pass 11185: E6 Yukawa sectors of two-qutrit splits and the gate-class transitions."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11185_e6_yukawa_sectors as P  # noqa: E402


def test_sectors():
    r = P.summarize()
    assert r['factorisations_are_frame_triangles'] and r['decomposition'] == [1, 10, 16]
    assert r['sectors'] == {'S1': 5, 'S16': 40} and r['sector_rule_holds']
    tr = r['transitions_per_split']
    assert tr['S1->S1|collinear'] == 4 and tr['S1->S16|collinear'] == 8 and tr['S1->S16|noncollinear'] == 32
    assert 'S1->S1|noncollinear' not in tr
    assert tr['S16->S1|collinear'] == 1 and tr['S16->S1|noncollinear'] == 4
    g = r['gate_level_from_S1']
    assert g == {'noncollinear->S16': 36864, 'collinear->S1': 4608, 'collinear->S16': 9216, 'equal->S1': 1152}
    assert r['frame_stabiliser_order'] == 1920
    assert not any('noncollinear' in k for k in r['frame_stabiliser_profile'])
