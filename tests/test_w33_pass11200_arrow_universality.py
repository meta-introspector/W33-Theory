"""Regression for Pass 11200: universality of A = n - c (qubits, 5-, 7-dits; five qutrits on all 940 classes) and the
large-n law.  Frozen data plus a live check of the law on all 11 classes of Sp(4,2) and 34 of PSp(4,5)."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11200_arrow_universality as U  # noqa: E402

D = json.loads((ROOT / "data" / "w33_pass11200_arrow_universality.json").read_text())


def test_frozen():
    seen = {(c['q'], c['n']) for c in D['classes']}
    assert {(2, 2), (2, 3), (2, 4), (2, 5), (5, 2), (7, 2)} <= seen
    for c in D['classes']:
        assert c['law_holds'] and c['total_ok'] and c['burnside_mean_invariant_planes'] == '1' and c['A_max'] == c['n']
    n5 = D['qutrit_n5_classes']
    assert n5['classes'] == 940 and n5['total_ok'] and n5['exact_everywhere']
    assert D['qutrit_burnside'] == {'n2': '1', 'n3': '1', 'n4': '1'}
    assert D['law_holds_everywhere']
    for v in D['qutrit_large_n'].values():
        assert 0.6 < v['E_c'] < 0.85 and 0.38 < v['P_A_equals_n'] < 0.53


def test_live_small_groups():
    for q, n in ((2, 2), (5, 2)):
        r = U.class_check(q, n, U.GAP_SMALL)
        assert r['law_holds'] and r['total_ok'] and r['burnside_mean_invariant_planes'] == '1'
