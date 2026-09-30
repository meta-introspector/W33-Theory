"""Regression for Pass 11192: time reversals (anti-symplectic maps) and the splits they respect."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11192_time_reversal_mereology as T  # noqa: E402

D = json.loads((ROOT / "data" / "w33_pass11192_time_reversal_mereology.json").read_text())


def test_frozen():
    assert D['one_qutrit'] == {'1': 12, '0': 12}
    tab = D['n2_exhaustive']['table']
    assert tab['sq-1|local15|fixing0'] == 72 and tab['sq1|local7|fixing6'] == 1080 and tab['sq0|local0|fixing0'] == 21888
    n2, n3 = D['n2'], D['n3']
    assert n2['outer_classes'] == 10 and n2['total_ok'] and n2['kramers_classes'] == 1
    assert n2['local_nowhere_fraction'] == '19/45'
    assert n3['outer_classes'] == 26 and n3['total_ok'] and n3['kramers_classes'] == 0
    assert n3['local_nowhere_fraction'] == '8377/12285'
    for d in (n2, n3):
        assert d['all_anti_symplectic'] and d['every_involution_is_local_conjugation'] and d['kramers_only_pairing']


def test_live_two_qutrit_classes():
    Bs = np.array(T.M.factorisations(2))
    rows = T.load_outer(2)
    assert sum(r['size'] for r in rows) == 25920
    sizes = {}
    for r in rows:
        assert T.anti_symplectic(r['T'], 2)
        if r['order'] == 2:
            sizes[r['size']] = (T.square_sign(r['T'], 2), T.locality(r['T'], Bs, 2)[:2])
    assert sizes == {36: (-1, (15, 0)), 540: (1, (7, 6))}
