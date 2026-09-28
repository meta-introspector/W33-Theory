"""Regression for Pass 11108: tachyon-freedom over the Wilson-line moduli, 104 models (+ rescan when frozen)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11108_tachyon_free_moduli_rescan as P  # noqa: E402


def test_counts_104():
    r = P.summarize()
    assert r['A_tachyon_free_counts_104'] == {'rho': 1, '1.0': 55, '1.5': 55, '2.0': 104}
    assert r['A_free_at_rho'] == ['18'] and r['A_survivors_free_at_rho'] == []
    assert r['A_levels_at_rho'] == {'-0.0556': 76, '-0.1667': 27, 'None': 1}


def test_enumerator_on_model_18_rows():
    import json
    rows = json.loads((ROOT / 'data' / 'w33_pass11108_model_shifts_104.json').read_text())['models']['18']['rows']
    assert P.levels(rows)['rho'] is None


def test_rescan_candidates_all_untwisted():
    r = P.summarize()
    assert r['B_new_models'] >= 123 and r['B_tachyon_free_at_rho'] >= 3
    assert r['rho_safe_models_all_untwisted_quarks']
    assert r['A_twisted_up_models'] == 31 and r['A_twisted_up_free_at_rho'] == []
    for v in r['B_candidates'].values():
        assert v['down_tree_level'] == 0 and all(e[1] == 0 for e in v['up_exponents'])     # top = charm
