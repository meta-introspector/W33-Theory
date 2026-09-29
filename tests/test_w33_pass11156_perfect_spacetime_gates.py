"""Regression for Pass 11156: det S_AA + det S_BA = 1; 13824 perfect gates, one local orbit; tetracode gate."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11156_perfect_spacetime_gates as P  # noqa: E402


def test_perfect():
    r = P.summarize()
    assert r['sp43_order'] == 51840 and r['formula_checked_on_random_cliffords']
    assert r['det_law_detAA_plus_detBA_eq_1'] and r['antiunitary_det_law_eq_minus1']
    assert r['perfect_count'] == 13824 == 24 ** 3 and r['perfect_all_blocks_invertible'] == 13824
    assert r['perfect_block_dets'] == {'(2, 2)': 13824} and r['perfect_local_double_cosets'] == [13824]
    assert r['tetracode_gate_is_antiunitary'] and r['tetracode_gate_is_perfect']
    assert r['tetracode_times_local_time_reversal_is_symplectic'] and r['tetracode_times_local_time_reversal_is_perfect']
    assert r['rank_census'] == {'space0_diag2': 576, 'space1_diag2': 18432, 'space2_diag0': 576, 'space2_diag1': 18432, 'space2_diag2': 13824}
