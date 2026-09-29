"""Regression for Pass 11143: space and time are two GQ(3,3) on the same 40 points sharing the 16 product lines."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11143_space_time_two_quadrangles as P  # noqa: E402


def test_two_quadrangles():
    r = P.summarize()
    assert r['spatial_lines'] == r['temporal_lines'] == 40 and r['common_lines'] == 16 and r['common_are_products']
    assert list(r['spatial_srg']) == list(r['temporal_srg']) == [12, [2], [4]]
    assert r['spatial_only_are_detminus1_graphs'] and r['temporal_only_are_SL23_graphs']
    assert r['clifford_unitaries_mod_phase'] == 216 and r['clifford_symplectic_classes'] == 24
    assert r['pdm_all_supported_on_temporal_lines'] and len(r['pdm_negativity']) == 1
    assert abs(r['pdm_negativity'][0] - 1) < 1e-6
    assert r['partial_transpose_lines_are_temporal'] and r['pt_states_equal_to_a_clifford_pdm'] == 216
    assert r['bijection_states_to_cliffords']
    assert r['gsp_spatial_order'] == 103680 and r['common_similitudes'] == 2304
