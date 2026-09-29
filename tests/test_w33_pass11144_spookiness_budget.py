"""Regression for Pass 11144: W(3,3) resources are classical in space and time; temporal qutrit CGLMP exceeds spatial."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11144_spookiness_budget as P  # noqa: E402


def test_budget():
    r = P.summarize()
    assert abs(r['lhv_max'] - 2) < 1e-12
    assert abs(r['spatial_omega_cglmp_bases'] - 2.8729340512) < 1e-8
    assert abs(r['spatial_optimum_adgl'] - 2.9148542155) < 1e-8
    assert abs(r['temporal_mixed_equals_spatial'] - r['spatial_omega_cglmp_bases']) < 1e-12
    assert abs(r['spatial_stabilizer_max'] - 2) < 1e-9 and abs(r['temporal_stabilizer_max'] - 2) < 1e-9
    assert abs(r['pauli_point_value'] - 2) < 1e-9 and abs(r['magic_path_curvature'] - r['pi2_over_6']) < 1e-5
    assert 3.162 < r['temporal_optimum'] < 3.1635 and r['temporal_optimum'] > r['spatial_optimum_adgl']
    assert r['classical_invasive_max_1state'] == 2.0 and r['classical_invasive_max_2states'] == 4.0
