"""Regression for Pass 11106: one-loop potential of the SO(16)xSO(16) A8 survivors at the orbifold point."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11106_one_loop_orbifold_point as P  # noqa: E402


def test_engine_consistency():
    c = P.checks()
    assert c["gamma1_6_max_dev"] < 1e-10                      # seed block Gamma_1(6)-invariant
    assert abs(c["S_phase"][0] - 1) < 1e-9 and abs(c["S_phase"][1]) < 1e-9
    assert abs(c["T_phase"][0] + 1) < 1e-9                    # Z[1,1] = -H[1,1]
    assert c["beta_S_invariance_dev"] < 1e-9


def test_potential_and_tachyons():
    r = P.summarize()
    assert r["family_modulus_min_at_rho"] and r["wilson_line_monotone"] and r["lambda_positive"]
    assert r["tachyon_free_at_2"] == 12
    assert len(r["tachyonic_at_some_radius_below_2"]) == 12
    m = r["massless_check_model2"]
    assert m["untwisted"] == [122, 122] and m["witten_twisted"] == [-68, -68] and m["theta_no_beta"] == [174, 174]
    assert m["corrected_model2"] == m["total_corrected"] == -252
    assert r["nB_minus_nF_corrected_range"] == [-480, -48] and r["nB_minus_nF_zero"] == 0
