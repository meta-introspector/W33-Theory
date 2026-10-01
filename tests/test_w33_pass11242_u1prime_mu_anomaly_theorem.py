"""Regression for Pass 11242: a mu-protecting anomaly-free U(1)' always leaves a colored state massless."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_one_model_and_anomaly():
    import w33_pass11242_u1prime_mu_anomaly_theorem as T
    import w33_pass10960_matter_even_dflat_closure as P
    L, _ = P.load_ledger()
    n = "Z6-II|Z6II_06__SM_20260917_1204"
    r = T.check_model(n, L[n])
    assert r["protected"] == 6 and r["theorem_violations"] == 0
    A, _ = T.su3_anomaly_vector(L[n])
    assert all(a == 0 for a in A[1:])


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11242_u1prime_mu_anomaly_theorem.json").read_text())
    s = d["summary"]
    assert s["su3_anomaly_zero_on_all_nonanomalous_U1s"] and s["theorem_violations"] == 0
    assert s["protected_cases"] == 812
