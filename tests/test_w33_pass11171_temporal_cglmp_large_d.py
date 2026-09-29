"""Regression for Pass 11171: temporal CGLMP to d = 24 (frozen scan), below 4 and monotone; optimiser gradient check."""
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11171_scan_temporal_cglmp_large_d as S  # noqa: E402
import w33_pass11171_temporal_cglmp_large_d as P  # noqa: E402


def test_summary():
    r = P.summarize()
    assert r['monotone'] and r['below_4']
    assert abs(r['merged']['24'] - 3.8528) < 1e-3 and abs(r['merged']['6'] - 3.4552) < 1e-3
    assert all(abs(v) < 5e-3 for v in r['validation_vs_committed'].values())
    assert 0.9 < r['fit_d_ge_6']['to4_alpha'] < 1.05


def test_gradient_and_d3():
    rng = np.random.default_rng(1)
    d = 5
    C = S.coeff(d)
    U = [np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))[0] for _ in range(4)]
    f, Om = S.value_grad(U, C)
    gn = sum(np.linalg.norm(o) ** 2 for o in Om)
    e = 1e-6
    fp, _ = S.value_grad([U[i] @ expm(e * Om[i]) for i in range(4)], C)
    fm, _ = S.value_grad([U[i] @ expm(-e * Om[i]) for i in range(4)], C)
    assert abs((fp - fm) / (2 * e) - 2 * gn) < 1e-5
    assert abs(S.run((3, 0, 3000))[2] - 3.162806507765822) < 1e-8
