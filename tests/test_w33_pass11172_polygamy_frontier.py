"""Regression for Pass 11172: the polygamy frontier -- corner (1,0) up to lambda_c = 0.92143, family = global."""
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11172_polygamy_frontier as P  # noqa: E402


def test_frontier():
    assert abs(P.family_max(0.5)[0] - 1) < 1e-12
    assert abs(P.family_max(1.0)[0] - 4 / math.sqrt(15)) < 1e-12
    lc = P.critical_lambda()
    assert abs(lc - 0.92143357) < 1e-6
    rng = np.random.default_rng(3)
    for lam in (0.5, 0.96, 1.0):
        g = max(P.seesaw_weighted(rng, lam) for _ in range(12))
        assert g <= P.family_max(lam)[0] + 1e-9
