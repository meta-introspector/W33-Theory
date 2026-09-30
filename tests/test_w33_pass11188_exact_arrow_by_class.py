"""Regression for Pass 11188: the exact intrinsic arrow on every conjugacy class (two qutrits recomputed, three
qutrits from the frozen certificate; the full three-qutrit scan needs the 110565 factorisations, ~5 minutes)."""
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_two_qutrit_exact_recomputed():
    import numpy as np
    import w33_pass11180_mereology as M
    import w33_pass11183_intrinsic_arrow as AR
    import w33_pass11188_exact_arrow_by_class as E
    cls = E.parse_classes()
    assert len(cls[2]) == 20 and len(cls[3]) == 74
    Bs = np.array(M.factorisations(2))
    J = M.form(2)
    dist = {}
    for c in cls[2]:
        S = c["S"]
        assert np.array_equal((S.T @ J @ S) % 3, J % 3)
        assert M.profile(S, Bs, 2)["local"] == c["fixed"]
        A = int(AR.exports(S, Bs, 2).sum(1).min())
        dist[A] = dist.get(A, 0) + c["size"]
    assert dist == {0: 9576, 2: 16344} and Fraction(16344, 25920) == Fraction(227, 360)


def test_three_qutrit_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11188_exact_arrow_by_class.json").read_text())["3"]
    dist = {int(k): v["elements"] for k, v in d["A_distribution"].items()}
    assert dist == {0: 286369344, 2: 2466749376, 3: 1832232960} and sum(dist.values()) == 4585351680
    for r in d["classes"]:
        if r["order"] % 7 == 0 or r["order"] % 13 == 0:
            assert r["A"] == 3
        if r["order"] % 5 == 0:
            assert r["A"] == 2
