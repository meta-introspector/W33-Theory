"""Regression for Pass 11098: tree-level Yukawa textures of the SO(16)xSO(16) A8 models -- no single heavy top."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11098_yukawa_textures as P  # noqa: E402

C = json.loads(P.OUT.read_text())


def test_epsilon_matrix_has_two_equal_singular_values():
    rng = np.random.default_rng(7)
    for _ in range(10):
        sv = sorted(P.epsilon_singular_values(rng.standard_normal(3) + 1j * rng.standard_normal(3)))
        assert np.allclose(sv, [0, 1, 1])


def test_no_single_heavy_top():
    assert C["models"] == 104 and C["up_rank_exactly_one"] == 0
    assert C["up_top_charm_degenerate_untwisted"] == 73 and C["up_twisted_models"] == 31
    assert C["o1_rank_distribution"]["down"] == {"0": 77, "3": 27}
    assert C["o1_rank_distribution"]["lepton"] == {"0": 53, "1": 24, "3": 27}


def test_sample_models_recomputed():
    for v in C["sample_recheck"].values():
        assert v["recomputed"] == v["frozen"]
