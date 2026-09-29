"""Regression for Pass 11120: rho-safe models have the sparsest Wilson-line cell structure (necessary, not sufficient)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11120_rho_cell_criterion as P  # noqa: E402


def test_cells():
    r = P.summarize()
    assert r['models'] == 491 and len(r['safe_models']) == 5
    assert r['safe_all_sparse'] and r['safe_norm1_cells_mixed_single_kappa']
    assert r['sparse_models'] == 20 and r['sparse_tachyonic'] == 15
    assert r['patterns']['(1.0, 1.667, 1.0)|safe=True'] == 5
