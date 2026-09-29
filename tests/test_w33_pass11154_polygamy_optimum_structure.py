"""Regression for Pass 11154: invariants of the N_AB + N_AC optimum are fifteenths (N_AB = N_AC = 2/sqrt15)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11154_polygamy_optimum_structure as P  # noqa: E402


def test_structure():
    assert P.summarize()['all_closed_forms_within_1e_7']
