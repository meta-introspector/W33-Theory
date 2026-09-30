"""Regression for Pass 11184: independent orbit recomputation and the reversal-pair structure (frozen data; the full
recomputation takes ~25 minutes)."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11184_merged_orbitals.json").read_text())
    assert d['orbits'] == 20 and d['matches_gap_subdegrees'] and d['fine_invariant_constant_on_orbits']
    assert d['self_paired_orbitals'] == 14 and d['non_self_paired_sizes'] == [256, 256, 2304, 2304, 6912, 6912]
    assert d['remaining_ties'] == 2 and d['remaining_ties_are_reversal_pairs']
