"""Regression for Pass 11153: in time, Pauli measurements with any start give no violation."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11153_temporal_magic_budget as P  # noqa: E402


def test_temporal_magic():
    r = P.summarize()
    assert abs(r['temporal_pauli_measurements_any_start'] - 2) < 1e-9 and r['in_time_magic_must_be_in_measurements']
