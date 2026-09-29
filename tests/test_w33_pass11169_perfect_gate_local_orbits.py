"""Regression for Pass 11169: six local orbits of perfect three-qutrit gates, labelled by pi in S3."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11169_perfect_gate_local_orbits as P  # noqa: E402


def test_orbits():
    r = P.summarize()
    assert r['n_classes'] == 6 and set(r['stabiliser_orders'].values()) == {4}
    assert r['each_class_one_orbit'] and r['class_size'] == 47775744
    assert r['two_qutrit_stabiliser'] == 24 and r['two_qutrit_orbit'] == 13824
    assert set(r['f9_pattern_counts'].values()) == {2048} and r['f9_total'] == 12288
