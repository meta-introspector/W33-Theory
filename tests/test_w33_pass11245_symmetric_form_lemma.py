"""Regression for Pass 11245: alternating / symmetric nondegenerate forms over F_2 = 2^-k (k even)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_lemma():
    import w33_pass11245_symmetric_form_lemma as L
    r = L.run()
    assert r["all_k_le_12"] and all(b["matches_formula"] for b in r["brute"])
