"""Regression for Pass 11247: no R-charges recorded for Z6-I; all Z6-II carry them."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_data_gap():
    import w33_pass11247_z6i_r_symmetry_data as Z
    r = Z.run()
    assert r["z6i_models"] == 87 and r["z6i_with_R"] == 0
    assert r["z6ii_with_R"] == r["z6ii_models"] == 128
