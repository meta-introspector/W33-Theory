"""Regression for Pass 11174: the forced pattern is Sym^2 of the det -1 swap; 36 four-clock potentials are perfect."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11174_lightcone_reflection_time_reversal as P  # noqa: E402


def test_reflection():
    r = P.summarize()
    assert r['lorentz_images'] == 24 and r['lightcone_reflection_is_sym2_swap'] and r['vkv_pi'] == [2, 1, 0]
    assert len(r['coordinate_permutations']) == 4
    assert {p['det_g'] for p in r['coordinate_permutations'] if p['L'] != [[1, 0, 0], [0, 1, 0], [0, 0, 1]]} == {2}
    assert r['four_clock_perfect'] == 36 and r['four_clock_rule_exact']
    assert r['single_oblique_clock_perfect'] == [True, True] and r['single_lightcone_clock_perfect'] == [False, False]
