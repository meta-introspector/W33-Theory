"""Replay the exact global intersection, not a numerical search for minima."""
import importlib.util,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location('global_poly',ROOT/'analysis/w33_20261001_g26_global_polynomial_vacuum.py')
M=importlib.util.module_from_spec(S);S.loader.exec_module(M)

def test_global_selector_complete_intersection_and_certificate():
    got=M.payload()
    assert got==json.loads((ROOT/'data/w33_20261001_g26_global_polynomial_vacuum.json').read_text())
    e=got['exact_complete_intersection']
    assert e['projective_solutions']==72
    assert e['all_reduced'] and e['no_points_at_infinity'] and e['no_coordinate_zero']
    assert got['potential']['unit_T_angular_gradient_Gram']==[[16,0],[0,100]]
