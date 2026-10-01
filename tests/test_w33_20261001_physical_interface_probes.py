"""Independent exact field moments and curved operator-response controls."""
from pathlib import Path
import importlib.util,json
ROOT=Path(__file__).resolve().parents[1]
def mod(name):
    p=ROOT/'analysis'/('w33_20261001_'+name+'.py')
    spec=importlib.util.spec_from_file_location(name,p)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m

def test_actual_signed_field_gram_and_balancing_covariance():
    m=mod('cartan_dflat_resource_bridge');got=json.loads(json.dumps(m.payload(),default=str))
    stored=json.loads((ROOT/'data/w33_20261001_cartan_dflat_resource_bridge.json').read_text())
    assert got==stored
    assert got['exact']['normalized_color_purity']=='92/225'
    assert got['SL3_balancing']['cubic_operator_rank_before_after']==[78,78]
    assert got['SL3_balancing']['cubic_operator_congruence_error']<1e-12

def test_curved_Dirac_heat_response_and_independent_matrix_control():
    m=mod('isovolume_dirac_gravity_response');got=json.loads(json.dumps(m.payload()))
    stored=json.loads((ROOT/'data/w33_20261001_isovolume_dirac_gravity_response.json').read_text())
    assert m.compare_certificate(got,stored)
    assert got['analytic']['Gaussian_integral_coefficients'][3]=='1/140'
    assert got['independent_finite_epsilon_control']['absolute_error']<1e-5
    row=got['numerical_curvature_replay'][-1]
    assert abs(row['normalized_curvature_coefficient']-1)<5e-6
