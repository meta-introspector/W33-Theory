"""Rebuild the actual gauge map; keep analytic/interval and numerical scopes separate."""
from pathlib import Path
import importlib.util,json
ROOT=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location('dflat_mass_loop',ROOT/'analysis/w33_20261001_dflat_cartan_mass_loop.py')
M=importlib.util.module_from_spec(S);S.loader.exec_module(M)

def test_full_plane_moments_mass_pencil_and_interval_saddle():
    got=M.payload();stored=json.loads((ROOT/'data/w33_20261001_dflat_cartan_mass_loop.json').read_text())
    assert M.compare_certificate(got,stored)
    assert got['normalization']['all_Cartan_cross_moment_max']<1e-9
    assert got['mass_interface']['mirror_spectrum_error']<1e-8
    assert got['exact_mirror_model']['trace_mass_fourth']==8
    w=got['loop_test']['sign_witnesses']
    assert w[0]['Rayleigh_bounds'][1]<0 and w[1]['Rayleigh_bounds'][0]>0

    d=got['dimension_five_light_mass']
    assert d['rank_before_after']==[78,79]
    assert abs(d['new_light_singular_value']-.001)<1e-10
    assert d['old_heavy_singular_value_error']<1e-10
    assert d['squared_mass_identity_error']<1e-10

    block=got['mass_to_gravity_Dirac']
    assert block['finite_Hilbert_dimension']==162
    assert block['self_adjoint_residual']==0
    assert abs(block['unit_radius_control_moments'][0]-40.000002)<1e-10
    assert abs(block['unit_radius_control_moments'][1]-(16+2e-12))<1e-10
