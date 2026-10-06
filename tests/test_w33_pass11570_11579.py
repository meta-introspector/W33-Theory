import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
R=json.loads((ROOT/'data/PART_W33_PASS11570_11579_EXECUTABLE_GAUGE_GRAVITY_OVERLAP.json').read_text())
P=R['passes']

def test_11570_exact_sm_intersection():
    a=P['11570']
    assert a['F4_H3C_stabilizer']['dimension']==16
    assert a['F4_H3C_stabilizer']['centroid_eigenspace_dimensions']==[8,8]
    assert a['spin9_intersection']['dimension']==12
    assert a['spin9_intersection']['center_dimension']==1
    assert a['spin9_intersection']['derived_centroid_eigenspace_dimensions']==[8,3]

def test_11571_curvature_is_nontrivial():
    a=P['11571']
    assert a['curvature_nonzero_site_plaquette_pairs']==81
    assert a['distinct_scalar_curvature_values']==11
    assert a['scalar_curvature_sum']=='-360725/209088'

def test_11573_11578_gw_firewall():
    assert max(abs(x['index']) for x in P['11573']['scans'])<1e-12
    assert P['11578']['max_GW_residual_across_executed_scans']<1e-12

def test_11574_selector():
    a=P['11574']
    assert a['su3_ideal']['fixed_subspace_dimension']==3
    assert a['su2_ideal']['fixed_subspace_dimension']==6

def test_11575_11576_matter_firewalls():
    assert P['11575']['exact_SM_commutant_on_Peirce16_dimension']==4
    assert P['11576']['trace_Q']==0
    assert P['11576']['trace_Q_cubed']==0
    assert P['11576']['symmetric_cubic_gauge_anomaly_nonzero_entries']==0

def test_11577_11579_chain():
    assert P['11577']['centralizer_dimension']==9
    assert 'u(3)' in P['11577']['identification']
    assert len(P['11579']['chain'])==4
