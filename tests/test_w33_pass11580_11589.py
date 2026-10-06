import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
P=json.loads((ROOT/'data/PART_W33_PASS11580_11589_SPIN10_SM_BREAKTHROUGH.json').read_text())['passes']

def test_11580_spin10_patisalam_hypercharge():
    a=P['11580']
    assert a['centralizer_dimension']==7
    assert a['centralizer_centroid_split']==[3,3,1]
    assert a['sixY_spectrum']=={'1':6,'-3':2,'-4':3,'2':3,'6':1,'0':1}

def test_11581_global_z6_and_anomalies():
    a=P['11581']
    assert a['kernel_order']==6
    assert all(v==0 for v in a['anomalies'].values())
    assert a['weak_doublet_count']==4

def test_11582_overlap_refinement():
    rows={(r['L'],r['q']):r for r in P['11582']['scans']}
    assert rows[(4,1)]['index']==-1
    assert rows[(4,2)]['index']==-4
    assert rows[(5,3)]['index']==-9
    assert rows[(5,4)]['index']==-16
    assert P['11582']['max_GW_residual_observed']<1e-12

def test_11583_clock_chirality_cycle():
    a=P['11583']
    assert a['spin10_chirality_eigenspaces']==[16,16]
    assert a['all_45_spin10_bivectors_commute_with_chirality']
    assert a['one_tick_anticommutator_error']<1e-12
    assert a['induced_complex_structure_square_minus_identity_error']<1e-12
    assert a['two_tick_flip_error']<1e-12

def test_11584_unique_10H_channel():
    a=P['11584']
    assert a['symmetric_square_multiplicities']==[10,126]
    assert a['exterior_square_multiplicity']==120

def test_11585_coarse_gravity_firewall():
    a=P['11585']
    assert a['exact_weighted_intR']=='0'
    assert a['exact_unweighted_sumR']=='-360725/209088'
    assert a['heat_delta_quadratic_coefficients']['0.1']>8

def test_11586_11588_charge_family_product():
    assert P['11586']['qBL_gcd']==1 and P['11586']['sixY_gcd']==1
    assert P['11586']['hypercharge_quantum']=='1/6'
    assert P['11587']['one_E6_27_dimension']==27
    assert P['11587']['three_independent_E6_families_dimension']==81
    assert P['11588']['external_overlap_index']==-1
    assert P['11588']['internal_spin10_weyl_rank']==16
    assert P['11588']['tensor_index_multiplicity']==-16
