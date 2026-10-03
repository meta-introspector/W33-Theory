"""Independent reconstructions of the stored operator witnesses."""
import json
import sys
from pathlib import Path
import numpy as np
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11390_11394_context_matter_clock import load

def packet():return json.loads((ROOT/'data/w33_pass11395_11402_joint_atlas_spin_portal.json').read_text())

def test_spin_double_cover_and_loop():
    p=packet()['spin'];qs=[list(map(sp.sympify,q))for q in p['binary_octahedral_quaternions']]
    assert len(qs)==48
    for q in qs:assert sp.simplify(sum(x*x for x in q))==1
    w,x,y,z=map(sp.sympify,p['loop']['quaternion'])
    assert w==0 and sp.simplify(x*x+y*y+z*z)==1
    U=sp.Matrix(p['loop']['su2']).applyfunc(sp.sympify)
    assert (U*U).applyfunc(sp.simplify)==-sp.eye(2)

def test_exact_one_portal_cubic():
    p=packet()['portal_polynomial'];g,x=sp.symbols('g x')
    assert sp.expand(sp.sympify(p['exact_characteristic_polynomial'])-x*(x-6)**2+3*g*g*(x-3)**2)==0
    for row in p['scans']:
        roots=np.array(row['squared_singular_values']);gg=row['g']
        assert np.max(abs(roots*(roots-6)**2-3*gg*gg*(roots-3)**2))<1e-5
        assert min(row['adjacent_singular_ratios'])<2
    y=sp.symbols('y');poly=y**4-6*y**3+108*y-324
    assert sp.Poly(poly,y).count_roots(3,6)==1

def test_selector_probabilities_from_full_graph():
    c=load();p=packet()['selector_portal'];A=np.array([[int(i!=j and bool(set(a)&set(b)))for j,b in enumerate(c['lines'])]for i,a in enumerate(c['lines'])])
    t=p['target_context']
    assert A[0,t]==0 and (A@A)[0,t]==4
    for row in p['scans']:
        _,v=np.linalg.eigh(np.diag([0]+[1]*39)-row['h_over_J']*A)
        assert abs(v[t,0]**2-row['target_probability'])<1e-14
    assert abs(p['scans'][-1]['probability_over_eta4']-16)<.3

def test_compiled_channels_reconstruct_operators():
    c=load();p=packet()['selector_portal']['three_channel_compiler'];blocks=[]
    for orb in p['orbit_pairs']:
        H=np.zeros((80,80))
        for a,b in orb:H[a,b]=H[b,a]=1
        blocks.append(c['plus'].conj().T@c['S'].conj().T@H@c['S']@c['minus'])
    channels=[sum(w*B for w,B in zip(weights,blocks))for weights in p['rank1_channel_weights']]
    for i,T in enumerate(channels):
        assert np.linalg.matrix_rank(T,tol=1e-10)==1
        assert abs(np.linalg.norm(T)-1)<1e-12
        for j,U in enumerate(channels):
            if i!=j:assert np.linalg.norm(T@U.conj().T)<1e-12 and np.linalg.norm(T.conj().T@U)<1e-12
    for row in p['scans']:
        H=sum(w*T for w,T in zip(row['probabilities'],channels))
        assert np.linalg.norm(np.linalg.svd(H,compute_uv=False)-row['singular_values'])<1e-12
    assert min(p['scans'][-1]['adjacent_ratios'])>10000

def test_dual_spread_reconstruction_and_dark_readout():
    c=load();p=packet()['chiral']['spread_extension'];Z=c['A'][:40,40:];B=np.zeros((40,36))
    for j,S in enumerate(p['spreads']):
        assert len(S)==10 and len(set().union(*(set(c['lines'][l])for l in S)))==40
        B[S,j]=1
    assert np.array_equal(3*Z.T@Z+B@B.T-3*np.ones((40,40)),18*np.eye(40))
    w,v=np.linalg.eigh(Z.T@Z);dark=v[:,w<1e-10]
    assert dark.shape[1]==15
    assert np.linalg.norm(dark.T@B@B.T@dark-18*np.eye(15))<1e-11

def test_extended_operator_has_nonabelian_kernel_multiplicity():
    c=load();p=packet()['chiral']['spread_extension'];B=np.zeros((40,36))
    for j,S in enumerate(p['spreads']):B[S,j]=1
    prior=json.loads((ROOT/'data/w33_pass11390_11394_context_matter_clock.json').read_text())['portal'];C=np.zeros((80,80))
    for a,b in prior['nonedge_orbit']:C[a,b]=C[b,a]=1
    H=np.zeros((116,116));H[:80,:80]=c['A']+.1*C;H[40:80,80:]=B;H[80:,40:80]=B.T
    si={tuple(S):i for i,S in enumerate(p['spreads'])};Ps=[]
    for g in c['x']['internal_h27']['group_permutations']:
        gp=list(g)+[80+si[tuple(sorted(g[40+l]-40 for l in S))]for S in p['spreads']]
        P=np.eye(116)[np.argsort(gp)];assert np.linalg.norm(P@H-H@P)<1e-12;Ps.append(P)
    w,v=np.linalg.eigh(H);K=v[:,abs(w)<1e-9];PN=sum(Ps)/27
    w,v=np.linalg.eigh(K.T@(np.eye(116)-PN)@K);K=K@v[:,w>.5]
    assert K.shape[1]==34 and abs(np.trace(K.T@Ps[1]@K)-7)<1e-10
    assert p['central_qutrit_multiplicities']==[3,3]
    assert p['linear_character_multiplicities']==[0]+[2]*8

def test_cycle_projector_and_control_overlap():
    p=packet();P=np.array(p['axis_bridge']['cycle_projector_160_numerator'],dtype=float)/160
    assert np.linalg.norm(P@P-P)<1e-12 and abs(np.trace(P)-81)<1e-12
    a,b=p['controls']['edge_control_pair'];u=P[:,a]/np.sqrt(P[a,a]);v=P[:,b]/np.sqrt(P[b,b])
    assert abs(u@v+1/3)<1e-12
    assert abs(2*(u@v)**2-1+7/9)<1e-12
    assert len(p['atlas']['cross_chart_histogram'])==2 and sum(p['atlas']['cross_chart_histogram'].values())==1600

def test_curvature_and_anomaly_boundaries():
    p=packet();W=sp.Matrix(p['curvature']['triangle']['loop']).applyfunc(sp.sympify)
    assert W.T*W==sp.eye(3) and W.det()==1 and sp.trace(W)==-1
    q,h=sp.symbols('q h');c={k:sp.sympify(v)for k,v in p['chiral']['supplied_SM_plus_nu_charge_family'].items()}
    assert sp.expand(sum(n*c[k]**3 for k,n in [('Q',6),('u_conjugate',3),('d_conjugate',3),('L',2),('e_conjugate',1),('nu_conjugate',1)]))==0
    assert p['chiral']['commuting_nonabelian_weak_gauge_factor']is False
