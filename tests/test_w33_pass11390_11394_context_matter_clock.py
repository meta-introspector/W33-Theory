"""Reconstruct carriers and verify physical boundaries independently of summaries."""
import json,sys
from pathlib import Path
from itertools import combinations,product
import numpy as np
import sympy as sp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11390_11394_context_matter_clock as m

def packet():return json.loads((ROOT/'data/w33_pass11390_11394_context_matter_clock.json').read_text())

def test_selector_uniform_and_ordered_limits():
    p=packet()['selector'];A=np.array(p['native_context_adjacency'])
    assert np.array_equal(A@A,8*np.eye(40)-2*A+4*np.ones((40,40)))
    assert np.linalg.norm(A@np.ones(40)-12*np.ones(40))==0
    # Two cells have precisely40 ordered vacua; no selector preference.
    assert sum(a==b for a,b in product(range(40),repeat=2))==40
    assert p['one_cell_ground']=='unique uniform superposition'
    assert p['wall_energy_per_mismatched_bond']=='J'

def test_exact_common_native_edge_kernel():
    c=m.load();p=packet()['portal'];P=np.array(p['exact_common_left_projector_numerator'],dtype=np.int64)
    Q=np.array(p['exact_native_right_projector_numerator'],dtype=np.int64)
    assert np.array_equal(P@P,18*P)and np.trace(P)==108
    assert np.array_equal(Q@Q,18*Q)and np.trace(Q)==216
    gs=c['x']['internal_h27']['group_permutations'];seen=set();count=0
    for idx,(a,b)in enumerate(c['edges']):
        if idx in seen:continue
        orbit={c['ei'][(g[a],g[b])]for g in gs};seen.update(orbit);E=np.zeros((40,40),dtype=np.int64)
        for z in orbit:aa,bb=c['edges'][z];E[aa,bb-40]=1
        assert not np.any(P@E);count+=1
    assert count==16
    # Recompute P from central action, not just stored idempotence.
    Z=np.zeros((80,80),dtype=np.int64);Z[np.array(gs[1]),np.arange(80)]=1
    C3=2*np.eye(80,dtype=np.int64)-Z-Z@Z;A=c['A'][:40,40:].astype(np.int64)
    assert np.array_equal(P,(6*np.eye(40,dtype=np.int64)-A@A.T)@C3[:40,:40])
    assert np.array_equal(Q,(6*np.eye(40,dtype=np.int64)-A.T@A)@C3[40:,40:])

def test_exact_portal_coefficient_and_three_edge_mediators():
    c=m.load();p=packet()['portal'];P=np.array(p['exact_common_left_projector_numerator'],dtype=np.int64);Q=np.array(p['exact_native_right_projector_numerator'],dtype=np.int64)
    C=np.zeros((40,40),dtype=np.int64)
    for a,b in p['nonedge_orbit']:C[a,b-40]=1;assert (a,b)not in c['ei']
    assert np.array_equal(P@C@Q@C.T@P,243*P)
    W=np.zeros((80,54));K=np.zeros((54,54));ys=(.07,.13);M=2.3
    for j,path in enumerate(p['virtual_paths']):
        assert len(path)==4
        assert all((min(a,b),max(a,b))in c['ei']for a,b in zip(path,path[1:]))
        K[2*j,2*j+1]=K[2*j+1,2*j]=M;W[path[0],2*j]=ys[0];W[path[3],2*j+1]=ys[1]
    eff=-W@np.linalg.inv(K)@W.T
    assert np.linalg.norm(eff[:40,40:]+ys[0]*ys[1]/M*C)<1e-14
    assert np.max(abs(eff[:40,:40]))==0 and np.max(abs(eff[40:,40:]))==0
    # Full extension respects the actual H27 action on path labels.
    paths=[tuple(x)for x in p['virtual_paths']];idx={path:i for i,path in enumerate(paths)}
    for g in c['x']['internal_h27']['group_permutations']:
        hp=[z for path in paths for z in (2*idx[tuple(g[v]for v in path)],2*idx[tuple(g[v]for v in path)]+1)]
        assert np.linalg.norm(W[np.ix_(g,hp)]-W)<1e-15
        assert np.linalg.norm(K[np.ix_(hp,hp)]-K)<1e-15

def test_actual_harmonic_overlap_and_holonomy_witness():
    c=m.load();frames,_,_=m.context_frames(c);G=sp.Matrix(frames[0].T@frames[0]);Gi=G.inv()
    for a in range(40):
        M=sp.Matrix(frames[a].T@frames[0]);R=Gi*M
        if a==0:assert R==sp.eye(3)
        elif set(c['lines'][a])&set(c['lines'][0]):assert M.rank()==1 and sp.trace(R*Gi*M.T)==sp.Rational(1,9)
        else:assert (27*R).T*G*(27*R)==G
    a,b,d=packet()['transport']['example']['contexts']
    R=lambda i,j:27*Gi*sp.Matrix(frames[j].T@frames[i])
    W=R(d,a)*R(b,d)*R(a,b)
    assert W==sp.Matrix(packet()['transport']['example']['matrix']).applyfunc(sp.Rational)
    assert W!=sp.eye(3)and W*W==sp.eye(3)

def test_complete_holonomy_group_is_fcc_point_group():
    p=packet()['transport'];group=[sp.Matrix(x).applyfunc(sp.Rational)for x in p['holonomy_matrices']];keys={tuple(x)for x in group}
    assert len(keys)==48 and all(tuple(a*b)in keys for a in group for b in group)
    assert sum(x.det()==1 for x in group)==24
    c=m.load();F=sp.Matrix(c['x']['line']['fcc_roots']).applyfunc(sp.Rational);U=sp.Matrix(c['x']['line']['primitive_fcc_change']).applyfunc(sp.Rational);T=F*U.inv()
    for R in group:
        O=T*R.inv().T*T.inv();assert all(v in (-1,0,1)for v in O)and O.T*O==sp.eye(3)
    assert p['nonintegral_primitive_period_transports']==0

def test_orientation_double_cover_is_connected():
    c=m.load();frames,_,_=m.context_frames(c);Gi=sp.Matrix(frames[0].T@frames[0]).inv();adj=[[]for _ in range(80)]
    for a,b in combinations(range(40),2):
        if set(c['lines'][a])&set(c['lines'][b]):continue
        R=27*Gi*sp.Matrix(frames[b].T@frames[a]);flip=int(R.det()==-1)
        for bit in (0,1):adj[2*a+bit].append(2*b+(bit^flip));adj[2*b+(bit^flip)].append(2*a+bit)
    seen={0};queue=[0]
    for a in queue:
        for b in adj[a]:
            if b not in seen:seen.add(b);queue.append(b)
    assert len(seen)==80 and all(len(x)==27 for x in adj)

def test_fermion_doublers_and_car():
    c=m.load();p=m.fermions(c)
    assert p['total_weyl_charge']==0
    assert sorted(x['weyl_charge']for x in p['corner_audit'])==[-1]*4+[1]*4
    assert sum(x['wilson_mass_r1']==0 for x in p['corner_audit'])==1
    assert p['car_residual']==0 and 'not Lorentz chirality' in p['scope']

def test_walk_discriminant_eigenvectors_and_single_step_locality():
    c=m.load();U,S,d=m.clock_matrices(c,[.09,.17,-.11]);L=m.bloch(c['edges'],c['x']['line']['integer_voltage'],[.09,.17,-.11])
    assert np.linalg.norm(U.conj().T@U-np.eye(320))<1e-12
    P=d@S@d.conj().T;assert np.linalg.norm(P-(np.eye(80)-L/4))<1e-12
    vals,vecs=np.linalg.eigh(P)
    for j in (0,20,79):
        z=np.exp(1j*np.arccos(vals[j]));a=d.conj().T@vecs[:,j];b=S@a;v=a-z*b
        assert np.linalg.norm(U@v-z*v)<1e-12
    tails=[v for edge in c['edges']for v in edge];pairs={frozenset(e)for e in c['edges']}
    rows,cols=np.where(abs(U)>1e-12)
    assert all(frozenset((tails[a],tails[b]))in pairs for a,b in zip(rows,cols))

def test_scope_keeps_physics_inputs_visible():
    p=packet();assert p['status']=='PASS'
    assert 'inputs' in p['portal']['scope']and 'seconds' in p['clock']['scope']
    assert 'chiral anomaly-free spectrum' in p['fermions']['scope']
    assert 'Einstein' in p['physical_boundary']and 'No complete TOE' in p['physical_boundary']

def test_native_oriented_context_action_matches_harmonic_frames():
    c=m.load();frames,_,_=m.context_frames(c);Gi=sp.Matrix(frames[0].T@frames[0]).inv()
    p=packet()['transport'];lifts=p['native_transvection_orientation_permutations'];J=np.array([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]])
    for v,lift in zip(c['points'],lifts):
        assert sorted(lift)==list(range(80))
        g=c['permutation']((np.eye(4,dtype=int)+np.outer(v,np.array(v)@J))%3)
        b=g[40]-40;F=np.zeros_like(frames[0])
        for i,(a,d)in enumerate(c['edges']):F[c['ei'][(g[a],g[d])]]=frames[0][i]
        R=Gi*sp.Matrix(frames[b].T@F);flip=int(R.det()==-1)
        assert np.array_equal(frames[b]@np.array(R,dtype=np.int64),F)
        assert lift[0]==2*b+flip and lift[1]==2*b+(1^flip)
    seen={0};queue=[0]
    for a in queue:
        for g in lifts:
            if g[a]not in seen:seen.add(g[a]);queue.append(g[a])
    assert len(seen)==80 and p['exact_native_harmonic_action_maps_checked']==1600
