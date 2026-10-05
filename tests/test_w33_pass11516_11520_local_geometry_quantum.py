"""Independent algebra, normal-direction and Choi checks for11516-11520."""
import json,hashlib,sys
from pathlib import Path
import numpy as np
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11516_11520_local_geometry_quantum as P

def load():return json.loads(P.OUT.read_text())

def test_inputs_and_five_scopes():
    d=load();assert d['status']=='PASS' and d['passes']==list(range(11516,11521))
    for name,sha in d['source_sha256'].items():assert hashlib.sha256(json.dumps(P.read(name),sort_keys=True,separators=(',',':')).encode()).hexdigest()==sha
    assert all(d[k]['status']=='PASS' for k in ['vacuum','flavor','metric','ward','quantum'])

def test_full_native_gradient_and_directional_replay():
    import w33_pass11438_finite_native_model as F
    d=load()['vacuum'];x=np.array(d['coordinates']);v=np.array(d['minimum_direction'])
    with F.native_context():
        val,g=F.fun(x);assert abs(val-d['value'])<1e-10 and np.linalg.norm(g)<1e-8
        row=d['directional_controls'][0];step=row['step'];up=F.fun(x+step*v)[0];down=F.fun(x-step*v)[0]
    assert abs(up-row['plus_energy'])<1e-10 and abs(down-row['minus_energy'])<1e-10
    assert len(d['full_Hessian_spectra'])==2 and all(len(row['eigenvalues'])==324 for row in d['full_Hessian_spectra'])
    assert d['normal_dimension']+d['gauge_rank']==324

def test_actual_adjoint_one_insertion_Ward_identity():
    import w33_pass11506_11510_five_physics_targets as Q
    _,t,_=Q.tensors();y=list(map(sp.Rational,Q.read(Q.SOURCES[1])['hypercharge_diagonal']))
    for i,j,k in np.argwhere(t):assert y[i]+y[j]==-y[k]
    d=load()['flavor'];assert d['charge_blocks']['down']['single_fermion_insertion']==d['charge_blocks']['lepton']['single_fermion_insertion']

def test_two_adjoint_split_and_engineered_targets():
    d=load()['flavor'];down=sp.Rational(d['charge_blocks']['down']['one_on_each_fermion']);lep=sp.Rational(d['charge_blocks']['lepton']['one_on_each_fermion'])
    assert lep/down==-9
    F0=sp.Matrix(d['exact_F0']);F2=sp.Matrix(d['exact_F2'])
    assert F0+down*F2==sp.diag(*map(sp.Rational,d['exact_target_down']))
    assert F0+lep*F2==sp.diag(*map(sp.Rational,d['exact_target_lepton']))

def test_metric_blindness_and_two_hop_rank():
    d=load()['metric'];src=P.read(P.SOURCES[1])['metric'];h=np.array(src['displacements_times80'],int)
    assert np.all(h*h==h[:,0,None]**2)
    rows=lambda v:sp.Matrix([[int(x[0]**2),int(x[1]**2),int(x[2]**2),2*int(x[0]*x[1]),2*int(x[0]*x[2]),2*int(x[1]*x[2])] for x in v])
    a=rows(h);assert a.rank()==4 and all(a*sp.Matrix(v)==sp.zeros(160,1) for v in d['exact_blind_directions'])
    allh=np.vstack([h,d['two_hop_displacements_times80']]);w=rows(allh)[d['rank_witness_indices'],:]
    assert w.det()==sp.Rational(d['rank_witness_determinant'])!=0

def test_two_hop_vectors_are_actual_oriented_paths():
    d=load()['metric'];src=P.read(P.SOURCES[1])['metric'];adj=[{} for _ in range(80)]
    for (a,b),v in zip(src['edges'],src['displacements_times80']):adj[a][b]=np.array(v);adj[b][a]=-np.array(v)
    for (a,c),v in zip(d['two_hop_edges'],d['two_hop_displacements_times80']):
        assert any(np.array_equal(adj[a][b]+adj[b][c],v) for b in adj[a] if c in adj[b])

def test_metric_observability_in_a_second_basis():
    d=load()['metric'];src=P.read(P.SOURCES[1])['metric'];S=np.array([[1,1,0],[0,1,1],[1,0,2]],int)
    assert round(np.linalg.det(S))!=0
    rows=lambda h:sp.Matrix([[int(v[0]**2),int(v[1]**2),int(v[2]**2),2*int(v[0]*v[1]),2*int(v[0]*v[2]),2*int(v[1]*v[2])] for v in h])
    h=np.array(src['displacements_times80']);a=np.vstack([h,d['two_hop_displacements_times80']])
    assert rows(h@S).rank()==4 and rows(a@S).rank()==6

def test_flux_winding_base_second_jet():
    z=sp.symbols('z0:6',real=True);S=[sp.eye(3),sp.diag(1,-1,0),sp.diag(1,1,-2)]
    for i,j in [(0,1),(0,2),(1,2)]:s=sp.zeros(3);s[i,j]=s[j,i]=1;S.append(s)
    Z=sum((x*s for x,s in zip(z,S)),sp.zeros(3));v=1+sp.trace(Z)/2+sp.trace(Z)**2/8;iv=1-sp.trace(Z)/2+sp.trace(Z)**2/8
    potential=v+sp.Rational(3,2)*iv+v*sp.trace(sp.eye(3)-Z+Z*Z/2)/2
    zero=dict.fromkeys(z,0);assert all(sp.diff(potential,x).subs(zero)==0 for x in z)
    H=sp.hessian(potential,z).subs(zero);assert H==sp.diag(6,1,3,1,1,1)

def test_flux_winding_volume_and_mode_bounds():
    rho=sp.symbols('rho',positive=True);V=rho**3+sp.Rational(3,2)/rho**3+sp.Rational(3,2)*rho
    assert sp.expand(2*rho**4*sp.diff(V,rho)/3-(2*rho**6+rho**4-3))==0
    lam=sp.symbols('lambda',nonnegative=True);r=3*lam/((4+lam)*(1+lam))
    assert sp.cancel(sp.Rational(1,3)-r-(lam-2)**2/(3*(lam+1)*(lam+4)))==0
    d=load()['metric']['flux_winding_control'];assert all(sp.Rational(x)>0 for x in d['exact_Hessian_lower_diagonal'])
    assert sp.Rational(d['exact_counterterm_lower_bound'])>0

def test_projector_pair_identity_independent():
    rng=np.random.default_rng(11519);a=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8));u,_=np.linalg.qr(a);P0=u[:,:4]@u[:,:4].conj().T
    a=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8));K=a-a.conj().T;dot=K@P0-P0@K;transport=dot@P0-P0@dot
    assert np.linalg.norm(transport@P0-P0@transport-dot)<1e-12
    j=np.empty((2,2))
    for x in range(2):
        for y in range(2):j[x,y]=2*np.trace(transport[4*x:4*x+4,4*y:4*y+4]@P0[4*y:4*y+4,4*x:4*x+4]).real
    assert np.linalg.norm(j+j.T)<1e-12
    assert np.max(abs(j.sum(axis=1)-[np.trace(dot[4*x:4*x+4,4*x:4*x+4]).real for x in range(2)]))<1e-12

def test_uniform_free_gap_and_current_replay():
    a=sp.symbols('a0:4',nonnegative=True);gap=sum(2*x-x*x for x in a)+(sum(a)-1)**2
    assert sp.expand(gap-1-2*sum(a[i]*a[j] for i in range(4) for j in range(i+1,4)))==0
    d=load()['ward'];J=np.array(d['pair_current']);assert np.linalg.norm(J+J.T)<1e-10
    assert np.linalg.norm(J.sum(axis=1)-d['anomaly_difference'])<1e-10
    assert d['uniform_bound']['gap']>0 and d['relative_covariance_error']<1e-10

def test_all_native_logical_entries_and_Choi_positivity():
    import w33_pass11476_11480_cubic_geometry_correlated as Q
    d=load()['quantum']
    old=next(r for r in P.read(P.SOURCES[2])['recovery']['rows'] if r['p']==.03 and r['relay']==0)
    assert np.linalg.norm(np.array(d['results'][0]['rows'][0]['logical_transfer'])-old['decoded_transfer'])<1e-9
    for case in d['results']:
        for r in case['rows']:
            T=np.array(r['logical_transfer']);J=P.dec(r['logical_Choi'])
            assert np.linalg.norm(J-Q.transfer_choi(T))<1e-10 and min(np.linalg.eigvalsh(J))>-1e-9
            assert abs(np.trace(T)/16-r['core_entanglement_fidelity'])<1e-10
            assert 0<=r['complete_flagged_entanglement_fidelity']<=r['correct_sector_probability']+1e-10

def test_effect_weighted_fidelity_against_explicit_Kraus():
    # Independent Choi/Kraus comparison, including a non-diagonal complex effect.
    from scipy.linalg import expm
    p=.13;k0=np.diag([1,np.sqrt(1-p)]);k1=np.array([[0,np.sqrt(p)],[0,0]])
    kraus=[np.kron(a,b) for a in [k0,k1] for b in [k0,k1]]
    J=sum(np.outer(k.T.ravel(),k.T.ravel().conj()) for k in kraus)/4
    U=expm(1j*np.array([[.3,.2j],[-.2j,-.1]]));e=U@np.diag([.4,.9])@U.conj().T;E=np.kron(e,e)
    vec=E.T.reshape(16)/2
    assert abs(np.vdot(vec,J@vec)-sum(abs(np.trace(E@k))**2 for k in kraus)/16)<1e-12
    assert abs(np.vdot((np.eye(4)-E).T.ravel()/2,J@((np.eye(4)-E).T.ravel()/2))-sum(abs(np.trace((np.eye(4)-E)@k))**2 for k in kraus)/16)<1e-12

def test_complete_cell_instrument_and_fault_layer():
    q=load()['quantum'];src=P.read(P.SOURCES[0])['recovery'];K=P.dec(src['logical_cell_subchannel']);E=np.kron(K.conj().T@K,K.conj().T@K);v=src['verification_flip']
    for case in q['results']:
        for row in case['rows']:
            J=P.dec(row['logical_Choi']);fe=lambda x:np.vdot(x.T.ravel()/2,J@(x.T.ravel()/2)).real
            assert abs((1-v)*((1-v)*fe(E)+v*fe(np.eye(4)-E))-row['complete_flagged_entanglement_fidelity'])<1e-12
        assert case['roundoff_seven_channel_diamond_bound']<1e-10
    assert q['results'][1]['rows'][0]['complete_flagged_entanglement_fidelity']<q['results'][0]['rows'][0]['complete_flagged_entanglement_fidelity']

def test_actual_fidelity_optimal_decoder_decisions():
    case=load()['quantum']['results'][0];scores=np.array(case['complete_fidelity_branch_scores']);row=case['rows'][2]
    choices=np.array(row['fidelity_optimal_report_to_logical_syndrome']);logical=choices//4096;s=(choices%4096)//64;t=choices%64
    report_a=np.repeat(np.arange(64),64);report_b=np.tile(np.arange(64),64);error=row['bit_error']
    bitdistance=lambda a,b:np.array([int(x^y).bit_count() for x,y in zip(a,b)])
    da=bitdistance(report_a,s);db=bitdistance(report_b,t)
    probability=error**(da+db)*(1-error)**(12-da-db)
    value=np.sum(scores[logical,s,t]*probability)
    assert abs(value-row['fidelity_optimal_complete_flagged_fidelity'])<1e-10
    assert value>=row['complete_flagged_entanglement_fidelity']-1e-10
    for index in [0,17,511,1027,4095]:
        a,b=divmod(index,64);ds=np.array([(a^x).bit_count() for x in range(64)]);dt=np.array([(b^x).bit_count() for x in range(64)])
        likelihood=error**(ds[:,None]+dt[None,:])*(1-error)**(12-ds[:,None]-dt[None,:])
        assert choices[index]==np.argmax(scores*likelihood)
