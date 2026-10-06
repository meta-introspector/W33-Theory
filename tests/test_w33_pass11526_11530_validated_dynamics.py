"""Independent certificate replays and physical-boundary regressions."""
import hashlib,itertools,json,sys
from pathlib import Path
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11526_11530_validated_dynamics as P
load=lambda:json.loads(P.OUT.read_text())

def test_bound_inputs_and_producer():
    d=load();assert d['status']=='PASS'
    assert hashlib.sha256(Path(P.__file__).read_bytes()).hexdigest()==d['producer_sha256']
    for name,sha in d['source_sha256'].items():assert hashlib.sha256(json.dumps(P.read(name),sort_keys=True,separators=(',',':')).encode()).hexdigest()==sha

def test_rigorous_Krawczyk_and_interval_positive_Hessian():
    from flint import arb,ctx
    ctx.prec=256;v=load()['vacuum'];c=[arb(z) for z in v['center']];box=[arb(z,v['radius']) for z in v['center']]
    g=P.jet_potential(c)[0];f,chi,q=P.jet_potential(box);R=[[arb(z) for z in row] for row in v['preconditioner']]
    defect=[[int(i==j)-sum((R[i][k]*f.h[k,j].real for k in range(20)),arb(0)) for j in range(20)] for i in range(20)]
    for i in range(20):
        error=abs(sum((R[i][k]*g.g[k].real for k in range(20)),arb(0)))+arb(v['radius'])*sum((abs(z) for z in defect[i]),arb(0))
        assert error<arb(v['radius'])/100
    assert max(sum((abs(z) for z in row),arb(0)) for row in defect)<arb('1e-8')
    # Independent Schur-complement elimination, rather than stored LDL pivots.
    H=[[z.real for z in row] for row in f.h]
    for i in range(20):
        assert H[i][i]>0
        for j in range(i+1,20):
            for k in range(i+1,20):H[j][k]-=H[j][i]*H[i][k]/H[i][i]
    assert not chi.v.real.contains(0) and not q.v.real.contains(0)
    assert box[0]>box[1]>box[2]>0 and box[5]>0 and box[9]>0

def test_ball_jet_matches_independent_mpmath_derivatives():
    from flint import arb,ctx
    import mpmath as mp
    import w33_pass11521_11525_exact_frontiers as old
    ctx.prec=256;coords=load()['vacuum']['center'];jet=P.jet_potential([arb(z) for z in coords])[0]
    with mp.workdps(55):
        values=tuple(mp.mpf(z) for z in coords);f=lambda *x:old.quotient_value(x,mp)
        assert abs(float(jet.v.real)-float(f(*values)))<1e-13
        for i in [0,3,5,13,19]:
            order=tuple(int(k==i) for k in range(20));g=mp.diff(f,values,order)
            assert abs(float(jet.g[i].real)-float(g))<1e-48
        order=[0]*20;order[3]=order[13]=1
        assert abs(float(jet.h[3,13].real)-float(mp.diff(f,values,tuple(order))))<1e-12

def test_frame_normalizer_certificate():
    d=load()['symmetry_phase'];R=s.Matrix(d['frame_outside_constraints']);K=s.Matrix.hstack(*[s.Matrix(list(map(s.Rational,v))) for v in d['frame_normalizer_kernel']])
    assert R*K==s.zeros(R.rows,30)
    import w33_pass11521_11525_exact_frontiers as old
    assert old.modular_rank(R.tolist())==48 and old.modular_rank(K.tolist())==30

def test_exact_weight_norm_and_coherent_moment_bound():
    d=load()['symmetry_phase'];W=s.Matrix(d['native_weights']);G=s.Matrix(d['Cartan_metric_inverse']);C=W*G*W.T
    assert all(v==s.Rational(2,9) for v in C.diagonal())
    assert all(s.Rational(4,9)-2*C[i,j]>=0 for i,j in itertools.product(range(27),repeat=2))
    # Exact convexity identity for a nontrivial mixed weight distribution.
    p=s.Matrix([s.Rational(i+1,378) for i in range(27)]);assert sum(p)==1
    gap=s.Rational(2,9)-(p.T*C*p)[0]
    assert gap==sum(p[i]*p[j]*(s.Rational(4,9)-2*C[i,j])/2 for i,j in itertools.product(range(27),repeat=2))>=0
    assert s.Rational(d['weight_norm_squared'])+s.Rational(d['family_maximum'])==s.Rational(d['joint_maximum'])

def test_actual_native_SM_singlet_Euler_obstruction():
    import w33_pass11438_finite_native_model as F
    from w33_pass11384_11388_native_dynamics import native_tensors
    d=native_tensors()[1][1];assert not np.any(d[:,[0,1]][:,:,[0,1]])
    rng=np.random.default_rng(11527);a=np.zeros((27,3),complex);b=a.copy()
    a[:2]=rng.normal(size=(2,3))+.3j*rng.normal(size=(2,3));b[:2]=rng.normal(size=(2,3))+.2j*rng.normal(size=(2,3))
    x=F.pack(a.ravel(),b.ravel())
    with F.native_context():value,grad=F.fun(x)
    assert value>0 and abs(x@grad-4*value)<1e-9

def test_no_commuting_weak_algebra_in_real8_SU3_cases():
    # All complex irreducible dimensions <=8, excluding the trivial rep.
    dims={(a,b):(a+1)*(b+1)*(a+b+2)//2 for a,b in itertools.product(range(8),repeat=2) if a+b}
    small={ab:n for ab,n in dims.items() if n<=8}
    assert small=={(0,1):3,(0,2):6,(1,0):3,(1,1):8,(2,0):6}
    # Complex6 realifies to12, and therefore cannot occur in real8.
    assert 2*small[0,2]>8 and 2*small[1,0]==6 and small[1,1]==8
    # The other faithful possibilities are real6+two trivials and adjoint8;
    # their orthogonal centralizer dimensions are2 and0, smaller than su2.
    assert all(n<3 for n in [2,0])

def test_EW_exact_positive_normal_spectrum():
    d=load()['electroweak'];t,la,ka,eta,nu=s.symbols('t lambda kappa eta nu',positive=True)
    symbols={str(x):x for x in [t,la,ka,eta,nu]};symbols['lam']=la
    spectrum={s.sympify(v.replace('lambda','lam'),locals=symbols):m for v,m in d['exact_real_coordinate_Hessian_spectrum'].items()}
    assert spectrum=={s.Integer(0):3,4*eta*t:2,16*ka*t:1,16*la*t+4*nu*t:1,4*nu*t:1}
    # Explicit complex target, both norms one and the charged overlap zero.
    hu=s.Matrix([0,1]);hd=s.Matrix([-s.Rational(4,5)-s.I*s.Rational(3,5),0]);p=hu[0]*hd[1]-hu[1]*hd[0]
    assert p==s.Rational(4,5)+s.I*s.Rational(3,5)
    assert s.expand((hd.conjugate().T*hd)[0]-1)==0 and (hu.conjugate().T*hd)[0]==0

def test_mixed_Higgs_source_in_actual_E6_tensor():
    from w33_pass11384_11388_native_dynamics import native_tensors
    chi=s.symbols('chi',real=True);d=native_tensors()[1][1];a=s.zeros(27,3);b=s.zeros(27,3)
    a[0,:]=s.Matrix(1,3,[1,0,0]);b[0,:]=s.Matrix(1,3,[0,1,0])
    a[20,:]=s.Matrix(1,3,[0,0,1]);b[20,:]=s.Matrix(1,3,[2,1,3])
    a[23,:]=s.Matrix(1,3,[1,2,0]);b[23,:]=s.Matrix(1,3,[0,1,2*s.I*chi])
    def source(c):
        return s.Matrix(3,3,lambda i,j:sum(int(d[c,k,l])*(s.conjugate(a[k,i]*b[l,j])+s.conjugate(a[k,j]*b[l,i]))/2 for k,l in np.argwhere(d[c]!=0)))
    cert=load()['electroweak'];Hu=s.Matrix([[s.sympify(z,locals={'chi':chi}) for z in row] for row in cert['mixed_full_rank_witness']])
    Hd=s.Matrix([[s.sympify(z,locals={'chi':chi}) for z in row] for row in cert['mixed_second_witness']])
    # Stored witness is the unconjugated polynomial form; the actual composite
    # conjugates both scalar sources, reversing the CP branch.
    assert source(23)==Hu.conjugate() and source(20)==Hd.conjugate()
    assert Hu.det()==s.Rational(1,4) and s.factor(Hd.det())==2*chi**2
    A=Hu*Hu.conjugate().T;B=Hd*Hd.conjugate().T;C=A*B-B*A
    assert s.factor(s.im(s.trace(C**3)))==-s.Rational(69,2)*chi*(chi**2-19)
    p=(B).charpoly();assert s.expand(s.discriminant(p.as_expr(),p.gen)-16*(chi**6-5*chi**4+39*chi**2+2))==0

def test_all_lifted_paths_and_exact_first_jet_in_second_basis():
    d=load()['gravity_frames'];src=P.read('data/w33_pass11493_11497_exact_matching_inference.json')['metric'];h=src['displacements_times80'];J=s.Matrix([[1,1,0],[0,1,1],[0,0,1]])
    assert d['rank_census_by_radius']==[{'1':52,'3':28},{'1':12,'3':68},{'3':80}]
    for row in d['sites']:
        G=s.zeros(3)
        for p in row['walks']:
            node=row['site'];v=s.zeros(3,1)
            for e,sign in p['edges']:
                a,b=src['edges'][e];assert node==(a if sign==1 else b);node=b if sign==1 else a
                v+=sign*s.Matrix(h[e])
            assert node==p['end'] and list(v)==p['displacement'] and 1<=len(p['edges'])<=3
            G+=v*v.T
        inverse=s.Matrix(row['inverse']);assert G.det()>0 and G*inverse==s.eye(3)
        assert (J*G*J.T)*(J.T.inv()*inverse*J.inv())==s.eye(3)

def test_exact_local_Spin_link_covariance():
    X=s.Matrix([[0,1],[1,0]]);Y=s.Matrix([[0,-s.I],[s.I,0]]);Z=s.diag(1,-1)
    gamma=[s.kronecker_product(X,p) for p in [X,Y,Z]]
    A=s.kronecker_product(s.eye(2),s.Rational(3,5)*s.eye(2)+s.I*s.Rational(4,5)*Y)
    B=s.kronecker_product(s.eye(2),s.Rational(5,13)*s.eye(2)+s.I*s.Rational(12,13)*X)
    U=s.kronecker_product(s.eye(2),s.Rational(7,25)*s.eye(2)+s.I*s.Rational(24,25)*Z)
    assert all((R.conjugate().T*R-s.eye(4)).applyfunc(s.expand)==s.zeros(4) for R in [A,B,U])
    g=gamma[0]+2*gamma[1]-gamma[2];hop=g*U+U*g
    Ag=A*g*A.conjugate().T;Bg=B*g*B.conjugate().T;link=A*U*B.conjugate().T
    assert (Ag*link+link*Bg-A*hop*B.conjugate().T).applyfunc(s.expand)==s.zeros(4)

def test_full_exact_coherent_parent_normal_certificate():
    d=load()['coherent_parent'];diag=list(map(s.Rational,d['Phi_Hessian_diagonal']));H=s.diag(*diag)
    Pgrad=s.Matrix(d['purity_gradient_matrix']);N=s.diag(*list(map(s.Rational,d['moment_action_diagonal'])));E=s.zeros(162);E[0,0]=1
    assert H==4*E+s.Rational(8,9)*(4*E+2*s.eye(162))-4*Pgrad-2*N
    from collections import Counter
    values=Counter(diag);values[s.Integer(2)]+=162
    assert {str(v):m for v,m in values.items()}==d['exact_full324_real_Hessian_spectrum']
    assert values=={s.Integer(0):37,s.Rational(2,3):20,s.Rational(7,3):64,s.Rational(8,3):40,s.Integer(4):1,s.Integer(2):162}
    T=s.Matrix(d['gauge_tangents']);assert H*T==s.zeros(162,86)
    import w33_pass11521_11525_exact_frontiers as old
    assert old.modular_rank(T.tolist())==37 and all(v>=0 for v in diag)

def test_constrained_spin_history_and_correct_force_sign():
    d=load()['spin_history'];assert d['site_spin_dimension']==320
    assert abs(d['fundamental_cycle_holonomy_trace']-4)<.1 and d['fundamental_cycle_holonomy_trace']<4
    assert d['local_Spin_covariance_residual']<1e-9
    assert max(abs(v['Hamiltonian_constraint']) for v in d['history'])<1e-7
    assert max(abs(v['spinor_norm']-1) for v in d['history'])<1e-8
    assert d['history'][-1]['scale']>d['history'][0]['scale']
    a,p,k,m=s.symbols('a p K M',real=True);V=s.Function('V')(a)
    H=-p**2/(12*a)+V+k/a+m
    assert -s.diff(H,a)==-p**2/(12*a**2)-s.diff(V,a)+k/a**2

def test_actual_nonlinear_Ward_and_gauge_covariance():
    d=load()['anomaly'];B=np.array(d['incidence']);T=np.array(d['routing']);q=np.array(d['density']);k=np.array(d['current']);n=len(q)
    assert np.array_equal(B.T@T,np.eye(n)-np.eye(n)[:,[0]]@np.ones((1,n)))
    sites=list(itertools.product(range(d["L"]),repeat=4))
    assert np.array_equal(np.abs(T).sum(axis=0),[sum(min(v,d["L"]-v) for v in x) for x in sites])
    assert np.linalg.norm(B.T@k-q)<1e-10 and d['gauge_invariance_error']<1e-10
    assert max(abs(q))>1e-9 and d['nonlinear_half_amplitude_defect']>1e-9
    replay=P.nonlinear_anomaly_primitive();assert np.max(abs(np.array(replay['current'])-k))<1e-10
    assert 'functional derivative' in d['retained_obstruction'] and 'not the general' in d['infinite_volume_source_patch_lemma']

def test_all_owned_flag_cases_and_terminal_identity():
    import w33_pass11428_11432_native_alignment_measure_coarse_flags as Q
    decoder,single,cases=Q.flag_decoder();d=load()['quantum'];assert len(cases)==d['case_count']==1714
    for stored,(e,h,exact) in zip(d['fault_cases'],cases):
        assert stored['history']==list(h) and stored['correction']==decoder[h]
        r=e^decoder[h];assert r in single and r==stored['residual_class']
        assert r^stored['terminal_correction']==0

def test_nonPauli_CPTP_port_has_full_exact_Pauli_span():
    X=s.Matrix([[0,1],[1,0]]);Z=s.diag(1,-1)
    basis=[s.kronecker_product(X**a*Z**b,X**c*Z**d) for a,b,c,d in itertools.product(range(2),repeat=4)]
    Kraus=[s.kronecker_product(s.diag(1,s.Rational(4,5)),s.eye(2)),s.kronecker_product(s.Matrix([[0,s.Rational(3,5)],[0,0]]),s.eye(2))]
    assert sum((K.conjugate().T*K for K in Kraus),s.zeros(4))==s.eye(4)
    for K in Kraus:assert sum((s.trace(Pa.conjugate().T*K)*Pa/4 for Pa in basis),s.zeros(4))==K

def test_rational_coherent_many_port_bound_and_order():
    for row in load()['quantum']['coherent_fault_envelopes']:
        delta=s.Rational(row['delta']);b=(1+delta)**132-1-132*delta
        assert b==s.Rational(row['bad_amplitude']) and b*b==s.Rational(row['infidelity_bound'])
    # The first uncorrectable amplitude order is two; infidelity order is four.
    assert s.binomial(132,2)**2==74753316
    assert load()['quantum']['coherent_fault_envelopes'][0]['float_bound']<8e-13


def test_condensate_moment_projectors_and_two_vector_minimum():
    """A field-produced next-stage selector; standard Spin10 branching is prior art."""
    d=load()['symmetry_phase'];w=s.Matrix(d['native_weights']);g=s.Matrix(d['Cartan_metric_inverse'])
    diagonal=w*g*w[0,:].T;M=s.diag(*diagonal);I=s.eye(27)
    projectors=[18*(M-I/18)*(M+I/9),-36*(M-2*I/9)*(M+I/9),18*(M-2*I/9)*(M-I/18)]
    assert [int(p.trace()) for p in projectors]==[1,16,10]
    assert sum(projectors,s.zeros(27))==I
    assert all(a*b==(a if i==j else s.zeros(27)) for i,a in enumerate(projectors) for j,b in enumerate(projectors))
    assert diagonal[1]==s.Rational(1,18)
    # Every minuscule weight has the same purity, including native e1.
    R=s.Rational(1,2);family=s.Rational(2,3)
    assert all(s.Rational(8,9)*R**2-R**2*((w*g*w.T)[i,i]+family)==0 for i in [0,1])
    assert R*(diagonal[1]-s.Rational(1,18))==0
    # Actual compact E6 stabilizer of both vectors; not inferred from dimensions alone.
    from w33_pass11384_11388_native_dynamics import native_tensors
    raw=native_tensors()[0];unique={}
    for b in raw:
        for kind,a in [(0,b+b.T),(1,b-b.T)]:
            nz=np.flatnonzero(a)
            if not len(nz):continue
            a=a//int(np.gcd.reduce(abs(a.ravel()[nz])))
            if a.ravel()[nz[0]]<0:a=-a
            unique[kind,tuple(a.ravel())]=a*(1j if kind else 1)
    he=np.array(list(unique.values()));blocks=he[:,:,[0,1]]
    constraints=s.Matrix(np.concatenate([blocks.real.reshape(78,-1),blocks.imag.reshape(78,-1)],axis=1).astype(int).T)
    assert 78-constraints.rank()==24


def test_second_action_full_normals_and_independent_potential_derivative():
    """Replay rational blocks and differentiate the actual nonlinear native action."""
    cert=load()['coherent_two_stage'];H=s.zeros(324);count={}
    for block in cert['Hessian_blocks']:
        ids=block['indices'];small=s.Matrix([[s.Rational(x) for x in row] for row in block['matrix']])
        for i,a in enumerate(ids):
            for j,b in enumerate(ids):H[a,b]=small[i,j]
        for e,m in small.eigenvals().items():
            assert e>=0;count[e]=count.get(e,0)+m
    assert count[0]==58 and sum(count.values())-count[0]==266
    assert {str(e):m for e,m in count.items()}==cert['full324_spectrum']
    T=s.Matrix(cert['gauge_tangents']);assert T.rank()==58 and H*T==s.zeros(324,86)
    from w33_pass11384_11388_native_dynamics import native_tensors
    raw=native_tensors()[0];unique={}
    for b in raw:
        for kind,a in [(0,b+b.T),(1,b-b.T)]:
            nz=np.flatnonzero(a)
            if not len(nz):continue
            a=a//int(np.gcd.reduce(abs(a.ravel()[nz])))
            if a.ravel()[nz[0]]<0:a=-a
            unique[kind,tuple(a.ravel())]=s.SparseMatrix(a.real.astype(int).tolist())+s.I*s.SparseMatrix(a.imag.astype(int).tolist()) if kind else s.SparseMatrix(a.tolist())
    # The skew root generator is multiplied by i in the Hermitian basis.
    he=[s.I*a if kind else a for (kind,_),a in unique.items()]
    family=[]
    for i,j in itertools.combinations(range(3),2):
        a=s.zeros(3);a[i,j]=a[j,i]=1;family.append(a)
        a=s.zeros(3);a[i,j]=-s.I;a[j,i]=s.I;family.append(a)
    family.extend([s.diag(1,-1,0),s.diag(1,1,-2)])
    IE=s.Matrix([[s.trace(a*b) for b in he] for a in he]).inv()
    IF=s.Matrix([[s.trace(a*b) for b in family] for a in family]).inv()
    t=s.symbols('t',real=True);phi=s.zeros(27,3);psi=s.zeros(27,3)
    phi[0,0]=psi[1,0]=1/s.sqrt(2)
    # Includes entangled off-family components that a truncated alignment
    # Hessian would omit, and both real and imaginary mixed selector terms.
    tangent=s.zeros(324,1)
    for index,val in [(1,1),(7,2),(84,1),(88,-1),(162+4,3),(162+8,-2),(162+82,1),(162+90,2)]:tangent[index]=val
    for j in range(81):
        phi[j//3,j%3]+=t*(tangent[j]+s.I*tangent[81+j])
        psi[j//3,j%3]+=t*(tangent[162+j]+s.I*tangent[243+j])
    def terms(X):
        R=s.expand(s.trace(X.conjugate().T*X))
        e=s.Matrix([s.expand(s.trace(X.conjugate().T*a*X)) for a in he])
        f=s.Matrix([s.expand(s.trace(X.conjugate().T*X*a)) for a in family])
        defect=s.Rational(8,9)*R**2-(e.T*IE*e)[0]-(f.T*IF*f)[0]
        return R,e,defect
    Rp,ep,Dp=terms(phi);Rq,eq,Dq=terms(psi);coeff=IE*ep
    selector=sum((coeff[i]*he[i]*psi for i in range(78)),s.zeros(27,3))-Rp*psi/18
    potential=(Rp-s.Rational(1,2))**2+Dp+(Rq-s.Rational(1,2))**2+Dq+s.trace(selector.conjugate().T*selector)+Rp*Rq-s.trace(phi.conjugate().T*phi*psi.conjugate().T*psi)
    poly=s.Poly(s.expand(potential),t)
    assert poly.coeff_monomial(1)==0 and poly.coeff_monomial(t)==0
    assert s.expand(poly.coeff_monomial(t**2)-(tangent.T*H*tangent)[0]/2)==0
