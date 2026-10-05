"""Independent exact reductions, physical boundary and recovery-dual checks."""
import sys,json,hashlib,itertools
from pathlib import Path
import numpy as np
import sympy as s
from functools import lru_cache
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11521_11525_exact_frontiers as P
@lru_cache(maxsize=1)
def load():return json.loads(P.OUT.read_text())

def test_binding_and_scopes():
    d=load();assert d['passes']==list(range(11521,11526))
    for path,sha in d['source_sha256'].items():assert hashlib.sha256(json.dumps(P.read(path),sort_keys=True,separators=(',',':')).encode()).hexdigest()==sha
    assert hashlib.sha256(Path(P.__file__).read_bytes()).hexdigest()==d['producer_sha256']
    assert all(d[key]['status']=='PASS' for key in ['vacuum','flavor','gravity','ward','quantum'])

def test_exact_row_penalty_against_full_native_moments():
    import w33_pass11438_finite_native_model as F
    d=load()['vacuum'];C=np.array(s.Matrix(d['row_penalty_matrix']),float)
    tri=F.TENSORS[0][0];E=np.zeros((27,3));E[tri[0],0]=1;E[tri[1],1]=1;E[tri[2],2]=tri[3]
    A=np.array([[1,1j,0],[0,2,1],[1,0,3]],complex);field=(E@A).ravel()
    mu=np.einsum('i,aij,j->a',field.conj(),F.H,field).real
    left=np.diag(A@A.conj().T).real;G=A.conj().T@A
    predicted=left@C@left+np.trace((G-np.trace(G)/3*np.eye(3))@(G-np.trace(G)/3*np.eye(3))).real
    assert abs(mu@mu-predicted)<1e-10
    assert s.Matrix(d['row_penalty_matrix'])*s.ones(3,1)==s.zeros(3,1)

def test_common_unitary_invariants_on_exact_hollow_locus():
    A=s.Matrix([[1,s.I,0],[0,1,s.I],[s.I,0,1]]);B=A+2*s.eye(3)
    U=s.Matrix([[s.Rational(3,5),s.Rational(4,5),0],[-s.Rational(4,5),s.Rational(3,5),0],[0,0,1]])
    assert U.conjugate().T*U==s.eye(3) and U.det()==1
    for D in [A,B]:
        assert (D.conjugate().T*D-(U*D).conjugate().T*(U*D)).applyfunc(s.expand)==s.zeros(3)
        assert s.expand(D.det()-(U*D).det())==0
        assert all(s.expand(a-b)==0 for a,b in zip((D*D.conjugate().T).diagonal(),((U*D)*(U*D).conjugate().T).diagonal()))
        assert len(set((D*D.conjugate().T).diagonal()))==1
    assert (A.T*B.conjugate()-(U*A).T*(U*B).conjugate()).applyfunc(s.expand)==s.zeros(3)

def test_spin8_schur_block_controls():
    d=load()['vacuum'];assert d['real_commutant_ranks']==[63]*3 and d['intertwiner_constraint_ranks']==[64]*3
    for witness in d['rank_witnesses']:assert P.modular_rank(witness['minor'],witness['prime'])==witness['rank']
    assert d['gradient_norm']<1e-8 and abs(d['value']+2.9357432396668512)<1e-9
    for block in d['normal_blocks']:
        B=np.array(block['block']);assert B.shape==(12,12) and np.linalg.norm(B-B.T)<1e-10
        assert np.linalg.norm(np.linalg.eigvalsh(B)-block['eigenvalues'])<1e-10
        assert block['Schur_residual']<1e-3
    assert d['cross_block_residual']<1e-3

def test_neutral_composite_Higgs_obstruction():
    from w33_pass11384_11388_native_dynamics import native_tensors
    _,(_,d,_)=native_tensors();data=load()['flavor'];y=list(map(s.Rational,P.read(P.SOURCES[4])['hypercharge_diagonal']))
    for a,b,c in np.argwhere(d):assert y[a]+y[b]+y[c]==0
    for h in data['Higgs_indices']:
        assert y[h]!=0
        assert all(d[h,a,b]==0 for a,b in itertools.product(data['hypercharge_zero_indices'],repeat=2))

def test_symmetric_insertion_ring_actual_charge_pairs():
    data=load()['flavor']['charge_pair_invariants'];assert s.Rational(data['sum'])==s.Rational(1,2)
    down=s.Rational(data['down_product']);lepton=s.Rational(data['lepton_product']);assert lepton/down==-9
    # A higher symmetric insertion is reduced to the same two generators.
    x,y=s.symbols('x y');assert s.expand(x**3+y**3-((x+y)**3-3*x*y*(x+y)))==0

def test_static_lapse_stationarity_no_go_exact():
    z=s.symbols('z',real=True);a,b,la,g,w0,wz=s.symbols('alpha beta Lambda gamma W0 Wz',real=True)
    V=la*s.exp(3*z/2)+a*s.exp(-3*z/2)+3*b*s.exp(z/2)/2+g*(w0+wz*z)
    laeq=s.solve(s.diff(V,z).subs(z,0),la)[0]
    E=s.expand(V.subs({z:0,la:laeq}));assert E==2*a+b+g*w0-s.Rational(2,3)*g*wz
    both=s.solve([V.subs(z,0),s.diff(V,z).subs(z,0)],[a,la]);assert s.expand(both[a]+b/2+g*w0/2-g*wz/3)==0
    assert s.Rational(load()['gravity']['inherited_energy_lower_bound'])>0
    t=s.symbols('t',positive=True)
    assert s.cancel(s.diff(s.log(1+t)-t/(1+t),t)-t/(1+t)**2)==0

def test_incidence_Dirac_zero_mode_factor():
    # Independent rectangular rank-two graph, no spectral dimension analogy.
    incidence=s.Matrix([[-1,1,0],[0,-1,1],[1,0,-1],[0,0,0]])
    D=s.zeros(7);D[:3,3:]=incidence.T;D[3:,:3]=incidence
    mass=s.symbols('mass');assert s.expand((mass*s.eye(7)+s.I*D).det()-mass*(mass**2*s.eye(3)+incidence.T*incidence).det())==0

def test_projector_vertical_curvature_identity():
    proj=s.diag(1,0);dot=s.Matrix([[0,1+s.I],[1-s.I,0]]);omega=s.Matrix([[2,3*s.I],[-3*s.I,4]])
    comm=lambda a,b:a*b-b*a
    assert s.trace(proj*comm(comm(omega,proj),dot))==-s.trace(omega*dot)
    assert proj*dot*proj==s.zeros(2)
    d=load()['ward'];assert d['exact_charge_sums']['1']==d['exact_charge_sums']['3']==0
    assert d['exact_charge_sums']['5']!=0 and all(n%2==0 for n in d['absolute_odd_charge_counts'].values())

def test_missing_measure_term_is_exact_one_form():
    # Check the chain-rule coefficient and curl independently with a polynomial k.
    x,y,t=s.symbols('x y t');k=s.Matrix([x*x+y,x*y]);S=s.integrate((s.Matrix([x,y]).T*k.subs({x:t*x,y:t*y},simultaneous=True))[0],(t,0,1))
    j=s.Matrix([s.diff(S,x),s.diff(S,y)])
    assert s.diff(j[1],x)-s.diff(j[0],y)==0
    assert s.expand(j[0]-s.integrate(t**2*(3*x*x+y*y)+t*y,(t,0,1)))==0

def test_recovery_objective_against_explicit_Kraus_composition():
    import w33_pass11476_11480_cubic_geometry_correlated as M
    pauli=[np.kron(a,b) for a,b in itertools.product(M.PAULI,repeat=2)]
    rng=np.random.default_rng(11525);N=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));N/=8
    R=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));R/=8
    transfer=np.array([[np.trace(a@N@b@N.conj().T).real/4 for b in pauli] for a in pauli])
    Q=P.recovery_objective(transfer);E,v=P.cell_effect()
    direct=(1-v)*((1-v)*abs(np.trace(E@R@N))**2+v*abs(np.trace((np.eye(4)-E)@R@N))**2)/16
    assert abs(np.vdot(R.ravel(),Q@R.ravel()).real-direct)<1e-13

def test_exact_recovery_primal_dual_certificate():
    c=load()['quantum']['exact_certificate'];Q=s.Matrix(c['objective']);X=s.Matrix(c['primal']);Y=s.Matrix(c['dual'])
    # Rank-one Pauli primal: exact PSD and TP, no tolerance deletion.
    assert X==X.conjugate().T and X.rank()==1 and s.trace(X)==4
    assert s.Matrix(4,4,lambda b,d:sum(X[4*a+b,4*a+d] for a in range(4)))==s.eye(4)
    pivots=P.exact_positive_pivots(s.kronecker_product(s.eye(4),Y)-Q)
    assert pivots==list(map(s.Rational,c['dual_pivots']))
    lower=s.re(s.trace(Q*X));upper=s.trace(Y)
    assert lower==s.Rational(c['lower']) and upper==s.Rational(c['upper']) and upper>=lower
    assert upper-lower<s.Rational(1,1000000)

def test_all_branch_exact_Bell_dual_bound():
    d=load()['quantum'];den=int(d['exact_Bell_denominator']);lower=upper=0
    assert len(d['exact_all_branch_Bell_objectives'])==4096
    for row in d['exact_all_branch_Bell_objectives']:
        re=np.array(row['real_numerators'],dtype=np.int64);im=np.array(row['imag_numerators'],dtype=np.int64)
        assert np.array_equal(re,re.T) and np.array_equal(im,-im.T)
        diag=np.diag(re);radius=(abs(re)+abs(im)).sum(axis=1)-abs(diag)
        lower+=4*int(diag.max());upper+=4*int((diag+radius).max())
    assert s.Rational(lower,den)==s.Rational(d['exact_all_branch_Pauli_lower'])
    assert s.Rational(upper,den)==s.Rational(d['exact_all_branch_coherent_upper'])
    assert upper>=lower

def test_single_ancilla_hook_is_logical_Z():
    # Binary Steane argument independent of the stored128x128 contraction.
    hook=load()['quantum']['extraction_hook'];cols=[7-i for i in range(7)]
    labels=[cols[i] for i in hook['data_Z_hook']]
    assert labels[0]^labels[1]==hook['X_check_syndrome']
    residual=hook['data_Z_hook']+[hook['standard_Z_correction']]
    assert len(set(residual))==3 and np.bitwise_xor.reduce([cols[i] for i in residual])==0
    H=np.array([[(j>>k)&1 for j in range(1,8)] for k in range(3)])
    stabilizer_weights={int(((np.array(v)@H)%2).sum()) for v in itertools.product(range(2),repeat=3)}
    assert stabilizer_weights=={0,4} and len(residual) not in stabilizer_weights
    assert np.linalg.norm(P.dec(hook['decoded_logical_operator'])-np.diag([1,-1]))<1e-12

def test_flagged_circuit_corrects_entire_named_hook_family():
    d=load()['quantum']['flagged_extraction_repair'];gates=d['gates'];ancilla=d['syndrome_ancilla'];flag=d['flag_ancilla']
    assert gates[0]==gates[-1]==[flag,ancilla] and len(d['rows'])==7
    H=np.array([[((7-j)>>k)&1 for j in range(7)] for k in range(3)])
    words={tuple((np.array(v)@H)%2) for v in itertools.product(range(2),repeat=3)}
    for row in d['rows']:
        z=[0]*9;z[ancilla]=1
        for control,target in gates[row['fault_after_gate']:]:z[control]^=z[target]
        assert z[flag]==row['flag_X_minus']
        assert [i for i in range(7) if z[i]]==row['data_Z_support']
        for i in row['correction_Z_support']:z[i]^=1
        assert tuple(z[:7]) in words

def test_explicit_native_Clifford_lift():
    d=load()['gravity']['explicit_flat_spin_lift'];gg=[s.Matrix(x) for x in d['spatial_gamma']];time=s.Matrix(d['Lorentzian_time_gamma'])
    assert time*time==-s.eye(4)
    for i,j in itertools.product(range(3),repeat=2):assert gg[i]*gg[j]+gg[j]*gg[i]==(2 if i==j else 0)*s.eye(4)
    assert all(time*g+g*time==s.zeros(4) for g in gg)
    h=P.read(P.SOURCES[1])['metric']['displacements_times80']
    for row in h:
        operator=sum((int(a)*g for a,g in zip(row,gg)),s.zeros(4));assert (operator*operator-sum(a*a for a in row)*s.eye(4)).applyfunc(s.expand)==s.zeros(4)
    assert d['numerical_square_replay_error']<1e-8 and d['zero_modes']==4*82

def test_all_orbit_exact_recovery_duals():
    d=load()['quantum']['coherent_recovery_orbit_bound'];assert d['orbit_count']==66
    import w33_pass11476_11480_cubic_geometry_correlated as M
    pauli=[np.kron(a,b) for a,b in itertools.product(M.PAULI,repeat=2)]
    lower=upper=s.Integer(0);seen=[]
    for row in d['certificates']:
        Q=s.Matrix(row['objective']);Y=s.Matrix(row['dual']);assert all(v>0 for v in P.exact_positive_pivots(s.kronecker_product(s.eye(4),Y)-Q))
        v=s.Matrix(list(P.rational_matrix(pauli[row['Pauli_choice']],0)));lo=s.re((v.conjugate().T*Q*v)[0]);up=s.trace(Y)
        assert lo==s.Rational(row['lower']) and up==s.Rational(row['upper']) and up>=lo
        assert len(row['members'])==row['multiplicity']
        assert all(P.syndrome_relation_key(i)==row['relation_key'] for i in row['members'])
        lower+=row['multiplicity']*lo;upper+=row['multiplicity']*up;seen.extend(row['members'])
    assert sorted(seen)==list(range(4096)) and lower==s.Rational(d['exact_Pauli_lower']) and upper==s.Rational(d['exact_CPTP_upper'])
    assert upper-lower<s.Rational(1,100000) and d['numerical_orbit_replay_error']<1e-10

def test_syndrome_kernel_really_classifies_GL3_orbits():
    matrices=[]
    for entries in itertools.product(range(2),repeat=9):
        g=np.array(entries).reshape(3,3)
        if P.modular_rank(g,2)==3:matrices.append(g)
    assert len(matrices)==168
    d=load()['quantum']['coherent_recovery_orbit_bound']
    def image(g,x):
        bits=(g@np.array([(x>>i)&1 for i in range(3)]))%2
        return sum(int(a)<<i for i,a in enumerate(bits))
    for row in d['certificates']:
        left,right=divmod(row['representative'],64);vectors=[left//8,left%8,right//8,right%8];orbit=set()
        for g in matrices:
            a,b,c,e=[image(g,x) for x in vectors];orbit.add(64*(8*a+b)+8*c+e)
        assert orbit==set(row['members'])

def test_high_precision_quotient_stationarity():
    import mpmath as mp
    q=load()['vacuum']['quotient_control']
    with mp.workdps(55):
        z=tuple(map(mp.mpf,q['coordinates']));fun=lambda *a:P.quotient_value(a,mp)
        gradient=[mp.diff(fun,z,tuple(int(i==j) for j in range(20))) for i in range(20)]
        assert max(abs(v) for v in gradient)<mp.mpf('1e-40')
        assert abs(fun(*z)-mp.mpf(q['value']))<mp.mpf('1e-45')
        assert all(mp.mpf(v)>0 for v in q['Hessian_eigenvalues'])
        assert all(z[i]>0 for i in range(3)) and z[5]>0 and z[9]>0
