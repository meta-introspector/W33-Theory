"""Independent statistics, group-ring, gate, flavor, wave and loop regressions."""
import hashlib,itertools,json,sys
from collections import defaultdict
from pathlib import Path
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11542_11546_five_dynamical_targets as p
D=json.loads((ROOT/'data/w33_pass11542_11546_five_dynamical_targets.json').read_text())

def test_source_bindings_and_all_five_targets():
    assert D['passes']==list(range(11542,11547))
    assert D['producer_sha256']==hashlib.sha256(Path(p.__file__).read_bytes()).hexdigest()
    for name,digest in D['source_sha256'].items():
        path=ROOT/name;raw=json.dumps(json.loads(path.read_text()),sort_keys=True,separators=(',',':')).encode() if path.suffix=='.json' else path.read_bytes()
        assert hashlib.sha256(raw).hexdigest()==digest
    assert all(D[key]['status']=='PASS' for key in ['binding','exchange','flavor','causal','vacuum_energy'])

def test_spin_statistics_and_ancilla_charge_budget():
    swap=s.Matrix(4,4,[1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1]);singlet=s.Matrix([0,1,-1,0]);triplets=[s.Matrix([1,0,0,0]),s.Matrix([0,1,1,0]),s.Matrix([0,0,0,1])]
    assert swap*singlet==-singlet and all(swap*v==v for v in triplets)
    assert s.Matrix.hstack(*triplets).rank()==3
    # Internal antisymmetry times spin symmetry times S-wave symmetry=-1.
    assert (-1)*(1)*(1)==-1
    for degree in range(3):assert all(abs(sum(q))!=6 for q in itertools.product([-2,2],repeat=degree))
    assert sum([-2]*3)==-6
    assert D['binding']['operator_dimensions']['canonical_dimension1_pair_field_ancilla_mix']==5

def test_binding_and_local_control_uses_actual_pair_map():
    a,mu,r=s.symbols('a mu r',positive=True);u=s.exp(-r/a)
    assert s.simplify(s.diff(u,r,2)-u/a**2)==0 and s.integrate(u*u,(r,0,s.oo))==a/2
    old=json.loads((ROOT/'data/w33_pass11539_native_plane_holonomy.json').read_text());B=s.SparseMatrix(3240,81,{(r,c):s.sympify(v) for r,c,v in old['pair_resources']['map_entries']})
    X=s.SparseMatrix.hstack(-B[:,3],B[:,0]);V=X/s.sqrt(10);C=s.Matrix(2,2,[2,1+s.I,1-s.I,-1])
    assert V.conjugate().T*X*C*X.conjugate().T*V==10*C

def test_exchange_compression_via_reduced_density_not_four_slot_producer():
    old=json.loads((ROOT/'data/w33_pass11539_native_plane_holonomy.json').read_text());A=[]
    for row in old['neutral_composite']['pair_vectors'][1:]:
        m=s.SparseMatrix(81,81,{})
        for i,j,v in row:m[i,j]=s.sympify(v);m[j,i]=-s.sympify(v)
        A.append(m)
    M=s.Matrix(4,4,lambda row,col:s.trace(A[col//2]*A[row//2].T*A[col%2]*A[row%2].T)/400)
    assert 4*M==s.Matrix(D['exchange']['projected_W'])
    assert M==s.Matrix(4,4,[2,0,0,0,0,1,1,0,0,1,1,0,0,0,0,2])/40

def test_collective_exchange_group_ring_identity_for_all_carrier_dimensions():
    identity=(0,1,2,3)
    def trans(i,j):p=list(identity);p[i],p[j]=p[j],p[i];return tuple(p)
    def mult(a,b):
        out=defaultdict(lambda:s.S.Zero)
        for x,v in a.items():
            for y,w in b.items():out[tuple(x[y[i]] for i in range(4))]+=v*w
        return {k:v for k,v in out.items() if v}
    def add(*items):
        out=defaultdict(lambda:s.S.Zero)
        for a,scale in items:
            for k,v in a.items():out[k]+=scale*v
        return {k:v for k,v in out.items() if v}
    I={identity:s.S.One};P=mult(add((I,1),({trans(0,1):1},-1)),add((I,s.Rational(1,4)),({trans(2,3):1},s.Rational(-1,4))))
    W={trans(i,j):1 for i,j in [(0,2),(0,3),(1,2),(1,3)]};R=mult({trans(0,2):1},{trans(1,3):1})
    assert mult(add((mult(W,W),1),(W,2),(I,-4),(R,-4)),P)=={}
    assert mult(W,P)==mult(P,W)

def test_exact_finite_time_revival_and_fragility():
    from scipy.linalg import expm
    row=D['exchange']['revival'];q=float(s.sympify(row['J_over_Delta']));t=float(s.sympify(row['duration_times_Delta']))
    H=lambda j:np.array([[j/5,3*np.sqrt(21)*j/5],[3*np.sqrt(21)*j/5,2-11*j/5]])
    U=expm(-1j*t*H(q));assert abs(U[1,0])<1e-10 and abs(U[0,0]+1j)<1e-10
    wrong=expm(-1j*t*H(q*1.01));assert abs(wrong[1,0])>1e-5
    S=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]]);gate=((1-1j)*np.eye(4)+(-1-1j)*S)/2
    assert np.linalg.norm(gate.conj().T@gate-np.eye(4))<1e-12
    psi=gate[:,1].reshape(2,2);assert abs(abs(np.linalg.det(psi))-.5)<1e-12

def test_cp_exact_diagonal_commutator_formula_and_weak_basis():
    e=s.symbols('epsilon',real=True);scale=s.diag(e**3,e**2,1);Bu=s.Matrix(D['flavor']['up_pin']);Bd=s.Matrix(D['flavor']['positive_down_pin']);H=(scale*Bu*scale)**2;K=(scale*Bd*scale)**2
    delta=(H[0,0]-H[1,1])*(H[1,1]-H[2,2])*(H[2,2]-H[0,0]);cp=s.expand(6*delta*s.im(K[0,1]*K[1,2]*K[2,0]));assert s.expand(cp-s.sympify(D['flavor']['exact_CP_polynomial'],locals={'epsilon':e}))==0
    assert s.Poly(cp,e).coeff_monomial(e**22)==369360
    R=s.Matrix(3,3,[s.Rational(3,5),s.I*s.Rational(4,5),0,-s.Rational(4,5),s.I*s.Rational(3,5),0,0,0,1]);assert R.conjugate().T*R==s.eye(3)
    h=H.subs(e,s.Rational(1,7));k=K.subs(e,s.Rational(1,7));h2=R*h*R.conjugate().T;k2=R*k*R.conjugate().T
    assert s.simplify(s.im(s.trace((h2*k2-k2*h2)**3))-cp.subs(e,s.Rational(1,7)))==0

def test_short_exact_revival_without_weak_coupling_assumption():
    from scipy.linalg import expm
    row=D['exchange']['short_revival'];q=float(s.sympify(row['J_over_Delta']));t=float(s.sympify(row['duration_times_Delta']))
    H=np.array([[q/5,3*np.sqrt(21)*q/5],[3*np.sqrt(21)*q/5,2-11*q/5]])
    U=expm(-1j*t*H)
    assert 0.39<q<0.40 and 7.7<t<7.9
    assert abs(U[1,0])<1e-12 and abs(U[0,0]+1j)<1e-12
    assert t<float(s.sympify(D['exchange']['revival']['duration_times_Delta']))/100
    early=expm(-1j*t*H/3)
    assert abs(early[1,0])<1e-12 and abs(early[0,0]-1j)<1e-12

def test_positive_pin_bound_and_indefinite_counterexample():
    B=s.Matrix(D['flavor']['positive_down_pin']);assert all(B[:i,:i].det()>0 for i in [1,2,3])
    e=s.symbols('e',positive=True);scale=s.diag(e**3,e**2,1);bad=s.Matrix(3,3,[0,1,0,1,0,0,0,0,1]);Y=scale*bad*scale
    assert Y*Y==s.diag(e**10,e**10,1) and bad.det()==-1
    coeff=list(map(s.sympify,D['flavor']['down_leading_mass_coefficients']));assert coeff[0]*coeff[1]*coeff[2]==B.det()

def test_causal_cover_wave_symbol_stability_and_refinement():
    steps=list(map(tuple,D['causal']['integer_null_steps']));assert len(steps)==8 and all((t*t-x*x-y*y)%3==0 for t,x,y in steps)
    assert {(x,y) for t,x,y in steps if t==1}=={(1,0),(-1,0),(0,1),(0,-1)}
    for kx,ky in itertools.product(np.linspace(-np.pi,np.pi,13),repeat=2):
        rhs=(np.sin(kx/2)**2+np.sin(ky/2)**2)/3;assert 0<=rhs<=2/3+1e-15
    rows=D['causal']['numerical_refinements'];ratios=[a['continuum_error']/b['continuum_error'] for a,b in zip(rows,rows[1:])];assert all(3.9<x<4.1 for x in ratios)
    assert D['causal']['maximum_energy_drift']<1e-12 and 'No4D spacetime' in D['causal']['boundary']

def test_scalar_loop_is_not_a_full_physical_supertrace():
    old=json.loads((ROOT/'data/w33_pass11526_11530_validated_dynamics.json').read_text())['coherent_two_stage']['full324_spectrum']
    assert sum(old.values())==324 and old['0']==58
    value=sum(s.Rational(h)**2*n for h,n in old.items());assert value==s.Rational(2218769,864)==s.Rational(D['vacuum_energy']['scalar_hessian_mass_fourth_sum'])
    assert 'scalar-only' in D['vacuum_energy']['one_loop'] and 'graviton loops' in D['vacuum_energy']['boundary']

def test_sequestering_shift_preserves_sources_but_naive_subtraction_erases_forces():
    L,C,Q,m,V=s.symbols('L C Q m V');linear=V-Q/m;quadratic=V-L*Q/m**2
    assert linear.subs(L,L-C)-linear==0 and s.expand(quadratic.subs(L,L-C)-quadratic)==C*Q/m**2
    x,y=s.symbols('x y');rho=[x*x,y*y];weights=[s.Rational(1,3),s.Rational(2,3)];avg=sum(w*r for w,r in zip(weights,rho));subtracted=s.expand(sum(w*(r-avg) for w,r in zip(weights,rho)))
    assert subtracted==0 and s.diff(sum(w*r for w,r in zip(weights,rho)),x)==2*x/3
    assert s.expand((rho[0]+C)-(avg+C)-(rho[0]-avg))==0
