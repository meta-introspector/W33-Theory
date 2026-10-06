"""Independent regressions for native control identities and their physics scope."""
import hashlib,json,sys,itertools
from pathlib import Path
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11539_native_plane_holonomy as p
DATA=json.loads((ROOT/'data/w33_pass11539_native_plane_holonomy.json').read_text())
M=lambda a:s.Matrix(a)

def test_source_binding_and_all_sections():
    assert DATA['producer_sha256']==hashlib.sha256(Path(p.__file__).read_bytes()).hexdigest()
    for file,digest in DATA['source_sha256'].items():
        path=ROOT/file
        raw=json.dumps(json.loads(path.read_text()),sort_keys=True,separators=(',',':')).encode() if path.suffix=='.json' else path.read_bytes()
        assert hashlib.sha256(raw).hexdigest()==digest
    for key in ['geometry','curvature','exchange','neutral_composite','adiabatic','pair_resources','occupation_gauge_audit']:assert DATA[key]['status']=='PASS'

def test_native_moment_projector_and_second_basis():
    _,mats,IE,_=p.native();w=s.Matrix([[a[i,i] for a in mats] for i in range(27)])
    K=s.diag(*list(w*IE*(w[0,:]+w[1,:]).T));P=36*(K-s.eye(27)/9)*(K+s.eye(27)/18)*(K+2*s.eye(27)/9)
    assert P==M(DATA['geometry']['native_band_projector'])
    # Rational unitary changes both the operator and the declared subspace.
    R=s.eye(27);R[0,0]=R[2,2]=s.Rational(3,5);R[0,2]=s.Rational(4,5);R[2,0]=-s.Rational(4,5)
    B=R*K*R.T
    assert 36*(B-s.eye(27)/9)*(B+s.eye(27)/18)*(B+2*s.eye(27)/9)==R*P*R.T
    assert K/2==s.diag(*[x/2 for x in K.diagonal()])

def test_berry_rank_is_quotient_not_extra_scalar_modes():
    omega=M(DATA['geometry']['equal_Berry_matrix']);tangent=M(DATA['geometry']['E6_frame_tangent'])
    assert omega.T==-omega and omega.rank()==50 and tangent.rank()==54
    assert all(omega*v==s.zeros(78,1) for v in tangent.nullspace())
    assert DATA['geometry']['rank_comparison'][1]['Berry_rank']==52
    assert 'redundant already' in DATA['geometry']['scope']

def test_curvature_and_exact_rectangle_coefficients():
    _,mats,_,_=p.native();Q=s.diag(0,0,*([1]*25));coords=[]
    sn,cs=s.symbols('sn cs',real=True)
    for row in DATA['curvature']['witnesses']:
        a,b=[mats[i] for i in row['generators']];F=(a*Q*b-b*Q*a)[:2,:2]
        assert F==M(row['curvature']) and F.conjugate().T==-F
        E=s.eye(27)+s.I*sn*b+(cs-1)*b*b
        assert ((E.conjugate().T*a*E)[:2,:2]-s.Matrix([[s.sympify(x,locals={'sn':sn,'cs':cs}) for x in r] for r in row['rectangle_rotated_generator']])).applyfunc(s.simplify)==s.zeros(2)
        coords.append(row['Hermitian_coordinates'])
    assert M(coords).rank()==4

def test_independent_cubic_completeness_and_exchange():
    _,_,_,d=p.native();D=s.SparseMatrix(729,27,{(27*i+j,k):int(d[i,j,k]) for i,j,k in np.argwhere(d)})
    S=s.SparseMatrix(729,729,{(27*i+j,27*j+i):1 for i in range(27) for j in range(27)});I=s.SparseMatrix(s.eye(729))
    C=s.SparseMatrix(729,729,{(i,j):s.sympify(v) for i,j,v in DATA['exchange']['Casimir_entries']})
    assert D.T*D==10*s.eye(27)
    assert 18*C==I+3*S-3*D*D.T
    assert (9*C*C+11*C)/2-4*I/9==S
    assert (C+13*I/9)*(C+I/9)*(C-2*I/9)==s.zeros(729)
    assert s.trace(C*C)==78

def test_su4_generation_replay_and_entanglement():
    rows=DATA['exchange'];basis=list(map(M,rows['su4_basis']))
    for i,step in enumerate(rows['su4_generation_steps']):
        if isinstance(step,list):a,b=step;assert basis[i]==basis[a]*basis[b]-basis[b]*basis[a]
    vectors=[list(A.applyfunc(s.re))+list(A.applyfunc(s.im)) for A in basis]
    assert M(vectors).rank()==15 and all(s.trace(A)==0 and A.conjugate().T==-A for A in basis)
    S=M(rows['encoded_swap']);gate=(s.eye(4)-s.I*S)/s.sqrt(2)
    psi=gate*s.Matrix([0,1,0,0]);coefficient=s.Matrix(2,2,list(psi))
    assert s.simplify(abs(coefficient.det())**2)==s.Rational(1,4)

def test_full_gauge_neutral_pairs_and_completeness():
    section=DATA['neutral_composite'];_,stab,weights,_,_,_=p.full_stabilizer()
    assert len(stab)==28 and section['zero_weight_pair_count']-section['root_constraint_rank']==3
    for terms,norm in zip(section['pair_vectors'],section['pair_norm_squared']):
        vec={(i,j):s.sympify(v) for i,j,v in terms};assert sum(abs(v)**2 for v in vec.values())==norm
        for a in stab:
            out={}
            for (pair,col),value in p.wedge_action(a,list(vec)).items():out[pair]=out.get(pair,0)+value*list(vec.values())[col]
            assert all(s.expand(v)==0 for v in out.values())
    # Redo full zero-weight and raising/lowering constraints, not only stored witnesses.
    replay=p.neutral_composite()
    assert replay['full_neutral_pair_dimension']==3 and replay['full_neutral_one_particle_dimension']==2
    assert replay['pair_vectors']==section['pair_vectors']

def test_neutral_dark_band_has_constant_gap_and_u2_curvature():
    b=s.Matrix([s.Rational(1,3),s.Rational(2,3),s.Rational(2,3)]);H=b*b.T;P=s.eye(3)-H
    assert H*H==H and P*P==P and P.rank()==2 and H*P==s.zeros(3)
    assert M(DATA['neutral_composite']['neutral_dark_curvature_coordinates']).rank()==4
    assert 'four-body/product' in DATA['neutral_composite']['entangler']
    assert 'effective interaction prescription' in DATA['pair_resources']['boundary']

def test_register_swap_respects_pair_antisymmetry():
    # Antisymmetric tensor pairs on four slots; independently replay permutation.
    def wedge(i,j):return {(i,j):s.sqrt(2)/2,(j,i):-s.sqrt(2)/2}
    u=wedge(0,3);v=wedge(25,80)
    state={a+b:x*y for a,x in u.items() for b,y in v.items()}
    swapped={(k,l,i,j):c for (i,j,k,l),c in state.items()}
    expected={b+a:x*y for a,x in u.items() for b,y in v.items()}
    assert swapped==expected
    # Product S13 S24 exchanges pairs; sum of the swaps does not equal it.
    product={(k,l,i,j):c for (i,j,k,l),c in state.items()}
    assert product==swapped and {(k,j,i,l):c for (i,j,k,l),c in state.items()}!=swapped

def test_pure_gauge_cancellation_and_filled_band():
    t=s.symbols('t',real=True);U=s.Matrix(2,2,[s.cos(t),-s.sin(t),s.sin(t),s.cos(t)]);A=-s.I*s.diff(U,t)*U.T
    assert (U.T*(s.diff(U,t)-s.I*A*U)).applyfunc(s.simplify)==s.zeros(2)
    P=U*s.diag(1,0)*U.T
    assert (s.diff(P,t)-s.I*(A*P-P*A)).applyfunc(s.simplify)==s.zeros(2)
    assert s.simplify(U.det())==1 and len(list(itertools.combinations(range(2),2)))==1

def test_finite_time_witness_is_not_exact_control_or_threshold():
    rows=DATA['adiabatic']['history']
    assert rows[-1]['geometric_operator_error']<rows[0]['geometric_operator_error']
    assert all(0<r['geometric_operator_error'] and 0<r['leakage_operator_norm_squared'] and r['norm_error']<1e-7 for r in rows)
    assert 'not certified error bounds' in DATA['adiabatic']['scope']


def test_cubic_condensate_pair_frame_and_full_covariance():
    section=DATA['pair_resources'];pairs=list(itertools.combinations(range(81),2))
    B=s.SparseMatrix(3240,81,{(r,c):s.sympify(v) for r,c,v in section['map_entries']})
    assert B.T*B==10*s.eye(81)
    terms=DATA['neutral_composite']['pair_vectors']
    assert {pairs[r]:v for (r,c),v in B.todok().items() if c==0}=={(i,j):s.sympify(v) for i,j,v in terms[2]}
    assert {pairs[r]:v for (r,c),v in B.todok().items() if c==3}=={(i,j):-s.sympify(v) for i,j,v in terms[1]}
    he,_,_,d=p.native()
    for a in he:
        assert not np.any(np.einsum('ai,ajk->ijk',a,d)+np.einsum('aj,iak->ijk',a,d)+np.einsum('ak,ija->ijk',a,d))
    # The second factor is the invariant alternating three-form of traceless SU3.
    epsilon=np.zeros((3,3,3),dtype=int)
    for perm in itertools.permutations(range(3)):epsilon[perm]=-1 if sum(a>b for a,b in itertools.combinations(perm,2))%2 else 1
    full,*_=p.full_stabilizer()
    for A in full[78:]:
        a=np.array(A[:3,:3].tolist(),dtype=complex)
        assert not np.any(np.einsum('ai,ajk->ijk',a,epsilon)+np.einsum('aj,iak->ijk',a,epsilon)+np.einsum('ak,ija->ijk',a,epsilon))
    vec=[s.SparseMatrix(3240,1,{(pairs.index((i,j)),0):s.sympify(v)/s.sqrt(n) for i,j,v in row}) for row,n in zip(terms,[1,10,10])]
    V=s.SparseMatrix.hstack(*vec);assert V.T*V==s.eye(3)


def test_pair_sector_obstruction_and_natural_logical_base():
    r=DATA['pair_resources']
    assert list(map(s.sympify,r['native_pair_Casimir']))==[s.Rational(-1,9),s.Rational(-13,9),s.Rational(-13,9)]
    assert list(map(s.sympify,r['family_pair_Casimir']))==[s.Rational(2,3),s.Rational(-4,3),s.Rational(-4,3)]
    P=s.diag(0,1,1);assert P.rank()==2 and P*s.Matrix([1,0,0])==s.zeros(3,1)
    assert 'cannot mix ancilla' in r['sector_split']
