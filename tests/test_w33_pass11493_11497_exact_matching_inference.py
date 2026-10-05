"""Independent checks for identities, actual objects, and excluded over-reads."""
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11493_11497_exact_matching_inference as M


def certificate():
    return json.loads(M.OUT.read_text())


def test_source_binding_and_scopes():
    c=certificate()
    assert c['passes']==list(range(11493,11498)) and c['reservation']=='2cda9139b'
    assert len(c['source_sha256'])==3 and 'data/w33_pass11293_yukawa_uv_budget.json' in c['source_sha256']
    for path,digest in c['source_sha256'].items():
        assert hashlib.sha256(json.dumps(M.read(path),sort_keys=True,separators=(',',':')).encode()).hexdigest()==digest
    for key in ['stationary','metric','portal','ward','joint']:
        assert c[key]['status']=='PASS' and len(c[key]['scope'])>180


def test_hollowisation_and_final_full_gradient_replay():
    import w33_pass11438_finite_native_model as F
    r=certificate()['stationary'];old=M.read(M.OLD)['vacuum']
    U=M.dec(r['hollowising_unitary']);E=M.dec(old['isotope']['native_tripotents'])
    assert np.linalg.norm(U.conj().T@U-np.eye(3))<1e-12 and abs(np.linalg.det(U)-1)<1e-12
    for key in ['A_at_vacuum','B_at_vacuum']:
        A=M.dec(old['polydisc'][key]);G=U@A@A.conj().T@U.conj().T
        assert np.linalg.norm(np.diag(G)-np.trace(G)/3)<1e-9
    z=np.array(r['final_reduced_coordinates'])
    A=(z[:9]+1j*z[9:18]).reshape(3,3);B=(z[18:27]+1j*z[27:]).reshape(3,3)
    with F.native_context():
        value,g=F.fun(F.pack((E@A).ravel(),(E@B).ravel()))
    assert abs(value-r['history'][-1]['value'])<1e-10
    assert abs(np.linalg.norm(g)-r['history'][-1]['full_gradient_norm'])<1e-10
    assert all(b['value']<=a['value']+1e-12 for a,b in zip(r['history'],r['history'][1:]))


def test_schur_complement_is_curvature_after_relaxation():
    # Exact mixed two-light/two-heavy polynomial, not producer's scalar control.
    x=sp.Matrix(sp.symbols('x0:2')); y=sp.Matrix(sp.symbols('y0:2'))
    A=sp.Matrix([[5,1],[1,4]]); B=sp.Matrix([[1,2],[3,-1]]); C=sp.Matrix([[4,1],[1,3]])
    V=(x.T*A*x)[0]/2+(x.T*B*y)[0]+(y.T*C*y)[0]/2
    solved=-C.inv()*B.T*x
    reduced=sp.expand(V.subs(dict(zip(y,solved)),simultaneous=True))
    assert sp.hessian(reduced,list(x))==A-B*C.inv()*B.T


def determinant_mod(matrix,prime):
    a=[[int(v)%prime for v in row] for row in matrix]; result=1
    for j in range(len(a)):
        pivot=next((i for i in range(j,len(a)) if a[i][j]),None)
        if pivot is None:return 0
        if pivot!=j:a[j],a[pivot]=a[pivot],a[j];result=-result
        value=a[j][j];result=result*value%prime; inverse=pow(value,-1,prime)
        for i in range(j+1,len(a)):
            factor=a[i][j]*inverse%prime
            for k in range(j+1,len(a)):a[i][k]=(a[i][k]-factor*a[j][k])%prime
    return result%prime


def test_native_cofactor_independent_modular_reconstruction():
    r=certificate()['metric'];det=Fraction(r['exact_reduced_Laplacian_determinant'])
    for prime in [1000003,1000033,1000037]:
        L=[[0]*79 for _ in range(79)]
        for (a,b),length in zip(r['edges'],r['squared_integer_lengths']):
            w=6400*pow(length,-1,prime)%prime
            if a<79:L[a][a]+=w
            if b<79:L[b][b]+=w
            if a<79 and b<79:L[a][b]-=w;L[b][a]-=w
        assert determinant_mod(L,prime)==det.numerator*pow(det.denominator,-1,prime)%prime


def test_metric_scale_and_second_coordinate_basis():
    r=certificate()['metric'];edges=np.array(r['edges']);h=np.array(r['displacements_times80'])/80
    B=np.zeros((160,80));B[np.arange(160),edges[:,0]]=1;B[np.arange(160),edges[:,1]]=-1
    def W(G,h,cell=1):
        w=cell*np.sqrt(np.linalg.det(G))/np.einsum('ei,ij,ej->e',h,G,h)
        return np.linalg.slogdet(B[:,:-1].T@(w[:,None]*B[:,:-1]))[1]/2
    base=W(np.eye(3),h)
    for scale in [.7,1.4,2.]:assert abs(W(scale*np.eye(3),h)-base-79/4*np.log(scale))<1e-10
    J=np.array([[1,.2,0],[.1,1.2,.3],[0,.1,.9]]);inv=np.linalg.inv(J)
    assert abs(W(inv.T@inv,h@J.T,np.linalg.det(J))-base)<1e-10
    assert abs(sum(r['edge_tree_probabilities'])-79)<1e-9


def test_portal_exact_matching_stability_and_SM_obstruction():
    import w33_pass11384_11388_native_dynamics as N
    d=np.array(N.native_tensors()[1][1]);r=certificate()['portal']
    flat=d.reshape(27,729);c=int(r['exact_cubic_contraction_constant'])
    assert np.array_equal(flat@flat.T,c*np.eye(27)) and not np.any(d[:2,:2,:])
    assert r['SM_pair_source_zero'] and r['operator_dimension']==5
    rng=np.random.default_rng(48);phi=rng.normal(size=(27,3))+1j*rng.normal(size=(27,3))
    source=np.einsum('abc,ai,bj->cij',d,phi.conj(),phi.conj())
    assert np.linalg.norm(source)**2<=c*np.linalg.norm(phi)**4+1e-8
    H=rng.normal(size=(27,3,3))+1j*rng.normal(size=(27,3,3));kappa=.3+.2j;mass2=4.
    V=mass2*np.linalg.norm(H)**2-2*np.vdot(H,kappa*source).real
    matched=mass2*np.linalg.norm(H-kappa*source/mass2)**2-abs(kappa)**2*np.linalg.norm(source)**2/mass2
    assert abs(V-matched)<1e-10


def test_induced_yukawa_actual_carrier_by_blocks():
    import w33_pass11384_11388_native_dynamics as N
    d=np.array(N.native_tensors()[1][1]);r=certificate()['portal']
    for source,matrix in zip(r['source_matrices'],r['induced_Yukawas']):
        source=M.dec(source);Y=M.dec(matrix)
        blocks=np.block([[sum(d[:,:,k]*source[k,i,j] for k in range(27)) for j in range(3)] for i in range(3)])
        order=np.arange(81).reshape(3,27).T.ravel()
        assert np.linalg.norm(Y-blocks[np.ix_(order,order)])<1e-12
        assert np.linalg.norm(Y-Y.T)<1e-12


def test_support_obstruction_independent_cut_projection():
    # Connected path with a three-site active ball and an unavoidable outside tail.
    target=np.array([1.,2.,-1.,3.,-5.]);D=np.zeros((5,2));D[0,0]=1;D[1,0]=-1;D[1,1]=1;D[2,1]=-1
    flow=np.linalg.lstsq(D,-target,rcond=None)[0];residual=D@flow+target
    assert np.allclose(residual[:3],target[:3].mean())
    assert np.dot(residual,residual)==np.dot(target[3:],target[3:])+sum(target[:3])**2/3
    r=certificate()['ward'];assert r['charge_moments']['3']==0 and r['charge_moments']['5']==4320
    assert r['rows'][0]['exact_projection_formula_numeric']>1e-9


def test_joint_rational_MAP_and_noisy_flags_are_CPTP():
    r=certificate()['joint'];good=Fraction(r['logical_weight']);bad=Fraction(r['accepted_leakage_weight']);erasure=Fraction(r['rejected_erasure_weight'])
    assert good+bad+erasure==1 and bad>0
    assert Fraction(r['exact_joint_success'])>Fraction(r['exact_separate_success'])>=Fraction(r['exact_raw_success'])
    bell=np.zeros(24);bell[[0,7,14,21]]=1/2
    flags=np.zeros((6,6));flags[4,4]=float(bad);flags[5,5]=float(erasure)
    J=float(good)*np.outer(bell,bell)+np.kron(np.eye(4)/4,flags)
    assert min(np.linalg.eigvalsh(J))>-1e-12
    assert np.linalg.norm(np.einsum('iaja->ij',J.reshape(4,6,4,6))-np.eye(4)/4)<1e-12
    assert abs(np.vdot(bell,J@bell)-float(good))<1e-12
    # Distinct supplied rates exercise the exact Bayesian code, not stored floats.
    alt=M.joint_inference(Fraction(1,10),Fraction(1,20),Fraction(1,100))
    assert Fraction(alt['exact_joint_success'])>=Fraction(alt['exact_separate_success'])


def test_exact_shear_law_and_instability_on_actual_cover():
    r=certificate()['metric'];hi=np.array(r['displacements_times80']);edges=np.array(r['edges'])
    assert np.all(abs(hi)==abs(hi[:,0,None]))
    B=np.zeros((160,79))
    for k,(a,b) in enumerate(edges):
        if a<79:B[k,a]=1
        if b<79:B[k,b]=-1
    def action(diag):
        w=np.sqrt(np.prod(diag))*6400/((hi**2)@diag)
        return np.linalg.slogdet(B.T@(w[:,None]*B))[1]/2
    base=action(np.ones(3))
    for diagonal in [[.8,1.2,1.7],[np.exp(.4),np.exp(-.4),1.]]:
        change=action(np.array(diagonal))-base
        expected=79/2*np.log(3*np.sqrt(np.prod(diagonal))/sum(diagonal))
        assert abs(change-expected)<1e-10
    t=sp.symbols('t',real=True)
    law=-sp.Rational(79,2)*sp.log((sp.exp(t)+sp.exp(-t)+1)/3)
    assert sp.diff(law,t,2).subs(t,0)==-sp.Rational(79,3)


def test_joint_memory_odds_explain_a_specific_decoder_change():
    r=certificate()['joint']['two_report_witness']
    assert r['reports']==[3,3] and r['separate_guess']==0 and r['joint_guess']==3
    assert Fraction(r['identity_vs_fault_single_odds'])==Fraction(2079,361)>1
    assert Fraction(r['identity_vs_fault_joint_odds'])==Fraction(2079,130321)<1


def test_CP_capacity_witness_with_independent_exact_sympy_arithmetic():
    from sympy.polys.matrices import DomainMatrix
    r=certificate()['exact_CP_control'];d=M.dec(M.read(M.SOURCE)['yukawa']['native_cubic'])
    tri=r['canonical_frame'];E=np.zeros((27,3),int)
    for i,k in enumerate(tri[:3]):E[k,i]=tri[3] if i==2 else 1
    grams=[]
    for key in ['A','B']:
        phi=E@M.dec(r[key]);J=np.einsum('abc,ai,bj->cij',d.conj(),phi.conj(),phi.conj())
        Y=np.einsum('abc,cij->aibj',d,J).reshape(81,81)
        assert np.array_equal(Y.real,np.rint(Y.real)) and np.array_equal(Y.imag,np.rint(Y.imag))
        exact=sp.Matrix([[sp.Integer(int(z.real))+sp.I*sp.Integer(int(z.imag)) for z in row] for row in Y])
        grams.append(DomainMatrix.from_Matrix(exact)*DomainMatrix.from_Matrix(exact.conjugate().T))
    C=grams[0]*grams[1]-grams[1]*grams[0]
    CM=C.to_Matrix()
    assert CM.conjugate().T==-CM
    trace=sp.trace((C*C*C).to_Matrix())
    assert trace==sp.I*sp.Integer(r['exact_trace_commutator_cube_imag']) and trace!=0
    assert sp.conjugate(trace)==-trace


def test_SM_neutral_mixed_source_matching_and_full_CP_factorization():
    from sympy.polys.matrices import DomainMatrix
    r=certificate()['SM_mixed_portal'];d=M.dec(M.read(M.SOURCE)['yukawa']['native_cubic'])
    families=[sp.Matrix([[sp.sympify(x) for x in row] for row in f]) for f in r['family_sextets']]
    D=sp.Matrix(np.rint(d[:,:,0].real).astype(int).tolist())
    assert D.rank()==10 and all(f.rank()==3 and f==f.T for f in families)
    grams=[]
    for f in families:
        source=np.zeros((27,3,3),complex);source[0]=np.array(f,complex)
        native=np.einsum('abc,cij->aibj',d,source).reshape(81,81)
        exact=sp.kronecker_product(D,f)
        assert np.linalg.norm(native-np.array(exact,complex))<1e-12
        dm=DomainMatrix.from_Matrix(exact);dag=DomainMatrix.from_Matrix(exact.conjugate().T)
        grams.append(dm*dag)
    C=grams[0]*grams[1]-grams[1]*grams[0]
    trace=sp.trace((C*C*C).to_Matrix())
    assert trace==-216000*sp.I==sp.sympify(r['exact_full_cubic_invariant'])


def test_quadratic_source_selects_hollowisation_but_preserves_fibers():
    r=certificate()['source_norm_fibers'];raw=M.read(M.SOURCE);old=M.read(M.OLD)['vacuum']
    d=M.dec(raw['yukawa']['native_cubic']);tri=old['polydisc']['canonical_frame'];E=np.zeros((27,3),int)
    for i,k in enumerate(tri[:3]):E[k,i]=tri[3] if i==2 else 1
    # Direct full27-tensor contraction at two independent Gaussian-integer matrices.
    for A in [np.array([[1,1j,0],[0,2,1],[1,0,3]]),np.array([[2,1,1j],[1j,1,0],[1,2,2]])]:
        phi=E@A;J=np.einsum('abc,ai,bj->cij',d.conj(),phi.conj(),phi.conj());G=A@A.conj().T
        left=sum(int(z.real)**2+int(z.imag)**2 for z in J.ravel())
        right=np.trace(G).real**2+np.trace(G@G).real-2*sum(np.diag(G).real**2)
        assert left==right
        # Fourier after diagonalisation realizes the Cauchy lower bound for one Gram.
        _,V=np.linalg.eigh(G);F=np.exp(2j*np.pi*np.outer(np.arange(3),np.arange(3))/3)/np.sqrt(3)
        U=F@V.conj().T;hollow=U@G@U.conj().T
        assert abs(sum(np.diag(hollow).real**2)-np.trace(G).real**2/3)<1e-10
    assert r['symbolic_coordinates']==18
    assert max(abs(z['source_norm_shift']) for z in r['controls'])<1e-10
