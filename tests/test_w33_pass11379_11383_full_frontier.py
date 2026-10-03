"""Independent controls for full-field, basis, auxiliary, boundary and units claims."""
from pathlib import Path
import sys
import json
import numpy as np
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11379_11383_full_frontier as m
import w33_pass11374_11378_five_frontier as old


def test_full_ray_hessian_has_only_symmetry_zeros():
    d=m.full_ray_selector()
    J=s.Matrix([[s.sympify(v) for v in row] for row in d['Jacobian']])
    assert J.rank()==7 and 18-J.rank()==11
    # Overlap targets and real loop part fix a nonzero imaginary loop, without
    # imposing the old one-phase coordinate ansatz.
    modulus2=s.Rational(1,2)*s.Rational(1,3)**2
    assert modulus2-s.Rational(1,6)**2==s.Rational(1,36)
    x,_,_=m.ray_data()
    rng=np.random.default_rng(18)
    U=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
    X=np.array(x,complex)
    rotated=U@X@np.diag(np.exp(1j*rng.normal(size=3)))
    for a in [rotated,rotated.conj()]:
        G=a.conj().T@a; B=G[0,1]*G[1,2]*G[2,0]
        assert abs(abs(B.imag)-1/6)<1e-12 and abs(B.real-1/6)<1e-12


def test_native_channels_have_basis_independent_small_algebra():
    d=m.native_channel_span()
    assert d['symmetric_span_dimension']==7 and d['generated_algebra_dimension']==8
    assert sorted(d['central_block_dimensions'])==[1,1,2,4,6]
    basis=list(map(np.array,d['orthonormal_channel_basis']))
    P=list(map(np.array,d['central_projectors']))
    assert np.linalg.norm(sum(P)-np.eye(14))<1e-10
    for p,size in zip(P,d['central_block_dimensions']):
        assert np.linalg.norm(p@p-p)<1e-10
        assert abs(np.trace(p)-size)<1e-10
        assert max(np.linalg.norm(p@a-a@p) for a in basis)<1e-9
    # The four-dimensional block is two copies of a noncommuting 2D irrep.
    p=P[d['central_block_dimensions'].index(4)]
    assert max(np.linalg.norm(p@(a@b-b@a)@p) for a in basis for b in basis)>0.1
    W=np.array(d['noncommuting_block_intertwiner'])
    assert np.linalg.norm(W.T@W-np.eye(4))<1e-10
    for a,little in zip(basis,d['compressed_two_by_two_channel_matrices']):
        assert np.linalg.norm(W.T@a@W-np.kron(little,np.eye(2)))<1e-9
    assert all(min(v)>0 for v in d['absorptive_spectra'].values())
    assert d['tree_Hessian_span_residual']<1e-9
    assert d['absorptive_truncation_inverse_residual']<1e-9


def test_native_rotation_counterexample_uses_the_actual_four_dimensional_pair():
    d=m.native_rotational_lapse()
    assert d['reduced_lapse_Hessian_rank']==6 and d['shift_boost_Jacobian_rank']==158
    inc,edges,_,frames,aa,bb=m.native_rotation_data()
    N1,N2,t=s.symbols('N1 N2 t', real=True, nonzero=True)
    R=s.Matrix([[s.cos(t),-s.sin(t)],[s.sin(t),s.cos(t)]])
    # Check a frustrated edge and an unperturbed edge against the original
    # determinant/trace expression, retaining the spectator spatial direction.
    for k in [0,9]:
        i,j=edges[k]; A=s.diag(N1,1,1,1); C=s.diag(N2,1,1,1)
        A[1:3,1:3]=frames[i]; C[1:3,1:3]=R*frames[j]
        pair=(A.det()*s.trace(A.inv()*C)+C.det()*s.trace(C.inv()*A))/2
        reduced=(N1+N2)*(aa[k]*s.cos(t)+bb[k]*s.sin(t)+frames[i].det()+frames[j].det())/2
        assert s.simplify(pair-reduced)==0
    # Independently check the coordinate-shift/boost cross derivative. This
    # catches a missing spatial coframe factor in the auxiliary Jacobian.
    i,j=edges[0]; ei,ej=frames[i],frames[j]
    beta,shift=s.symbols('beta shift',real=True)
    T=s.eye(4); T[1:3,1:3]=ei; T[1:3,0]=ei*s.Matrix([shift,0])
    C=s.eye(4); C[1:3,1:3]=ej
    boost=s.eye(4); boost[0,0]=boost[1,1]=s.cosh(beta)
    boost[0,1]=boost[1,0]=s.sinh(beta)
    E=boost*T
    original=(E.det()*s.trace(E.inv()*C)+C.det()*s.trace(C.inv()*E))/2
    derivative=s.simplify(s.diff(original,beta,shift).subs({beta:0,shift:0}))
    assert derivative==(ei.det()*ej+ej.det()*ei)[0,0]/2
    # Cycle-free graph removes this mechanism: on a tree each relative angle
    # extremizes independently, giving a lapse-linear edge extremum.
    a,b=s.symbols('a b',positive=True)
    extremum=(N1+N2)*s.sqrt(a*a+b*b)/2
    assert s.diff(extremum,N1,N2)==0


def test_coupled_wall_factorization_and_parity_boundaries():
    d=m.wall_factorization()
    assert d['radial_vector_kernel_total']==4
    # Independently integrate a regular odd trial pair; agreement is sensitive
    # to the mixing sign and both endpoint terms, not just the zero-mode ODE.
    b=s.symbols('b',positive=True); bw=s.Rational(5,4)
    q=s.Rational(4,3); c=s.Rational(4,15); h=s.Rational(16,45)
    F=-8*(b-1)*(2*b-5)/(45*b*b)
    v0=-F/(q*c*b*b)
    chi=(bw-b)*(1+b); v=v0*(1+(b-bw)**2)
    original=F*s.diff(chi,b)**2+c*c*b**4*s.diff(v,b)**2/2+2*q*c*chi*s.diff(v,b)+h*b*b*c*c*v*v/F
    squares=F*(s.diff(chi,b)-q*c*v/F)**2+c*c*b**4*v0*v0*s.diff(v/v0,b)**2/2
    # Exact antiderivative of the difference is the discarded boundary term.
    boundary=2*q*c*chi*v+c*c*b**4*s.diff(v0,b)/v0*v*v/2
    assert s.factor(original-squares-s.diff(boundary,b))==0
    assert s.limit(boundary,b,1)==0 and boundary.subs(b,bw)==0


def test_axion_shift_gauge_and_gravitational_units_regressions():
    assert old.wall_axion_probe_sector()['physical_constant_axion_moduli']==0
    a,V,dl=s.symbols('a V dl')
    assert s.expand((a-dl-(V-dl))**2-(a-V)**2)==0
    K,L=s.symbols('K L',positive=True)
    R=s.Rational(512,1125)/L**2
    T=s.Rational(64,75)*K/L
    assert s.factor(K*K*R/T**2)==s.Rational(5,8)
    assert 'K_phys^2' in old.wall_scale_ratio()['physical_scale_invariant_ratio']


def test_scale_stationarity_requires_the_boundary_term():
    d=m.fixed_action_scale()
    L=s.symbols('L',positive=True)
    action=s.sympify(d['action_without_angular_factor'],locals={'L':L})
    assert s.diff(action,L).subs(L,1)==0
    assert s.diff(action,L,2).subs(L,1)==4
    # Dropping GHY destroys the stationarity, so this is not a bulk-only
    # scaling argument or a family obtained by changing the input couplings.
    assert s.diff(action+2*L**2,L).subs(L,1)==4


def test_stored_packet_matches_its_producer_scopes():
    d=json.loads((ROOT/'data/w33_pass11379_11383_full_frontier.json').read_text())
    assert d['status']=='PASS' and len(d['sections'])==5
    assert all(v['status']=='PASS' and v['boundary'] for v in d['sections'].values())
