"""Independent witnesses, Ward identities, nonlinear frame and fluctuation controls."""
from pathlib import Path
import json,sys
import numpy as np
import sympy as s
from scipy.linalg import expm
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
def cert(n):return json.loads(next((ROOT/'data').glob(f'w33_pass{n}_*.json')).read_text())

def test_all_larger_flavor_witnesses_replay_full_complex_matrix_betas():
 import w33_pass11335_larger_complex_flavor_search as F
 r=cert(11335);assert len(r['rows'])==28 and r['both_barriers_nonpositive_count']==2
 for row in r['rows']:
  A,B=F.unpack(np.array(row['coordinates']),*F.bases(row['n']));assert max(np.linalg.norm(x) for x in F.beta(A,B,73/3))<1e-6
  assert abs(F.barrier(A)-row['S_portal_barrier'])<1e-10

def test_relaxed_flavor_norm_inequality_in_independent_complex_samples():
 from w33_pass11335_larger_complex_flavor_search import bases,unpack,beta
 rng=np.random.default_rng(113351)
 for n in [3,5,8]:
  A,B=unpack(rng.normal(size=2*n*n),*bases(n));H=A.conj().T@A;Q=B.conj().T@B;fa,_=beta(A,B);lhs=np.trace(A.conj().T@fa).real;rhs=np.trace(H).real**2+6.2*np.trace(H@H).real+3.5*np.trace(H@Q).real;assert lhs>=rhs-1e-8

def test_exact_isotropic_ray_has_nilpotent_cross_and_two_nonpositive_bounds():
 from w33_pass11340_isotropic_flavor_portal import flavor,beta
 x,w,t,A,B=flavor();assert np.linalg.norm(B@A.conj().T@B)<1e-12
 assert max(np.linalg.norm(v) for v in beta(A,B,73/3))<1e-11
 r=cert(11340);assert s.Rational(r['exact_S_lower_bound'])<0 and s.Rational(r['exact_K_lower_bound'])<0
 # A CP conjugate is flavor-equivalent; complex isotropy is not asserted as CP violation.
 R=np.diag([1,1,-1,1]);assert np.max(abs(R@B@R-B.conj()))<1e-12 and np.max(abs(R@A@R-A.conj()))<1e-12

def test_actual_isotropic_box_replays_on_a_new_so_background():
 import w33_pass11313_symmetric_adjoint_portal_closure as P
 from w33_pass11340_isotropic_flavor_portal import flavor
 _,_,_,A,B=flavor();S,K=P.sample(np.random.default_rng(113401));M=np.kron(A,S)+np.kron(B,K);W=M.conj().T@M;r=cert(11340);assert abs(P.invariants(S,K)@r['fermion_box_coefficients']-np.trace(W@W).real)<1e-8
 assert r['root_count']==len(r['quartic_roots']) # No completeness inference from the finite search.

def test_many_healthy_gaussian_mediators_still_only_make_spectral_quartics():
 from w33_pass11336_gaussian_mediator_alignment import eliminate
 rng=np.random.default_rng(113361);n=4;Z=rng.normal(size=(5,5));M=Z@Z.T+np.eye(5);a=rng.normal(size=5);b=rng.normal(size=5);A=np.diag([1.,2.,3.,5.]);T=rng.normal(size=(n,n));T=(T+T.T)/2;U=expm(1j*T);B=U@np.diag([1.,3.,4.,7.])@U.conj().T;X,v=eliminate(A,B,M,a,b);A0=A-A.trace()*np.eye(n)/n;B0=B-B.trace()*np.eye(n)/n
 assert np.linalg.norm(np.einsum('ab,bij->aij',M,X)+a[:,None,None]*A0+b[:,None,None]*B0)<1e-10
 chi=a@np.linalg.solve(M,b);h=1e-5;U=expm(1j*h*T);V=U.conj().T;fd=(eliminate(A,U@B@U.conj().T,M,a,b)[1]-eliminate(A,V@B@V.conj().T,M,a,b)[1])/(2*h);assert abs(fd-(-1j*chi*np.trace((B@A-A@B)@T)).real)<1e-6

def test_native_boost_currents_and_lorentz_equations_under_changed_lapses():
 from w33_pass11338_native_boost_constraint_slice import solve_boosts,cycle_basis
 B,C=cycle_basis();rng=np.random.default_rng(113381);p=rng.normal(size=80)*.3;p-=p.mean();eta,j=solve_boosts(B,p);cos=np.cosh(B.T@eta);assert max(abs(C.T@np.arcsinh(j)))<1e-10
 for N in [np.ones(80),rng.uniform(.7,1.3,80)]:
  L=abs(B).T@N;shift=-np.linalg.pinv((B*cos)@B.T)@(B@(L*j));assert max(abs(B@(L*j+cos*(B.T@shift))))<1e-9

def test_pair_vielbein_trace_formula_in_actual_random_coframes():
 rng=np.random.default_rng(113382)
 for _ in range(12):
  ni,nj=rng.uniform(.5,1.5,2);si,sj=rng.normal(size=2)*.2;d=rng.normal()*.3;ei=np.eye(4);ei[0,0]=ni;ei[1,0]=si;ej=np.eye(4);ej[0,0]=nj;ej[1,0]=sj;L=np.eye(4);L[:2,:2]=[[np.cosh(d),np.sinh(d)],[np.sinh(d),np.cosh(d)]];ej=L@ej
  actual=(ni*np.trace(np.linalg.solve(ei,ej))+nj*np.trace(np.linalg.solve(ej,ei)))/2;expected=(ni+nj)*np.cosh(d)+(sj-si)*np.sinh(d)+ni+nj;assert abs(actual-expected)<1e-12

def test_hard_ward_columns_are_symmetric_compatible_and_annihilate_first_tadpole():
 r=cert(11337);T=np.array(r['Goldstone_tangents']);K=np.array(r['hard_Ward_columns']);W=np.array(r['symmetric_minimal_completion']);g=np.array(r['hard_gradient']);assert np.linalg.matrix_rank(T)==8
 assert max(abs(T.T@g))<1e-10 and max(abs(W@T-K).flat)<1e-8 and max(abs(W-W.T).flat)<1e-12
 assert r['hard_family_orbit_error']<1e-8

def test_holomorphic_slice_family_tangent_for_another_generator():
 from w33_pass11337_quotient_hard_ward_columns import Chart,orbit_coordinates,H
 C=Chart();F=H.Q.Q.L.generators()[82];h=.0004;xp=orbit_coordinates(C,F,h);xm=orbit_coordinates(C,F,-h);v=C.N.conj().T@(1j*F@C.p);t=np.sqrt(2)*np.r_[v.real,v.imag];assert np.linalg.norm((xp-xm)/(2*h)-t)<1e-7

def test_global_connection_constraint_keeps_constant_wilson_lines_flat():
 z,w=leggauss(400);b=1+(z+1)/8;w/=8;chi=np.sin(3*b)+.2*b;mu=np.sum(w*chi/b**4)/np.sum(w/b**4);hh=-2*(chi-mu)/b**4
 assert abs(np.sum(w*hh))<1e-12
 orig=np.sum(w*(b**4*hh**2/4+chi*hh));reduced=-np.sum(w*(chi-mu)**2/b**4);assert abs(orig-reduced)<1e-12
 assert np.max(abs(-2*(np.ones_like(b)-1)/b**4))==0

def test_coupled_odd_zero_satisfies_bulk_and_both_wall_conditions():
 b=s.symbols('b');F=s.Rational(8,15)/b-b*b/30-1/(2*b*b);u=1/b-s.Rational(4,5);v=s.Rational(5,2)/b**4-s.Rational(8,3)/b**3+s.Rational(1,6)
 assert s.simplify(-s.diff(F*s.diff(u,b),b)-2*u/b**4)==0 and u.subs(b,s.Rational(5,4))==0
 assert v.subs(b,1)==0 and s.diff(v,b).subs(b,s.Rational(5,4))==0
 assert s.simplify(s.diff(v,b)/5+2*u/b**4)==0 and cert(11339)['combined_shape_vector_zero_count']==6

def test_vector_operator_convergence_and_m1_supersolution_positive():
 from w33_pass11339_coupled_maxwell_metric_modes import spectrum
 for m,parity in [(0,'even'),(0,'odd'),(1,'even')]:
  a=spectrum(m,parity,16)[:4];b=spectrum(m,parity,24)[:4];assert min(b)>-1e-8 and max(abs(a-b))<1e-5
 t=s.symbols('t');b=1+t/4
 for p in [b**4-15*b*b-10*b-60,b**3+b*b+b-15]:
  # Upper bound on [0,1]: constant plus sum of all positive coefficients.
  coeff=s.Poly(s.expand(p),t);c=coeff.coeff_monomial(1);assert c+sum(max(0,coeff.coeff_monomial(t**i)) for i in range(1,5))<0

def test_backreacted_winding_wall_replays_at_another_winding_strength():
 from w33_pass11341_winding_lifted_wall import exact_data
 for h in [s.Rational(1,100),s.Rational(1,200)]:
  b,F,d=exact_data(h);q=s.sympify(d['q']);rho=s.sympify(d['rho']);c=s.sympify(d['c']);assert s.simplify(q*q-q-h*s.Rational(5,4))==0
  assert s.simplify(q/c*s.Rational(1,5)-1)==0 and float(rho)>0
  dp=s.diff(b*b*F,b);assert float(dp.subs(b,s.Rational(5,4)))>0

def test_winding_lift_has_positive_certified_rayleigh_bounds_and_convergence():
 from w33_pass11341_winding_lifted_wall import shape_spectrum
 r=cert(11341);lo,hi=map(s.Rational,r['exact_smallest_shape_bounds']);assert lo>0 and float(lo)<r['shape_eigenvalues'][0]<float(hi)
 assert max(abs(shape_spectrum(16)[:4]-shape_spectrum(24)[:4]))<1e-5

def test_stronger_physical_background_closes_exact_isotropic_af_branch():
 r=cert(11340)['stronger_exact_barrier'];k,w,c=[s.Rational(r[a]) for a in ['Riccati_coefficient','linear_coefficient','constant']];assert c-w*w/(4*k)==s.Rational(r['exact_positive_global_minimum'])>0
 assert all(s.Rational(x)>0 for x in r['positive_principal_minors'])
 a2=(1-1/s.sqrt(10))/6;b2=(1+2/s.sqrt(10))/6;assert s.simplify(4*a2+2*b2-1)==0 and s.simplify(4*a2*a2+2*b2*b2-s.Rational(1,5))==0
