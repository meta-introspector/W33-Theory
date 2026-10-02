"""Independent mathematical controls and boundaries for the seven physical investigations."""
from pathlib import Path
import sys,json
import numpy as np
import sympy as s
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11313_symmetric_adjoint_portal_closure as P
import w33_pass11314_untargeted_phase_vacuum as F
import w33_pass11315_regulator_free_goldstone_poles as G
import w33_pass11316_native_edge_fierz_pauli as T
import w33_pass11318_radiative_scale_on_af_trajectory as R

def cert(n):return json.loads(next((ROOT/'data').glob(f'w33_pass{n}_*.json')).read_text())

def test_portal_Hessians_against_independent_directional_second_differences():
 rng=np.random.default_rng(11319);S,K=P.sample(rng);v=rng.normal(size=99);v/=np.linalg.norm(v);ds=np.einsum('a,aij->ij',v[:54],P.BS);dk=np.einsum('a,aij->ij',v[54:],P.BK)
 target=np.einsum('i,aij,j->a',v,P.hessians(S,K),v);h=.0005
 got=(P.invariants(S+h*ds,K+h*dk)-2*P.invariants(S,K)+P.invariants(S-h*ds,K-h*dk))/h**2
 assert max(abs(got-target))<2e-7

def test_exact_beta_projection_with_an_independent_rotated_background():
 C,D,proof=P.exact_coefficients();X=np.zeros((10,10));X[1,8]=.27;X[8,1]=-.27;O=expm(X)
 S,K=P.sample(np.random.default_rng(1131319));S=O@S@O.T;K=O@K@O.T;I=P.invariants(S,K);H=P.hessians(S,K)
 assert np.max(abs(np.einsum('a,aij->ij',I,C)-np.einsum('aij,bji->ab',H,H)/2))<1e-10
 V=P.vector_gram(S,K);assert abs(I@D-1.5*np.trace(V@V))<1e-10
 assert s.Rational(proof['evaluation_matrix_determinant'])!=0

def test_old_K_barrier_is_an_exact_square_completion():
 x,y=s.symbols('x y');z=x+y/5
 f=106*x*x+38*x*y+3*y*y-s.Rational(272,3)*x+18
 g=24*x*y+19*y*y-s.Rational(272,3)*y+12
 low=s.Rational(6571,62)*z*z-s.Rational(272,3)*z+s.Rational(102,5)
 assert s.expand(f+g/5-low)==s.expand(s.Rational(62,25)*(y+s.Rational(5,62)*z)**2)
 assert s.Rational(102,5)-s.Rational(272,3)**2/(4*s.Rational(6571,62))==s.Rational(298418,295695)>0
 for row in cert(11313)['exact_portal_barriers']['barriers']:
  A=s.Matrix([[s.Rational(x) for x in r] for r in row['positive_matrix']]);assert all(A[:i,:i].det()>0 for i in range(1,6))

def test_dual_Weyl_matrix_is_Hermitian_realification():
 S,K=P.sample(np.random.default_rng(99));eps=np.array([[0.,1.],[-1.,0.]]);M=np.kron(np.eye(2),S)+np.kron(eps,K)
 assert np.max(abs(M-M.T))<1e-12
 assert np.max(abs(np.linalg.eigvalsh(M)-np.repeat(np.linalg.eigvalsh(S-1j*K),2)))<1e-12
 assert abs((8.2+3.5)*730/351-27+8/3)<1e-12

def test_phase_minimum_is_CP_even_and_has_positive_local_curvature():
 t=np.array([0.,0.,0.,np.pi,0.]);v=np.array([.1,-.2,.3,-.1,.7]);assert abs(F.val(v)-F.val(-v))<1e-12
 grad,H=F.derivatives(t,.001);assert max(abs(grad))<1e-8 and np.linalg.eigvalsh(H)[0]>1e-5
 for delta in [.2,.8,1.5]:
  q=t.copy();q[4]=delta;assert F.val(q)>F.val(t)
 assert 'spectra/angles remain supplied' in cert(11314)['inputs']

def test_massless_Goldstone_limit_is_finite_away_from_zero_momentum():
 a=.0001;s0=.00033;exact=G.soft_factor(s0,a,0)
 errors=[abs(G.soft_factor(s0,a,d)-exact) for d in [1e-9,1e-10,1e-11]]
 assert errors[2]<errors[1]<errors[0] and abs(exact.imag-1)<1e-12
 assert G.soft_factor(1e-20,a,0).real>G.soft_factor(1e-10,a,0).real
 d=cert(11315);assert d['massless_bubble_rank']==11 and sum(r['multiplicity'] for r in d['rows'])==14

def test_spin_two_gauge_map_and_physical_quotients_at_another_momentum():
 k=np.array([np.sqrt(14.),1.,2.,3.]);E=T.operator(k);g=T.gauge(k)
 assert max(abs(E@g).flat)<1e-12
 assert np.linalg.matrix_rank(E,tol=1e-9)==4 and np.linalg.matrix_rank(g,tol=1e-9)==4
 rest=T.operator(np.array([2.,0.,0.,0.]),4.);assert np.linalg.matrix_rank(rest,tol=1e-9)==5
 assert cert(11316)['total_linear_polarizations']==397

def test_old_internal_tensor_square_does_not_supply_a_massless_operator():
 from w33_pass11289_cycle_gram_gluing_flux import complex_data
 _,_,lines,_=complex_data();N=np.zeros((40,40))
 for j,line in enumerate(lines):N[list(line),j]=1
 A=N@N.T-4*np.eye(40);L=12*np.eye(40)-A;w,E=np.linalg.eigh(L);E=E[:,abs(w-10)<1e-8];assert E.shape[1]==24
 b=np.random.default_rng(299).normal(size=(24,24));b=(b+b.T)/2;b-=np.trace(b)*np.eye(24)/24;h=E@b@E.T
 assert abs(np.trace(h))<1e-12 and max(abs(L@h+h@L-20*h).flat)<1e-12

def test_history_Kunneth_and_flat_cover_obstruction():
 betti=np.convolve([1,81,81,1],[1,1]).tolist();assert betti==[1,82,162,82,1]
 assert sum((-1)**i*b for i,b in enumerate(betti))==0 and betti[1]>4
 assert cert(11317)['betti']==betti

def test_neck_stress_from_independent_warped_radius_formulas():
 x,a=s.symbols('x a',positive=True);r=s.sqrt(x*x+a*a);rp=s.diff(r,x);rpp=s.diff(rp,x)
 rho=s.simplify((1-rp**2-2*r*rpp)/r**2);pr=s.simplify((rp**2-1)/r**2);pt=s.simplify(rpp/r)
 assert s.simplify(rho+pr+2*a*a/r**4)==0 and s.simplify(pt-a*a/r**4)==0
 d=cert(11319);assert 'not the global81' in d['scope'][0]

def test_radiative_scale_matches_supertrace_and_original_AF_ray():
 d=R.payload();assert d['B']>0 and d['mass_supertrace_beta_error']<1e-12
 sm=np.array(d['scalar_squared_masses']);vm=np.array(d['vector_squared_masses']);fm=np.array(d['Weyl_squared_masses'])
 B=(sm@sm+3*vm@vm-2*fm@fm)/(64*np.pi**2);assert abs(B-d['B'])<1e-14
 assert d['scalon_mass_over_v']<.02 and abs(d['vacuum_energy_over_v_fourth']+B/2)<1e-14
 assert np.max(abs(np.array(d['UV_trajectory_terminal_ratios'])-[.31514549843936723,-.12087307962152648]))<1e-9
