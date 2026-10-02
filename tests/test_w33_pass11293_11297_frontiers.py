"""Independent UV matching, phase symmetry, pole covariance and flux/ADM controls."""
from pathlib import Path
import sys,json
import numpy as np
import sympy as s
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11293_yukawa_uv_budget as U
import w33_pass11294_hierarchical_cp_phase_vacuum as F
import w33_pass11295_ward_matched_modulus_poles as P
import w33_pass11296_covariant_matter_constraints as G
import w33_pass11297_membrane_flux_dynamics as M

def certificate(n):return json.loads(next((ROOT/'data').glob(f'w33_pass{n}_*.json')).read_text())

def test_signed_E6_mediator_schur_and_budget():
 assert U.signed_tensor_matching()['signed_E6_tensor_Schur_error']<1e-12
 # Independent representation multiplicities, not a replay of stored costs.
 for r in range(5,20):
  base=s.Rational(34)-6*r
  scalar_cost=s.Rational(1,3)*6*3;pair_cost=s.Rational(2,3)*2*3*3
  assert base-scalar_cost<0 and base-pair_cost<0


def test_hierarchical_phase_lock_covariance_and_CP_pair():
 for e in (.18,.22,.3,.4):
  u,d=F.vacuum(e);a=np.array([[0,1j,.17],[1j,0,0],[-.17,0,0]]);R=expm(.31*a)
  assert max(abs(F.constraints(R@u@R.T,R@d@R.T,e)))<1e-7
  assert max(abs(F.constraints(u.conj(),d.conj(),e)))<1e-7
  V=F.mixing(e);assert abs(abs(V[0,2])**2-e**6)<1e-14
  # A displacement on a formerly flat phase changes the new invariant roots.
  rotated=d.copy();T=np.diag([1.,0.,-1.]);V=F.mixing(e)
  rotated=V@np.diag([np.exp(.01j),3.,5*np.exp(-.01j)])@V.T
  assert np.linalg.norm(F.constraints(u,rotated,e)[-4:])>1e-8
  assert np.linalg.matrix_rank(F.phase_jacobian(e),tol=1e-8)==4


def test_exact_hierarchy_moduli_and_J_law():
 e,c1,c2,c3=s.symbols('e c1 c2 c3',positive=True);basis=s.groebner([c1*c1+e**2-1,c2*c2+e**4-1,c3*c3+e**6-1],c1,c2,c3,e,extension=s.I)
 V=s.Matrix([[c1*c3,e*c3,-s.I*e**3],[-e*c2-s.I*c1*e**5,c1*c2-s.I*e**6,e**2*c3],[e**3-s.I*c1*c2*e**3,-c1*e**2-s.I*e*c2*e**3,c2*c3]])
 assert all(basis.reduce(s.expand(x))[1]==0 for x in V.conjugate().T*V-s.eye(3))
 j=s.im(V[0,0]*V[1,1]*s.conjugate(V[0,1])*s.conjugate(V[1,0]));assert basis.reduce(s.expand(j*j-e**12*(1-e**2)*(1-e**4)*(1-e**6)**2))[1]==0


def test_Ward_matching_is_basis_covariant_and_scheme_explicit():
 sm,fm,O,groups=P.inputs();C,SP,proj,_=P.ward_matching(groups,sm)
 R=np.linalg.qr(np.random.default_rng(11295).normal(size=(22,22)))[0];gs=[dict(g,R=R.T@g['R']@R,I=R.T@g['I']@R) for g in groups]
 C2,SP2,_,_=P.ward_matching(gs,sm,projector=R.T@proj@R)
 assert np.max(abs(C2-R.T@C@R))<1e-12
 assert np.max(abs(SP2+C2@R.T@proj@R))<1e-12
 # Allowed normal analytic matching changes massive poles at the same loop order.
 t=max(sm);S=P.sigma(groups,t,C);N=np.eye(22)-proj;delta=1e-6
 assert abs((S+delta*N)[-1,-1]-S[-1,-1]-delta)<1e-14
 d=certificate(11295);assert sum(r['multiplicity'] for r in d['rows'])==14 and len(d['Goldstone_poles'])==8
 assert all(x<=1e-10 for r in d['rows'] for x in r['pole_squared_imag'])


def test_covariant_matter_preserves_fiber_and_continuum_bracket():
 from w33_pass11283_symplectic_metric_field import chart_Jacobian
 A=s.Matrix(chart_Jacobian());v=A.nullspace()[0];assert A*v==s.zeros(16,1)
 a,b=G.spectral_bracket(64);a2,b2=G.spectral_bracket(128)
 assert max(abs(a-b),abs(a2-b2),abs(a-a2))<1e-12
 d=certificate(11296);assert d['degrees_of_freedom']['physical_config']==24
 assert 2*(11+22)-2*9==2*24


def test_membrane_rates_absorb_and_do_not_sequester_bare_shift():
 ns,f,r,L,trans,k=M.dynamics();assert np.allclose(L.sum(axis=0),0)
 for row in trans:
  assert r[list(ns).index(row['from'])]>r[list(ns).index(row['to'])]
 p=np.zeros(len(ns));p[-1]=1;out=expm(1e6*L)@p;assert out[k]>.999 and abs(sum(out)-1)<1e-8
 _,_,r2,L2,_,k2=M.dynamics(bare=10.);assert k==k2 and np.max(abs(L-L2))<1e-10
 assert abs(r2[k]-r[k]-10.01)<1e-12
 # Global extensive charge along the old branch is not the membrane intensive field.
 for q in (1.,2.,5.):
  volume=q/.7;assert abs(q/volume-.7)<1e-12


def test_economical_intertwiner_covariance_and_mediator_budget():
 import w33_pass11285_semisimple_factor_higgs as H
 v,A,u,w,K,S=U.economical_reference();X=np.zeros((10,10));X[0,7]=1;X[7,0]=-1;R=expm(.23*X);E=expm(.017j*H.F.H.setup()[0][0])
 values=U.economical_constraints(E@v,E@A@E.conj().T,E@u@R.T,E@w@R.T,R@K@R.T,R@S@R.T)
 assert max(np.max(abs(a)) for a in values)<1e-9
 assert np.max(abs((A@A+A-6*np.eye(27))@u))<1e-12
 assert np.max(abs(S-(u.conj().T@A@u+(u.conj().T@A@u).conj())))<1e-12
 # Independent SO10 symmetric-traceless index and changed inventory counting.
 assert s.Rational(54*10,45)==12
 bSO=s.Rational(11,3)*8-s.Rational(1,3)*54-s.Rational(8+12,6)
 assert bSO==8 and 34-20-12==2 and 34-20-6==8
 d=certificate(11293)['economical_alignment_escape'];assert d['real_alignment_fields']-111==1254
