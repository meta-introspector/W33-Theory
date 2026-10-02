"""Independent algebraic, geometric and variational controls for five targets."""
from pathlib import Path
import json,sys
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
def cert(n,name):return json.loads((ROOT/f'data/w33_pass{n}_{name}.json').read_text())

def test_complex_phase_flow_uses_physical_invariant_not_fixed_matrix_phase():
 import w33_pass11325_complex_yukawa_phase_audit as Y
 a1,a2,b,theta=.8,1.1,.6,.73;A=np.diag([a1,a2]).astype(complex);B=b*np.exp(.5j*theta)*np.array([[0,1],[-1,0]]);fa,fb=Y.beta(A,B,73/3)
 actual=2*(fb[0,1]/B[0,1]).imag-(fa[0,0]/a1).imag-(fa[1,1]/a2).imag
 expected=np.sin(theta)*(12*a1*a2/5+b*b*(a1/a2+a2/a1));assert abs(actual-expected)<1e-12 and actual>0
 U=expm(1j*np.array([[.2,.3],[.3,-.1]]));aa,bb=Y.beta(U.T@A@U,U.T@B@U,73/3);assert max(np.max(abs(aa-U.T@fa@U)),np.max(abs(bb-U.T@fb@U)))<1e-11

def test_paired_larger_rays_have_exact_positive_scalar_bound():
 import sympy as s
 r=cert(11325,'complex_yukawa_phase_audit')
 for row in r['paired_larger_flavor_rays']:
  n=row['active_flavors'];x=s.Rational(row['a_squared']);w=s.Rational(row['b_squared']);k=s.Rational(73,3)
  assert (s.Rational(31,5)+n)*x+s.Rational(7,2)*w==k and s.Rational(21,5)*x+(s.Rational(11,2)+n)*w==k
  assert s.Rational(row['S_portal_barrier_exact_minimum'])>0

def test_noncommuting_local_vacuum_has_full_spectrum_and_cp_conjugate():
 import w33_pass11326_noncommuting_flavor_operator as F
 r=cert(11326,'noncommuting_flavor_operator')['vacua'][1];x=np.array(r['coordinates']);assert np.linalg.norm(F.gradient(x,1e9,h=5e-6))<1e-5
 assert r['family_orbit_rank']==8 and len(r['normal_Hessian_eigenvalues'])==16 and min(r['normal_Hessian_eigenvalues'])>.14
 assert min(np.diff(r['U_Gram_eigenvalues']))>.3 and min(np.diff(r['D_Gram_eigenvalues']))>.3
 assert r['commutator_norm']>1 and abs(r['CKM_type_Jarlskog'])>.09
 assert abs(r['CP_conjugate_J']+r['Im_trace_commutator_cubed'])<1e-12
 assert r['CP_conjugate_energy_error']<1e-10 and r['Hessian_step_control_error']<.01

def test_angular_bound_is_sharp_in_independent_fourier_control():
 import w33_pass11326_noncommuting_flavor_operator as F
 vals=np.array([0.,(3-np.sqrt(3))/6,(3+np.sqrt(3))/6]);A=np.diag(vals);W=np.exp(2j*np.pi*np.outer(np.arange(3),np.arange(3))/3)/np.sqrt(3);B=W@A@W.conj().T;C=A@B-B@A;J=np.trace(C@C@C).imag
 assert abs(F.J_NORMALIZATION*J*J/2**12-1)<1e-12
 # The cutoff denominator only decreases the normalized invariant.
 assert F.J_NORMALIZATION*J*J/2.1**12<1

def test_quotient_chart_covariance_under_another_coordinate_basis():
 from w33_pass11327_nonlinear_quotient_mass_maps import Chart
 C=Chart();x=np.random.default_rng(113270).normal(size=22)*.001;R=np.diag(np.exp(1j*np.linspace(.1,.8,11)));O=np.block([[R.real,-R.imag],[R.imag,R.real]]);D=Chart(N=C.N@R);y=O.T@x
 assert max(abs(C.retract(x)-D.retract(y)))<1e-11
 assert np.max(abs(D.metric(y)-O.T@C.metric(x)@O))<1e-8

def test_covariant_mass_maps_replay_exact_reference_and_nonzero_connection():
 r=cert(11327,'nonlinear_quotient_mass_maps');ev=np.linalg.eigvalsh(np.array(r['reference_scalar_mass_matrix']));expected=np.array([0.]*8+[162/49]*10+[486/49]*2+[486/7,648/7]);assert max(abs(ev-expected))<.001
 M=np.array(r['reference_Weyl_real'])+1j*np.array(r['reference_Weyl_imag']);fw=np.sort(np.linalg.svd(M,compute_uv=False)**2);target=np.array([81/49]*8+[324/49]*2+[81.]);assert max(abs(fw-target))<1e-4
 assert r['generic_moment_residual']<1e-12 and r['generic_connection_Hessian_correction_norm']>1e-6 and r['generic_Kahler_block_error']<1e-7

def test_rank_two_lapse_matrix_implies_many_constraints_when_sum_is_nonzero():
 X=np.array([.2,-.1,.7,1.3]);y=np.array([1.,2.,3.,4.]);M=X[:,None]-X[None,:];assert np.linalg.matrix_rank(M)==2
 assert np.max(abs((M@y)[1:]-(M@y)[0]-(X[1:]-X[0])*sum(y)))<1e-12
 r=cert(11328,'determinant_multivielbein_scope');assert r['native_case']['native_edges']==160 and r['native_case']['determinant_quadratic_edges']==3160
 assert (2560-2*14-1738)//2==397

def test_determinant_quadratic_form_in_random_symmetric_coframes():
 rng=np.random.default_rng(33);hs=rng.normal(size=(3,4,4));hs=(hs+hs.transpose(0,2,1))/2
 q=lambda h:np.trace(h@h)-np.trace(h)**2
 expected=9/2*sum(q(hs[i]-hs[j]) for i in range(3) for j in range(i))
 def v(t):return np.linalg.det(3*np.eye(4)+t*sum(hs))-27*sum(np.linalg.det(np.eye(4)+t*h) for h in hs)
 h=1e-4;actual=(v(h)+v(-h)-2*v(0))/(2*h*h);assert abs(actual-expected)<1e-5

def test_wall_shape_operator_variational_levels_and_exact_flat_modes():
 import sympy as s
 from w33_pass11329_wall_shape_fluctuations import spectrum
 r=cert(11329,'wall_shape_fluctuations');assert r['exact_zero_modes']['count']==2
 b=s.symbols('b');p=(16*b-b**4-15)/30;assert p.subs(b,1)==0 and s.diff(p,b).subs(b,s.Rational(5,4))>0
 for m in [0,1,2]:
  for parity in ['even','odd']:
   low=spectrum(m,parity,12)[:4];high=spectrum(m,parity,20)[:4];assert min(high)>-1e-9 and np.min(low-high)>-1e-7
 assert abs(spectrum(0,'even',20)[0])<1e-9 and spectrum(0,'odd',20)[0]>1

def test_nonlinear_maps_extend_prior_linear_jets_with_quadratic_remainders():
 r=cert(11327,'nonlinear_quotient_mass_maps')['prior_22_field_linear_jet_controls']
 assert [x['scale'] for x in r]==[1.,.5]
 for key in ['scalar_linear_jet_residual','Weyl_linear_jet_residual']:
  assert .18<r[1][key]/r[0][key]<.35
