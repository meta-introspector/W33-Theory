"""Independent controls for the five scoped physical frontier witnesses."""
import json,sys
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.linalg import sqrtm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
def cert(n,name):return json.loads((ROOT/f'data/w33_pass{n}_{name}.json').read_text())

def test_every_real_yukawa_support_solves_original_ratio_equations():
 d=cert(11320,'general_real_two_flavor_yukawas');assert len(d['branches'])==8
 M=sp.Matrix([[sp.Rational(41,20),sp.Rational(103,20),sp.Rational(7,2)],[sp.Rational(103,20),sp.Rational(41,20),sp.Rational(11,2)],[sp.Rational(21,20),sp.Rational(33,20),sp.Rational(15,2)]])
 for r in d['branches']:
  x=sp.Matrix([sp.Rational(r[k]) for k in ['Sigma_squared','Delta_squared','K_Yukawa_squared']]);f=M*x-sp.ones(3,1)*sp.Rational(73,3)
  assert all(x[i]*f[i]==0 for i in range(3)) and sp.Rational(r['exact_positive_Riccati_minimum'])>0
 assert d['new_noncommuting_branch']['commutator_norm']>1

def test_flavor_vacua_stationarity_normal_stability_and_relative_cp():
 import w33_pass11321_full_flavor_bifurcation as F
 d=cert(11321,'full_flavor_bifurcation')
 for r in d['vacua']:
  x=np.array(r['coordinates']);assert np.linalg.norm(F.gradient(x,r['sigma'],1e-5))<3e-6
  assert min(r['normal_Hessian_eigenvalues'])>.001 and r['Hessian_step_control_max_difference']<1e-5
  assert abs(r['Jarlskog_commutator_cubic'])<1e-18
  if r['sigma']>.36:assert abs(r['Im_trace_UdaggerD'])>.01 and r['orbit_rank']==7
 assert d['symmetric_branch_control_sigma037']['negative_modes_below_minus001']>=5

def test_cp_conjugate_vacuum_reverses_invariant_without_energy_bias():
 import w33_pass11321_full_flavor_bifurcation as F
 r=cert(11321,'full_flavor_bifurcation')['vacua'][1];x=np.array(r['coordinates']);U,D=F.fields(x);y=F.pack(U.conj(),D.conj())
 assert abs(F.val(x,r['sigma'])-F.val(y,r['sigma']))<1e-11
 assert abs(np.trace(U.T@D.conj()).imag+r['Im_trace_UdaggerD'])<1e-12

def test_hard_matching_is_regulator_free_but_does_not_erase_soft_ir():
 import mpmath as mp
 d=cert(11322,'hard_soft_invariant_matching');assert 'regulator' not in d
 mp.mp.dps=90
 for r in d['counterfunctions']:
  p=list(map(mp.mpf,r['P_coefficients_ascending']))
  if r['soft_derivatives_set_to_zero']:assert abs(p[1])<mp.mpf('1e-60') and abs(2*p[2])<mp.mpf('1e-60')
  c=mp.mpf(5)/6 if r['kind']=='vector_hard' else mp.mpf('1.5')
  for x in r['hard_nodes']:
   x=mp.mpf(str(x));v=sum(j*p[j]*x**(j-1) for j in range(1,len(p)));expected=2*x*(mp.log(x/mp.mpf('.01'))-c+mp.mpf('.5'));assert abs(v-expected)<mp.mpf('1e-40')
 assert d['soft_IR_log_coefficient_rank']==11

def test_pair_metric_square_root_matches_named_adm_potential():
 def metric(N,s):
  g=np.eye(4);g[0,0]=-N*N+s*s;g[0,1]=g[1,0]=s;return g
 Ni,Nj,si,sj=1.17,.93,.14,-.21
 actual=Ni*np.trace(sqrtm(np.linalg.solve(metric(Ni,si),metric(Nj,sj)))).real
 assert abs(actual-(np.sqrt((Ni+Nj)**2-(si-sj)**2)+2*Ni))<1e-12

def test_native_lapse_hessian_by_finite_difference_after_shift_elimination():
 import w33_pass11323_native_cycle_lapse_hessian as G
 r=cert(11323,'native_cycle_lapse_hessian');B,C=G.cycle_basis();N=np.array(r['lapses']);p=np.array(r['momentum_source']);H=np.array(r['lapse_Hessian']);v=np.random.default_rng(77).normal(size=80);v/=np.linalg.norm(v);h=2e-4
 def grad(N):
  j,_,_=G.solve_currents(B,C,N,p);return abs(B)@np.sqrt(1+j*j)
 diff=(grad(N+h*v)-grad(N-h*v))/(2*h)
 assert np.max(abs(diff-H@v))<1e-7
 assert np.linalg.matrix_rank(H,tol=1e-10)==78 and max(np.linalg.eigvalsh(H))<1e-10
 assert r['symmetry_obstruction']['native_flag_orbit_size']==160

def test_wall_quantization_and_both_junctions_are_exact():
 import w33_pass11324_quantized_warped_history_wall as W
 r=W.exact_data();b=sp.symbols('b',positive=True);F=sp.Rational(8,15)/b-b*b/30-1/(2*b*b);bw=sp.Rational(5,4)
 assert F.subs(b,1)==0 and sp.diff(F,b).subs(b,1)==sp.Rational(2,5)
 assert (sp.diff(F,b)-2*F/b).subs(b,bw)==0
 assert sp.Rational(r['tension_squared'])==16*F.subs(b,bw)/bw**2
 assert sp.simplify(sp.sympify(r['magnetic_flux'])/(4*sp.pi))==1
 assert sp.Rational(r['naive_radial_potential_second_derivative'])==-sp.Rational(512,625)
 assert sp.Rational(r['other_spatial_junction_linear_mismatch'])!=0

def test_one_handle_map_is_explicit_primitive_native_quotient():
 from w33_pass11289_cycle_gram_gluing_flux import cycle_basis
 r=cert(11324,'quantized_warped_history_wall')['native_one_handle_quotient'];_,C=cycle_basis();assert np.array_equal(np.array(r['cochain'])@C,np.r_[1,np.zeros(80,dtype=int)])
