#!/usr/bin/env python3
"""Actual mediator Yukawa contractions, coupled running and scalar-closure audit."""
from pathlib import Path
import sys,json
import numpy as np
import sympy as s
from scipy.sparse import csr_matrix, bmat
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_global_e6_cartan_covariants as I

def yukawas(y1,y2):
 _,d,_=I.tensors();z=csr_matrix((81,81),dtype=complex);out=[]
 for a in range(27):
  B=csr_matrix(y1*np.kron(d[:,:,a],np.eye(3))/np.sqrt(2));Y=bmat([[z,B,z],[B.T,z,z],[z,z,z]],format='csr');out.extend([Y,1j*Y])
 for i in range(3):
  for j in range(i,3):
   T=np.zeros((3,3));T[i,j]=T[j,i]=1 if i==j else 1/np.sqrt(2)
   B=csr_matrix(y2*np.kron(np.eye(27),T)/np.sqrt(2));Y=bmat([[z,z,B],[z,z,z],[B.T,z,z]],format='csr');out.extend([Y,-1j*Y])
 return out

def coefficient(y1,y2,a):
 Ys=yukawas(y1,y2);Y=Ys[a];W=sum((B.conj().T@B for B in Ys),start=csr_matrix((243,243),dtype=complex));beta=(W.conj()@Y+Y@W)/2
 for B in Ys:
  beta+=2*(B@Y.conj().T@B)+B*float((B.conj().multiply(Y)).sum().real)
 norm=float(Y.conj().multiply(Y).sum().real);c=complex(Y.conj().multiply(beta).sum()/norm)
 error=float(np.sqrt(abs((beta-c*Y).conj().multiply(beta-c*Y).sum())));assert error<1e-9
 return c.real,error

def flow(t,x):
 e,a,f,u,v=x;L=16*np.pi**2
 return np.array([-2*e**3,-8*a**3,(68/3)*f**3,u*(40*u*u+v*v-52*e*e-8*f*f),v*(5*u*u+29*v*v-52*e*e-8*f*f)])/L

def payload():
 checks=[]
 for u,v in [(1.,1.),(1.,2.),(2.,1.)]:
  a,err=coefficient(u,v,0);b,err2=coefficient(u,v,54);assert abs(a-(40*u*u+v*v))<1e-9 and abs(b-(5*u*u+29*v*v))<1e-9;checks.append({'y1':u,'y2':v,'beta1_over_y1':a,'beta2_over_y2':b,'tensor_closure_error':max(err,err2)})
 sol=solve_ivp(flow,(0,300),[.4,.3,.05,.05,.06],rtol=1e-10,atol=1e-12,dense_output=True);assert sol.success
 rows=[{'log_scale':t,'couplings':sol.sol(t).tolist()} for t in (0,10,100,300)]
 x,y=s.symbols('x y');f1=124*x*x+s.Rational(224,5)*x*y+s.Rational(159,50)*y*y-104*x+18;f2=24*x*y+s.Rational(127,5)*y*y-104*y+60
 resultant=s.resultant(f1,f2,x);real_roots=int(s.Poly(resultant,y).count_roots(-s.oo,s.oo));assert real_roots==0;roots=[]
 for yy in s.nroots(resultant,maxsteps=200):
  if abs(s.im(yy))>1e-10:continue
  yy=float(s.re(yy));xx=(104*yy-s.Rational(127,5)*yy*yy-60)/(24*yy);xx=float(xx)
  roots.append({'lambda1_over_gSO_squared':xx,'lambda2_over_gSO_squared':yy,'bounded_on_S_sector':bool(xx+yy/10>0 if yy>=0 else xx+yy*73/90>0)})
 return {'status':'PASS','result_scope':'PASS_COUPLED_GAUGE_YUKAWA_FLOW_AND_QUARTIC_CLOSURE_AUDIT',
 'channel_rank_obstruction':{'one_pair_matching':'C_ab=y1_a*y2_b/M for d(psi,psi,v_a) F_bdagger. Its2x2 channel matrix has rank at most1; the desired diagonal two-channel coefficient matrix has rank2. One vectorlike pair cannot independently supply both11294 sectors. The computed Yukawa running is the active single-channel truncation, with both sextets retained in the gauge inventory.', 'two_pair_E6_beta':-10, 'two_scalar_mediator_E6_beta':2, 'boundary':'A second vectorlike pair costs another12 and destroys E6 AF in this inventory. Two heavy(27,bar6) scalar mediators instead cost another6 and retain bE6=2; their complete Yukawa/portal flow is not computed here.'},
 'Yukawa_formula':'(16pi²)beta_Ya=1/2(W*Ya+YaW)+2sum_b Yb Ya† Yb+sum_b Yb ReTr(Yb†Ya)-3sum_g g²{C2F,Ya}; W=sum_b Yb†Yb. Real scalar/Weyl convention, actual243 fermions and66 real scalar Yukawa components.',
 'coupled_beta':'bE6=2,bSO10=8,bSU3_family=-68/3 for the two-sextet fermion-mediator inventory; beta_y1=y1(40y1²+y2²-52gE²-8gF²)/(16pi²), beta_y2=y2(5y1²+29y2²-52gE²-8gF²)/(16pi²).',
 'tensor_checks':checks,'flow_rows':rows,'family_Landau_log_scale':float(8*np.pi**2/((68/3)*.05**2)),
 'global_family_fixed_Yukawa_ratios':['40/33','50/33'],
 'pure_S_quartic_definition':'V_S=lambda1(TrS²)²/4+lambda2 TrS4/4, S real SO10 symmetric-traceless54 with canonical trace norm.',
 'necessary_pure_S_beta':'16pi² beta_l1=124l1²+(224/5)l1l2+(159/50)l2²-120gSO²l1+18gSO4; beta_l2=24l1l2+(127/5)l2²-120gSO²l2+60gSO4.',
 'gauge_quartic_identity':'sum_i<j(si-sj)^4=10TrS4+3(TrS²)² for tracelessS. Therefore a TrS4-only alignment coupling generates an independent double-trace counterterm.',
 'pure_S_fixed_ratio_controls':roots,'exact_real_resultant_root_count':real_roots,
 'closure_boundary':'The gauge/Yukawa contractions and their flow are computed for the named two Yukawas. Pure-S quartic contributions are mandatory but NOT the full quartic flow: K,S,U,V,A and family portals in11293 contribute additional counterterms. The original positive-square coefficient ansatz is not RG closed. No complete-coupling asymptotic-freedom claim is made.',
 'prior_owners':['analysis/w33_pass11293_yukawa_uv_budget.py','analysis/w33_pass11294_hierarchical_cp_phase_vacuum.py'],
 'primary_sources':['https://arxiv.org/abs/hep-ph/0211440']}
if __name__=='__main__':
 out=payload();(ROOT/'data/w33_pass11298_coupled_running_closure.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
