"""Finite-momentum Goldstone limit for the full22-field matched bubble operator."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11295_ward_matched_modulus_poles as P

def soft_factor(s,a,delta):
 if delta==0:return complex((-np.log(s/a)-1-s/a)/np.pi,1.)
 f=lambda x:x*(1-x)
 points=[(1-np.sqrt(1-4*delta/s))/2,(1+np.sqrt(1-4*delta/s))/2] if 4*delta<s else []
 real=-quad(lambda x:np.log(abs(delta-s*f(x))),0,1,points=points,epsabs=1e-11)[0]
 real+=quad(lambda x:np.log(delta+a*f(x)),0,1,epsabs=1e-11)[0]
 real-=(s+a)*quad(lambda x:f(x)/(delta+a*f(x)),0,1,epsabs=1e-11)[0]
 return complex(real/np.pi,np.sqrt(max(0,1-4*delta/s)))

def payload():
 sm,fm,O,groups=P.inputs();C,SP,proj,skip=P.ward_matching(groups,sm);soft=[g for g in groups if g['kind']=='scalar' and g['a']==g['b']==0]
 R=sum((g['R'] for g in soft),start=np.zeros((22,22)));rank=int(np.linalg.matrix_rank(R,tol=1e-13));assert rank==11
 rows=[];a=.0001;control_error=0.
 for s in sorted(set(round(float(x),12) for x in sm if x>1e-10)):
  ix=np.where(abs(sm-s)<1e-10)[0];base=P.sigma(groups,s,C);analytic=R*soft_factor(s,a,0)
  numeric=sum((P.Q.dispersion(g,s,-a)+1j*P.Q.rho(g,s) for g in soft),start=np.zeros((22,22),complex));control_error=max(control_error,float(max(abs(analytic-numeric).flat)))
  hard=base-numeric;exact=hard+analytic;poles=s-np.linalg.eigvals(exact[np.ix_(ix,ix)]);controls=[]
  for delta in (1e-9,1e-10,1e-11):
   reg=hard+R*soft_factor(s,a,delta);error=float(np.linalg.norm(reg-exact));controls.append({'soft_mass_squared':delta,'matrix_error':error})
  assert controls[-1]['matrix_error']<controls[0]['matrix_error']
  rows.append({'tree_mass_squared':s,'multiplicity':len(ix),'regulator_free_leading_pole_squared_real':poles.real.tolist(),'regulator_free_leading_pole_squared_imag':poles.imag.tolist(),'regulated_controls':controls})
 assert control_error<1e-10
 return {'status':'PASS','result_scope':'PASS_REGULATOR_FREE_FINITE_MOMENTUM_GOLDSTONE_LIMIT_IN_MATCHED_22_FIELD_EFT','massless_bubble_rank':rank,'analytic_vs_spectral_error':control_error,'rows':rows,'formula':'For spectral coefficientR and s0=-a, twice-subtracted massless bubble is R[-log(s/a)-1-s/a+i pi]/pi, s>0. Independently replayed by regulated Feynman-parameter integral.','scope':['This takes the on-shell soft mass to zero at finite external momentum and reproduces the prior conditional14 leading massive poles. It does not make the zero-momentum normal-normal Hessian finite.','Matching fixes the same input local two-point Ward columns as11295; an invariant tadpole/counterterm functional and hard Goldstone resummation are still required for a full physical pole prediction.','No heavy vector or full UV threshold matching is added here;11305 physical vector cuts remain a separate input.','The rank11 IR tensor is inherited from the actual quotient cubic, not a radial-only surrogate.'],'prior_owners':['analysis/w33_pass11295_ward_matched_modulus_poles.py','analysis/w33_pass11305_vector_matrix_and_total_jet_audit.py'],'primary_sources':['https://arxiv.org/abs/1609.06977']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11315_regulator_free_goldstone_poles.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['massless_bubble_rank'],d['analytic_vs_spectral_error'])
