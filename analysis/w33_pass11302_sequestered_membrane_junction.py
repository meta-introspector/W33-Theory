#!/usr/bin/env python3
"""A named three-form separation and closed two-cap junction/global-constraint laboratory."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import root
ROOT=Path(__file__).resolve().parents[1]

def caps(R,Lambda,G,bare=0.):
 rho=np.array([Lambda+bare+.5*1.2**2,Lambda+bare+.5*.8**2]);H2=rho/(3*G)
 if min(H2)<=0 or max(H2)*R*R>=1:raise ValueError('outside positive-curvature two-cap branch')
 c=np.array([-1.,1.])*np.sqrt(1-H2*R*R);vol=2*np.pi**2/H2**2*(2/3-c+c**3/3);area=2*np.pi**2*R**3
 average=(4*rho@vol+3*.4*area)/(G*sum(vol));israel=(c.sum()/R)-.4/(2*G)
 return vol,average,israel

def inputs():
 L=.1;G=1.;tau=.4/(2*G);hp=(L+.5*1.2**2)/(3*G);hm=(L+.5*.8**2)/(3*G)
 R=1/np.sqrt(hm+(hp-hm+tau*tau)**2/(4*tau*tau));vol,avg,j=caps(R,L,G);assert abs(j)<1e-12
 mu=.7;Q=mu*sum(vol);Qhat=-Q*avg/(2*mu)
 return R,L,G,mu,Q,Qhat

def equations(x,mu,Q,Qhat,bare=0.):
 R,L,G=x
 try:vol,avg,j=caps(R,L,G,bare)
 except ValueError:return np.array([1e3,1e3,1e3])
 return np.array([j,(sum(vol)-Q/mu)/(Q/mu),(avg+2*mu*Qhat/Q)/abs(2*mu*Qhat/Q)])

def payload():
 R,L,G,mu,Q,Qhat=inputs();sol=root(equations,[R*.99,L+.01,G*1.01],args=(mu,Q,Qhat),tol=1e-10);assert sol.success and max(abs(equations(sol.x,mu,Q,Qhat)))<1e-9
 vol,avg,j=caps(*sol.x);shift=.03;other=sol.x.copy();other[1]-=shift;error=max(abs(equations(other,mu,Q,Qhat,bare=shift)));assert error<1e-9
 return {'status':'PASS','result_scope':'PASS_DECLARED_SEPARATE_MEMBRANE_FORM_JUNCTION_AND_GLOBAL_SEQUESTERING_SOLUTION',
 'named_action':'Original rigid-Lambda,kappa² sequestering with topological Fseq and Fhat, plus a distinct Maxwell four-form Fmem and charged membrane T*area-q*integral A_mem. Never put the membrane charge on Fseq in the linear-sigma branch.',
 'local_form_equations':'Variation of local Lambda gives Fseq=mu4*volume_form for linear sigma. Its extensive Qseq imposes Vol=Qseq/mu4. Fseq cannot at the same time be a jumping homogeneous Maxwell membrane field. The added Fmem instead obeys fplus-fminus=q independently.',
 'laboratory_geometry':'Closed Euclidean two-cap S4 with outside large cap and inside small cap, both positive curvature. This added instanton geometry is NOT the genus81 M3xS1 W33 completion; no identification is made.',
 'junction':'Hplus²=rhoplus/(3kappa²), Hminus²=rhominus/(3kappa²); sqrt(1-Hminus²R²)-sqrt(1-Hplus²R²)=T R/(2kappa²).',
 'global_equations':'Vol=Qseq/mu4 and <R>=-2mu4 Qhat/(M0² Qseq), M0²=1. Integral curvature includes the membrane trace term3T*area/kappa² as well as bulk4rho/kappa².',
 'inputs':{'T':.4,'fplus':1.2,'fminus':.8,'charge':.4,'mu4':mu,'Qseq':float(Q),'Qhat':float(Qhat)},
 'solved':{'radius':float(sol.x[0]),'Lambda':float(sol.x[1]),'kappa_squared':float(sol.x[2]),'cap_volumes':vol.tolist(),'average_curvature':float(avg)},'residual':float(max(abs(equations(sol.x,mu,Q,Qhat)))),'bare_shift_control':{'shift':shift,'residual_with_Lambda_compensation':float(error)},
 'scope':'A nontrivial numerical solution for junction, volume and curvature-flux conditions with fixed input sectors generated from a physical reference and recovered from a perturbed start. A common bare vacuum shift is canceled by Lambda->Lambda-shift without changing geometry or flux sectors.',
 'boundaries':['Independent Fmem, its charge/offset/tension and the closed cap geometry are additional inputs.','Qseq,Qhat and the membrane extensive flux are not derived or quantized here; no observed CC selection.','No bounce-negative-mode, stability, nucleation-rate, Lorentzian cosmology or global genus81 instanton construction is proved.'],
 'prior_owners':['analysis/w33_pass11284_closed_history_flux.py','analysis/w33_pass11297_membrane_flux_dynamics.py'],
 'primary_sources':['https://arxiv.org/abs/1604.04000','https://arxiv.org/html/1505.01492v2']}
if __name__=='__main__':
 out=payload();(ROOT/'data/w33_pass11302_sequestered_membrane_junction.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
