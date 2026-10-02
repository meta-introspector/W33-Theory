#!/usr/bin/env python3
"""Declared membrane relaxation and its incompatibility with a uniform linear flux constraint."""
from pathlib import Path
import json,math
import numpy as np
import sympy as s
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1]

def dynamics(q=1.,offset=.23,tension=.2,bare=-.01,extent=6):
 ns=np.arange(-extent,extent+1);f=q*ns+offset;rho=bare+f*f/2;L=np.zeros((len(ns),len(ns)));trans=[]
 for j,n in enumerate(ns):
  near=int(n-np.sign(f[j]));i=np.where(ns==near)[0]
  if not len(i):continue
  i=int(i[0]);delta=rho[j]-rho[i]
  if delta<=0:continue
  R=3*tension/delta;B=27*np.pi**2*tension**4/(2*delta**3);rate=np.exp(-B)
  L[i,j]=rate;L[j,j]-=rate;trans.append({'from':int(n),'to':int(near),'Delta_rho':float(delta),'bubble_radius':float(R),'bounce_action':float(B),'rate':float(rate)})
 terminal=np.where(np.diag(L)==0)[0];assert len(terminal)==1;assert np.max(abs(L.sum(axis=0)))<1e-14
 return ns,f,rho,L,trans,int(terminal[0])

def payload():
 ns,f,rho,L,trans,k=dynamics();initial=np.zeros(len(ns));initial[-1]=1;rows=[]
 for t in (0.,10.,1000.,1e6):
  prob=expm(t*L)@initial;assert min(prob)>-1e-10 and abs(sum(prob)-1)<1e-8
  rows.append({'time':t,'terminal_probability':float(prob[k]),'mean_flux':float(prob@f),'mean_energy':float(prob@rho)})
 assert rows[-1]['terminal_probability']>.999
 Q,mu,e=s.symbols('Q mu e',positive=True);vol=Q/mu;intensive=s.simplify(Q/vol);density=s.simplify((Q**2/(2*e**2*vol))/vol);assert intensive==mu and s.diff(density,Q)==0
 scans=[]
 for bare in (-1.,-.01,.01,1.):
  _,_,r,_,_,j=dynamics(bare=bare);scans.append({'bare_energy':bare,'terminal_energy':float(r[j])})
 assert abs(scans[-1]['terminal_energy']-scans[0]['terminal_energy']-2)<1e-12
 # R extrema of the flat-space thin-wall action are classical, imported, not new.
 R,T,D=s.symbols('R T Delta',positive=True);bounce=2*s.pi**2*T*R**3-s.pi**2*D*R**4/2
 assert s.simplify(s.diff(bounce,R).subs(R,3*T/D))==0
 assert s.simplify(bounce.subs(R,3*T/D)-27*s.pi**2*T**4/(2*D**3))==0
 return {'status':'PASS','result_scope':'PASS_DECLARED_MEMBRANE_MARKOV_RELAXATION_AND_UNIFORM_FLUX_CONSTRAINT_AUDIT',
 'action':'Added compact3-form with intensive field f=*F4, Maxwell energy f²/2, charged membranes of chargeq and tensionT. Homogeneous sectors f_n=offset+nq. Local membrane jumps change f byq; these are additional inputs.',
 'transition_rule':'Zero-temperature flat-space thin-wall downward transitions with Delta_rho>0, R=3T/Delta_rho, B=27pi²T⁴/(2Delta_rho³), rate=A exp(-B), A=1 input. Upward/gravitational/cosmological-volume corrections are omitted.',
 'parameters':{'charge':1.,'offset':.23,'tension':.2,'bare_energy':-.01,'prefactor':1.},'transitions':trans,'evolution':rows,
 'terminal_sector':int(ns[k]),'terminal_energy':float(rho[k]),'terminal_intensive_flux':float(f[k]),'bare_shift_scan':scans,
 'selection':'The finite downward chain has one absorbing sector minimizing|offset+nq|. The actual evolution reaches it fromn=6. Residual energy=bare+distance(offset,qZ)²/2, with charge/offset/bare inputs. Bare vacuum-energy shifts move the terminal energy unchanged; this is not sequestering.',
 'extensive_intensive_audit':'Prior11284/11289 use extensiveQ=integralF and a uniform linear constraint Vol=Q/mu4. On that homogeneous branch f=Q/Vol=mu4 and Maxwell density=mu4²/(2e²), independent ofQ. VaryingQ along that branch is not the intensive ladder used by Brown-Teitelboim relaxation. Do not identify these two flux variables.',
 'sequestering_boundary':'This excludes naively combining the uniform linear-volume branch with the homogeneous membrane ladder. It does not forbid inhomogeneous bubbles with a global constraint, nonlinear sigma, or additional fluxes; those require solving their actual junction and global equations.',
 'boundaries':['Known thin-wall bounce and vacuum-relaxation mechanism applied to a declared model, not a novel mechanism or observed CC prediction.','No gravity/CDL correction, volume-weighted measure, prefactor derivation, upward transitions or membrane spectrum from W33.','Topology provides a flux class, not its charge, offset, tension or selected observed vacuum energy.'],
 'prior_owners':['analysis/w33_pass11284_closed_history_flux.py','analysis/w33_pass11289_cycle_gram_gluing_flux.py'],
 'primary_sources':['https://arxiv.org/abs/1804.09985','https://arxiv.org/abs/2306.09412','https://arxiv.org/html/1505.01492v2']}
if __name__=='__main__':
 out=payload();(ROOT/'data/w33_pass11297_membrane_flux_dynamics.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
