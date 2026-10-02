"""Untargeted one-loop phase minimization in the declared fixed-spectrum flavor restriction."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
import numpy as np
from scipy.optimize import minimize

def mix(e,d):
 a,b,c=e,e*e,e**3;A,B,C=np.sqrt(1-np.array([a,b,c])**2);z=np.exp(1j*d)
 return np.array([[A*C,a*C,c/z],[-a*B-A*b*c*z,A*B-a*b*c*z,b*C],[a*b-A*B*c*z,-A*b-a*B*c*z,B*C]])
def fields(t,e=.22):
 U=np.diag(np.array([1.,2.,3.])*np.exp(1j*np.r_[t[:2],-sum(t[:2])]))
 V=mix(e,t[4]);D=V@np.diag(np.array([1.,3.,5.])*np.exp(1j*np.r_[t[2:4],-sum(t[2:4])]))@V.T
 return U,D

def val(t):
 U,D=fields(t);v=0
 for k in (1.,2.):
  A=.2*(U+k*D).conj();M=np.block([[A,3*np.eye(3)],[3*np.eye(3),np.zeros((3,3))]]);w=np.linalg.eigvalsh(M.conj().T@M)
  v-=2*np.sum(w*w*(np.log(w/9)-1.5))/(64*np.pi**2)
 return v

def derivatives(t,h):
 n=len(t);I=np.eye(n);f=val(t);grad=np.array([(val(t+h*I[i])-val(t-h*I[i]))/(2*h) for i in range(n)])
 H=np.empty((n,n))
 for i in range(n):
  for j in range(n):
   H[i,j]=(val(t+h*I[i]+h*I[j])-val(t+h*I[i]-h*I[j])-val(t-h*I[i]+h*I[j])+val(t-h*I[i]-h*I[j]))/(4*h*h)
 return grad,H

def payload():
 rng=np.random.default_rng(11314);rows=[]
 for _ in range(64):
  q=minimize(val,rng.uniform(-np.pi,np.pi,5),method='BFGS',options={'gtol':1e-7,'maxiter':400});t=(q.x+np.pi)%(2*np.pi)-np.pi
  rows.append({'value':float(q.fun),'phases':t.tolist(),'optimizer_success':bool(q.success)})
 rows.sort(key=lambda x:x['value']);t=np.array([0.,0.,0.,np.pi,0.]);grad,H=derivatives(t,.002);grad2,H2=derivatives(t,.001)
 ev=np.linalg.eigvalsh(H2);assert max(abs(grad2))<1e-8 and ev[0]>1e-6
 U,D=fields(t);Z=np.linalg.det(U)*np.linalg.inv(U)@D;cp_error=max(abs(np.array([np.trace(Z),np.trace(Z@Z)]).imag));assert cp_error<1e-10
 cpcontrols=[]
 for delta in (0.,.3,1.,np.pi/2):
  v=t.copy();v[4]=delta;cpcontrols.append({'delta':delta,'potential':float(val(v))})
 assert val(t)<=rows[0]['value']+1e-8
 return {'status':'PASS','result_scope':'PASS_UNTARGETED_CP_CONSERVING_LOCAL_PHASE_MINIMUM_IN_FIXED_FLAVOR_RESTRICTION','action':'The11304 two real-coupling Weyl probe determinants, k=1,2, y=.2, M=3, mu²=9; no11294 phase-root potential included.','inputs':'Fu Takagi values1,2,3; Fd1,3,5; sin anglesepsilon,epsilon²,epsilon³ with epsilon=.22. Five coordinates are4 determinant-preserving Takagi phases and variableCKM delta. These spectra/angles remain supplied.','reference_phases':t.tolist(),'reference_value':float(val(t)),'phase_gradient_max':float(max(abs(grad2))),'phase_Hessian':H2.tolist(),'phase_Hessian_eigenvalues':ev.tolist(),'step_halving_Hessian_error':float(max(abs(H-H2).flat)),'CP_invariant_imaginary_error':float(cp_error),'CKM_Jarlskog':0.,'delta_displacements':cpcontrols,'multistart_results':rows,'scope':['Positive local5-coordinate Hessian is a finite-difference witness, not a certified global minimum or positive24-field vacuum.','The64-start search is not exhaustive. It finds no evidence here for radiatively selected nonzeroCKM CP; phase sensitivity alone did not deliver CP violation.','A specifiedMS matching potential is used. Other allowed finite flavor counterterms can change this vacuum.','Observed hierarchy and angles remain imposed; no mass or mixing prediction is asserted.'],'prior_owners':['analysis/w33_pass11304_singlet_phase_probes.py','analysis/w33_pass11294_hierarchical_cp_phase_vacuum.py'],'primary_sources':['https://doi.org/10.1103/PhysRevD.7.1888']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11314_untargeted_phase_vacuum.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['phase_Hessian_eigenvalues'])
