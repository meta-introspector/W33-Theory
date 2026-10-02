#!/usr/bin/env python3
"""A declared one-parameter hierarchical CP potential with relative phase locking."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11286_misaligned_family_vacua as F

def mixing(e):
 a,b,c=e,e*e,e**3;A,B,C=np.sqrt(1-np.array([a,b,c])**2);z=1j
 return np.array([[A*C,a*C,c/z],[-a*B-A*b*c*z,A*B-a*b*c*z,b*C],[a*b-A*B*c*z,-A*b-a*B*c*z,B*C]])

def vacuum(e):
 V=mixing(e);return np.diag([1.,2.,3.]).astype(complex),V@np.diag([1.,3.,5.])@V.T

def phase_invariants(U,D):
 Z=np.linalg.det(U)*np.linalg.inv(U)@D
 return np.array([np.trace(Z),np.trace(Z@Z)])

def constraints(U,D,e):
 out=list(F.constraints(U,D)[:10]);A=U@U.conj().T;B=D@D.conj().T;Ur,Dr=vacuum(e);Ar=Ur@Ur.conj().T;Br=Dr@Dr.conj().T
 for p in (1,2):
  for q in (1,2):out.append((np.trace(np.linalg.matrix_power(A,p)@np.linalg.matrix_power(B,q))-np.trace(np.linalg.matrix_power(Ar,p)@np.linalg.matrix_power(Br,q))).real)
 target=phase_invariants(Ur,Dr)
 for z,t in zip(phase_invariants(U,D),target):
  # A real-coefficient polynomial with CP-conjugate roots. Positive norm is CP even.
  v=z*z-2*t.real*z+abs(t)**2
  out.extend([v.real,v.imag])
 return np.array(out)

def phase_jacobian(e):
 U,D=vacuum(e);V=mixing(e);dirs=[np.diag([1.,0.,-1.]),np.diag([0.,1.,-1.])];cols=[];h=1e-5
 for sector in range(2):
  for T in dirs:
   if sector==0:
    dz=1j*U@T;du=dz;dd=np.zeros((3,3),complex)
   else:du=np.zeros((3,3),complex);dd=1j*V@np.diag([1.,3.,5.])@T@V.T
   plus=phase_invariants(U+h*du,D+h*dd);minus=phase_invariants(U-h*du,D-h*dd);der=(plus-minus)/(2*h);cols.append(np.r_[der.real,der.imag])
 return np.array(cols).T

def payload():
 rows=[]
 for e in (.18,.22,.3,.4):
  U,D=vacuum(e);V=mixing(e);assert np.max(abs(V.conj().T@V-np.eye(3)))<1e-12
  assert abs(np.linalg.det(V)-1)<1e-12
  residual=max(abs(constraints(U,D,e)));assert residual<1e-7
  j=F.jarlskog(V);exact=e**12*(1-e**2)*(1-e**4)*(1-e**6)**2;assert abs(j*j-exact)<1e-15
  P=phase_jacobian(e);sv=np.linalg.svd(P,compute_uv=False)
  cols=[];h=2e-6
  for sector in range(2):
   for a in F.directions():
    plus=(U+h*a,D) if sector==0 else (U,D+h*a);minus=(U-h*a,D) if sector==0 else (U,D-h*a)
    cols.append((constraints(*plus,e)-constraints(*minus,e))/(2*h))
  J=np.array(cols).T;sc=np.max(abs(J),axis=1);sc[sc<1e-12]=1;sing=np.linalg.svd(J/sc[:,None],compute_uv=False);rank=int(sum(sing>1e-7))
  assert rank==16 and sv[-1]>1e-8
  assert max(abs(constraints(U.conj(),D.conj(),e)))<1e-7
  rows.append({'epsilon':e,'constraint_residual':float(residual),'normal_rank':rank,'scaled_singular_values':sing.tolist(),'phase_jacobian_singular_values':sv.tolist(),'phase_jacobian_determinant':float(np.linalg.det(P)),'Jarlskog':j,'Jarlskog_squared':exact,'sin_squared_angles':[e**2,e**4,e**6],'phase_target_real':phase_invariants(U,D).real.tolist(),'phase_target_imag':phase_invariants(U,D).imag.tolist()})
 return {'status':'PASS','result_scope':'PASS_DECLARED_HIERARCHICAL_CP_VACUA_WITH_LOCAL_PHASE_STABILIZATION',
 'potential':'Prior spectral/determinant selectors; four mixed moment squares targeted to P(epsilon); plus |z_k²-2Re(t_k)z_k+|t_k|²|² for z1=Tr(adj(Fu)Fd), z2=Tr((adj(Fu)Fd)²). All coefficients real; positive squares are CP even.',
 'polynomial_scope':'Adj(Fu) is polynomial of degree2; inverse in evaluation is its nonsingular numerical representation. Maximum degree24. Reference moment and phase-root coefficients are declared inputs.',
 'invariant_proof':'For R in SU3, adj(R Fu R^T) R Fd R^T = R^-T adj(Fu)Fd R^T. Trace powers are invariant under this similarity. CP conjugates all z, and each phase polynomial has real coefficients.',
 'hierarchy':'s12=epsilon, s23=epsilon², s13=epsilon³, delta=pi/2; J²=epsilon^12(1-epsilon²)(1-epsilon^4)(1-epsilon^6)². This is a prescribed one-parameter deformation, not a W33-derived or measured CKM prediction.',
 'local_isolation':'24 real sextet fields, constraint Jacobian rank16 at all four scanned vacua;8 remaining directions are the SU3 orbit. Phase-invariant Jacobian on the previous four flat directions is nonsingular. Numerical controls, not an interval theorem for all epsilon.',
 'rows':rows,'boundaries':['Eigenvalue targets, epsilon powers and phase-root coefficients are imposed, not fitted or derived from W33.','Only local isolation at scanned points; other discrete minima and epsilon->0 singular limits are not classified.','Absolute masses, observed CP phase and RG stability remain open.'],
 'prior_owners':['analysis/w33_pass11286_misaligned_family_vacua.py','analysis/BT891_yukawa_texture_from_grading.md'],
 'primary_sources':['https://arxiv.org/abs/1103.2915']}
if __name__=='__main__':
 out=payload();(ROOT/'data/w33_pass11294_hierarchical_cp_phase_vacuum.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
