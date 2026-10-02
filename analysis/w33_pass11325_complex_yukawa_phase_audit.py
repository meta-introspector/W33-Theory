"""Complex symmetric/antisymmetric SO10 Yukawa RG, including physical phase flow."""
from pathlib import Path
import json,sys
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11313_symmetric_adjoint_portal_closure as P

def beta(A,B,k=0.):
 H=5.4*A.conj().T@A+4.5*B.conj().T@B
 fa=(H.T@A+A@H)/2+.8*A@A.conj().T@A+B@A.conj().T@B+np.trace(A.conj().T@A).real*A-k*A
 fb=(H.T@B+B@H)/2+B@B.conj().T@B+1.2*A@B.conj().T@A+np.trace(B.conj().T@B).real*B-k*B
 return fa,fb

def direct_control(n):
 rng=np.random.default_rng(11325+n);A=(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))*.1;A=(A+A.T)/2;B=(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))*.1;B=(B-B.T)/2
 Y=np.array([np.kron(A,x) for x in P.BS]+[np.kron(B,x) for x in P.BK]);H=np.einsum('aji,ajk->ik',Y.conj(),Y);fa,fb=beta(A,B);error=0.
 for ix in [0,17,54,65]:
  a=Y[ix];v=(H.T@a+a@H)/2+2*np.einsum('aij,jk,akl->il',Y,a.conj().T,Y)
  tr=np.real(np.einsum('aij,ij->a',Y.conj(),a));v+=np.einsum('a,aij->ij',tr,Y)
  target=np.kron(fa,P.BS[ix]) if ix<54 else np.kron(fb,P.BK[ix-54]);error=max(error,float(max(abs(v-target).flat)))
 assert error<1e-11;return error

def paired_rays():
 rows=[];k=s.Rational(73,3)
 for n in [2,4,6,8]:
  M=s.Matrix([[s.Rational(31,5)+n,s.Rational(7,2)],[s.Rational(21,5),s.Rational(11,2)+n]])
  x,w=M.inv()*s.Matrix([k,k]);gamma=n*x;box=n*x*x;minimum=s.factor(36-s.Rational(6,5)*box-(120-4*gamma-s.Rational(16,3))**2/(4*s.Rational(16580,159)));assert minimum>0
  J=np.kron(np.eye(n//2),np.array([[0,1],[-1,0]]));A=np.sqrt(float(x))*np.eye(n);B=np.sqrt(float(w))*J;fa,fb=beta(A,B,float(k));res=max(max(abs(fa).flat),max(abs(fb).flat));assert res<1e-11
  rows.append({'active_flavors':n,'a_squared':str(x),'b_squared':str(w),'S_portal_barrier_exact_minimum':str(minimum),'ratio_residual':float(res)})
 return rows

def payload():
 controls={str(n):direct_control(n) for n in [2,4]}
 return {'status':'PASS','result_scope':'PASS_COMPLEX_TWO_FLAVOR_PHASE_LOCKING_AND_LARGER_PAIRED_AF_BARRIERS','matrix_beta':'H=(27/5)A†A+(9/2)B†B; betaA=(H^T A+A H)/2+(4/5)AA†A+BA†B+Tr(A†A)A-kA; betaB=(H^T B+B H)/2+BB†B+(6/5)AB†A+Tr(B†B)B-kB; k=73/3 in ratio flow.','direct_tensor_controls':controls,'two_flavor_phase':'Takagi A=diag(a1,a2)>=0, B=b exp(i theta/2)epsilon. theta=arg[bcomplex²/(a1complex*a2complex)] is unchanged by unitary flavor congruence. For full rank, beta_theta=sin(theta)[(12/5)a1*a2+|b|²(a1/a2+a2/a1)]. The bracket is strictly positive. Thus any physical fixed phase has theta=0 orpi.','phase_classification':'theta0/pi are equivalent by flavor congruence to real A,B, hence reuse all eight real11320 supports. Rank-one A and nonzeroB cannot be stationary modulo flavor rotations: betaA in the zero-singular diagonal is nonzero. A0 orB0 have removable phases. This closes the complex two-active-flavor fixed-ray loophole, not every larger texture.','paired_larger_flavor_rays':paired_rays(),'scope':['All complex two-active Yukawa relative fixed rays reduce to the previously audited real family and retain scalar AF obstruction. This is a fixed-ray theorem, not an all-trajectory theorem.','Four/six/eight-active results cover A=aI and B=b blockepsilon only, with the unchanged eight-vector Weyl gauge budget. All displayed mixed rays retain positive S barriers. Arbitrary larger complex textures remain open.','A hypothesized complex rotating ray failed the producer assertion: its phase-flow sign was wrong. The corrected full tensor contraction and phase equation exclude it. No new complex fixed point is asserted.'],'prior_owners':['analysis/w33_pass11320_general_real_two_flavor_yukawas.py','analysis/w33_pass11313_symmetric_adjoint_portal_closure.py'],'primary_sources':['https://arxiv.org/abs/hep-ph/0211440']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11325_complex_yukawa_phase_audit.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['paired_larger_flavor_rays'])
