"""All real two-active-vector SO10 Yukawa fixed ratios and exact scalar barriers."""
from pathlib import Path
import sys,json,itertools
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11313_symmetric_adjoint_portal_closure as P

def flavor_beta(A,B,gauge=0.):
 A2=A@A;B2=B@B
 fa=6.2*A2@A-2.25*(B2@A+A@B2)+B@A@B+np.trace(A2)*A-gauge*A
 fb=-5.5*B2@B+2.7*(A2@B+B@A2)-1.2*A@B@A-np.trace(B2)*B-gauge*B
 return fa,fb

def contraction_control():
 A=np.array([[.5,.2],[.2,-.1]]);B=np.array([[0.,.4],[-.4,0.]])
 Y=np.array([np.kron(A,x) for x in P.BS]+[np.kron(B,x) for x in P.BK]);Y2=np.einsum('aji,ajk->ik',Y,Y);fa,fb=flavor_beta(A,B);error=0.
 for ix in [0,17,54,65]:
  a=Y[ix];val=(Y2.T@a+a@Y2)/2+2*np.einsum('aij,jk,akl->il',Y,a.T,Y)
  val+=np.einsum('a,aij->ij',np.einsum('ij,aji->a',a,Y),Y)
  target=np.kron(fa,P.BS[ix]) if ix<54 else np.kron(fb,P.BK[ix-54]);error=max(error,float(np.max(abs(val-target))))
 assert error<1e-11
 return error

def branches():
 k=s.Rational(73,3);M=s.Matrix([[s.Rational(41,20),s.Rational(103,20),s.Rational(7,2)],[s.Rational(103,20),s.Rational(41,20),s.Rational(11,2)],[s.Rational(21,20),s.Rational(33,20),s.Rational(15,2)]])
 rows=[]
 for mask in itertools.product([0,1],repeat=3):
  ix=[i for i,b in enumerate(mask) if b];sol=M.extract(ix,ix).inv()*s.ones(len(ix),1)*k if ix else []
  u,v,w=[next((sol[j] for j,i in enumerate(ix) if i==l),s.S(0)) for l in range(3)];assert all(x>=0 for x in [u,v,w])
  gamma=(u+v)/2;box=(u*u+6*u*v+v*v)/8
  bound=s.factor(36-s.Rational(6,5)*box-(120-4*gamma-s.Rational(16,3))**2/(4*s.Rational(16580,159))) if w else s.Rational(298418,295695)
  assert bound>0
  rows.append({'support':list(mask),'Sigma_squared':str(u),'Delta_squared':str(v),'K_Yukawa_squared':str(w),'barrier_sector':'S' if w else 'K','exact_positive_Riccati_minimum':str(bound),'minimum_numeric':float(bound)})
 return rows

def payload():
 rows=branches();error=contraction_control();new=next(r for r in rows if r['support']==[1,1,1]);u,v,w=[float(s.Rational(new[k])) for k in ['Sigma_squared','Delta_squared','K_Yukawa_squared']];A=np.diag([(np.sqrt(u)+np.sqrt(v))/2,(np.sqrt(u)-np.sqrt(v))/2]);B=np.array([[0.,np.sqrt(w)],[-np.sqrt(w),0.]])
 fa,fb=flavor_beta(A,B,27-8/3);res=max(np.max(abs(fa)),np.max(abs(fb)));assert res<1e-11
 return {'status':'PASS','result_scope':'PASS_EXHAUSTIVE_REAL_TWO_ACTIVE_YUKAWA_FIXED_RAY_BARRIERS','normal_form':'A=A^T real2x2, B=-B^T real2x2; flavorO2 diagonalizesA and leaves B=b eps up to orientation. Sigma=a1+a2, Delta=a1-a2. This is the entire real two-active-flavor texture family, with six spectator vector Weyls and bSO=8/3.','beta_A':'(31/5)A³-(9/4){B²,A}+BAB+TrA² A-27g² A','beta_B':'-(11/2)B³+(27/10){A²,B}-(6/5)ABA-TrB² B-27g² B','fixed_ratio_equations':'Sigma[(41/20)Sigma²+(103/20)Delta²+(7/2)b²-73/3]=0; Delta[(103/20)Sigma²+(41/20)Delta²+(11/2)b²-73/3]=0; b[(21/20)Sigma²+(33/20)Delta²+(15/2)b²-73/3]=0. Every support reduces to a nonsingular rational linear system.','branches':rows,'direct_99_tensor_contraction_error':error,'new_noncommuting_branch':{'A':A.tolist(),'B':B.tolist(),'commutator_norm':float(np.linalg.norm(A@B-B@A)),'fixed_ratio_residual':float(res)},'fermion_box':'TrM4=TrA4 TrS4+TrB4 TrK4+4Tr(A²B²)Tr(S²K²)+2Tr(ABAB)Tr(SKSK), M=A tensorS+B tensorK.','scope':['All eight real support patterns have a positive portal-independent scalar Riccati barrier. Therefore no complete seven-quartic AF fixed ray exists in this two-active real family.','Complex flavor matrices, more active flavors, different representations, strong completion and asymptotic safety remain open. No universal SO10 or all-trajectory theorem is asserted.','The additional unequal/noncommuting fixed Yukawa branch is genuine; it does not cure the scalar obstruction.'],'prior_owners':['analysis/w33_pass11313_symmetric_adjoint_portal_closure.py'],'primary_sources':['https://arxiv.org/abs/hep-ph/0211440']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11320_general_real_two_flavor_yukawas.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],[(r['support'],r['minimum_numeric']) for r in d['branches']])
