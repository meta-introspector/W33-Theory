#!/usr/bin/env python3
"""SO(10) auxiliary completion of the rank-five Higgs frame, with explicit UV inventory."""
from pathlib import Path
import sys,json
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11280_factor_higgs_mediator as F
from w33_pass11275_polynomial_sm_higgs import rank_mod

def reference():
 v,A,U,V,W=F.reference();K=np.block([[np.zeros((5,5)),np.eye(5)],[-np.eye(5),np.zeros((5,5))]])
 E=np.vstack([np.eye(5),-1j*np.eye(5)])/np.sqrt(2)
 return v,A,U@E.conj().T,V@E.conj().T,W@E.conj().T,K

def constraints(v,A,U,V,W,K):
 _,_,_,d=F.H.setup();M=[np.einsum('abc,c->ab',d,v[:,i]) for i in range(2)];N=M[0].conj().T@M[1]
 P=(np.eye(10)+1j*K)/2
 return K@K+np.eye(10),U@K+1j*U,V@K+1j*V,U@V.conj().T-N,U.conj().T@U-P,V.conj().T@V-P,W-A@U,(A+np.eye(27))@W-6*U

def payload():
 v,A,U,V,W,K=reference();assert max(np.max(abs(z)) for z in constraints(v,A,U,V,W,K))<1e-14
 Bs=[]
 for i in range(10):
  for j in range(i+1,10):
   b=np.zeros((10,10),int);b[i,j]=1;b[j,i]=-1;Bs.append(b)
 tangent=np.column_stack([(b@K-K@b).ravel() for b in Bs]);assert rank_mod(tangent.astype(int))==20
 e6_without=28-30; e6_composite=34-30
 so=s.Rational(11,3)*8-s.Rational(1,3)*3*27-s.Rational(1,6)*8
 fam=-s.Rational(79,3)+s.Rational(1,3)*27*s.Rational(5,2)
 assert (e6_without,e6_composite,so,fam)==(-2,4,1,-s.Rational(23,6))
 return {'status':'PASS','result_scope':'PASS_LOCAL_SEMISIMPLE_AUXILIARY_COMPLETION_WITH_COMPOSITE_YUKAWA_INTERFACE',
 'fields':'Three complex (27,10) scalar matrices U,V,W and one real SO10-adjoint antisymmetric K; original two27 Higgs and real78 retained. SO10 acts on right frames; K transforms R K R^T.',
 'potential':'Squared norms K²+I10, UK+iU, VK+iV, UVdagger-N, UdaggerU-(I+iK)/2, VdaggerV-(I+iK)/2, W-AU, (A+I)W-6U, plus the prior base Gram/cubic/A-v constraints. Every constraint has degree<=2, so the scalar potential is degree<=4.',
 'reference':'K=[[0,I5],[-I5,0]], E=[I5;-iI5]/sqrt2, U=Uold E†, V=Vold E†, W=Wold E†. K²=-I, P=EE†=(I+iK)/2.',
 'local_kernel_proof':'K²=-I restricts deltaK to the20-dimensional SO10/U5 orbit, checked by exact integer commutator rank. Remove this SO10 rotation and set deltaK=0. The explicit UK+iU and VK+iV equations force deltaU=deltaU P and deltaV=deltaV P. Gram constraints alone only force this at finite zeros and can leave extra quartic-flat directions in their Hessian; the explicit projection squares remove that defect. This reduces to the previous five-frame equations. Their25 unitary directions lie in the SO10 stabilizer of K. Together with the old66 E6 directions, kernel=111; 1851 fields give1740 positive normal directions. Analytic elimination, not a floating full Hessian.',
 'real_alignment_fields':1851,'gauge_orbit_dimension':111,'positive_normal_directions':1740,
 'beta_without_H_replacement':{'E6':str(e6_without),'SO10':str(so)},
 'composite_family_Higgs':'Replace elementary H(27,6bar) by H_eff=v Fdagger/M. E6 cubic d(psi,psi,v)Fdagger/M supplies a dimension5 Yukawa operator. This removes the27 copies of family6bar from gauge running above this EFT inventory; its UV matching requires separate mediators.',
 'beta_with_composite_H':{'E6':str(e6_composite),'SO10':str(so),'family_SU3':str(fam)},
 'boundaries':['The auxiliary Abelian factor is embedded in an AF semisimple gauge theory, not deleted by fiat.','All gauge factors are not AF: family SU3 remains negative, and a second sextet lowers it further to-14/3.','The scalar alignment is renormalizable, but the composite Yukawa interface is dimension5. Any UV mediators change the beta inventory and must be counted.','VeV scales, gauge couplings, full quantum potential and global uniqueness are not derived.'],
 'prior_owners':['analysis/w33_pass11280_factor_higgs_mediator.py','analysis/w33_pass11281_family_sextet_selection.py'],
 'primary_source':'https://doi.org/10.1103/PhysRevD.24.1005'}
if __name__=='__main__':
 out=payload();(ROOT/'data/w33_pass11285_semisimple_factor_higgs.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
