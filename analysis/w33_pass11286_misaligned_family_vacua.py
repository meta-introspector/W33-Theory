#!/usr/bin/env python3
"""CP-even spectral decorrelation selects nonzero CKM in a two-sextet toy EFT."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))

def fourier():
 return np.exp(2j*np.pi*np.outer(np.arange(3),np.arange(3))/3)/np.sqrt(3)

def vacuum():
 V=np.exp(1j*np.pi/6)*fourier();return np.diag([1.,2.,3.]).astype(complex),V@np.diag([1.,3.,5.])@V.T

def constraints(Fu,Fd):
 out=[]
 for F,sv in [(Fu,[1,2,3]),(Fd,[1,3,5])]:
  S=F@F.conj().T
  out.extend([(np.trace(np.linalg.matrix_power(S,k))-sum(x**(2*k) for x in sv)).real for k in (1,2,3)])
  z=np.linalg.det(F)-np.prod(sv);out.extend([z.real,z.imag])
 A=Fu@Fu.conj().T;B=Fd@Fd.conj().T
 for p in (1,2):
  for q in (1,2):
   a=np.linalg.matrix_power(A,p);b=np.linalg.matrix_power(B,q);out.append((np.trace(a@b)-np.trace(a)*np.trace(b)/3).real)
 return np.array(out)

def potential(Fu,Fd):return float(constraints(Fu,Fd)@constraints(Fu,Fd))

def directions():
 out=[]
 for i in range(3):
  for j in range(i,3):
   a=np.zeros((3,3),complex);a[i,j]=a[j,i]=1;out.extend([a,1j*a])
 return out

def jarlskog(V):return float((V[0,0]*V[1,1]*V[0,1].conjugate()*V[1,0].conjugate()).imag)

def payload():
 Fu,Fd=vacuum();c=constraints(Fu,Fd);assert max(abs(c))<1e-8
 ds=directions();cols=[];h=1e-5
 for sector in range(2):
  for a in ds:
   plus=(Fu+h*a,Fd) if sector==0 else (Fu,Fd+h*a);minus=(Fu-h*a,Fd) if sector==0 else (Fu,Fd-h*a)
   cols.append((constraints(*plus)-constraints(*minus))/(2*h))
 J=np.array(cols).T;sv=np.linalg.svd(J,compute_uv=False);rank=int(sum(sv>1e-5));assert rank==12
 V=fourier();j=jarlskog(V);assert abs(j*j-1/108)<1e-15
 A=Fu@Fu.conj().T;B=Fd@Fd.conj().T;comm=A@B-B@A
 gap=lambda a:np.prod([a[i]-a[k] for i,k in [(0,1),(1,2),(2,0)]])
 inv=np.trace(comm@comm@comm).imag/(6*gap([1,4,9])*gap([1,9,25]));assert abs(abs(inv)-abs(j))<1e-14
 assert potential(Fu,Fd.conj())<1e-15 and potential(Fu,np.diag([1.,3.,5.]))>1
 return {'status':'PASS','result_scope':'PASS_CP_EVEN_TWO_SEXTET_TOY_VACUUM_WITH_CALCULABLE_NONZERO_MIXING',
 'potential':'Two prior spectral/determinant sextet sum-of-squares selectors, plus sum_p,q=1,2 [Tr(Su^p Sd^q)-Tr(Su^p)Tr(Sd^q)/3]^2. CP-even, SU3-congruence invariant; maximal degree16. Positive weights arbitrary.',
 'exact_mixing_argument':'For distinct eigenvalues, row/column normalization and the four mixed moments invert the two3x3 Vandermonde matrices, forcing |Vij|²=1/3. Every3x3 complex Hadamard has J²=1/108; equivalently its unitarity triangle is equilateral. Thus no CP-conserving mixing matrix is a zero. CP-conjugate branches have equal energy.',
 'reference_singular_values':{'up':[1,2,3],'down':[1,3,5]},'CKM_modulus_squared':[[1/3]*3]*3,'Jarlskog_squared_exact':'1/108','Jarlskog_reference':j,
 'standard_angles':'sin²theta13=1/3, sin²theta12=sin²theta23=1/2; cosdelta=0 in the standard CKM convention.',
 'finite_step_constraint_rank':rank,'finite_step_singular_values':sv.tolist(),'real_field_dimension':24,'normal_directions':12,'SU3_gauge_orbit_dimension':8,'additional_flat_phase_directions':4,
 'scope':'The exact equations fix the CKM sector and exclude J=0, but four relative sextet phases remain flat in addition to gauge directions. Numerical rank is a control, not an interval proof. Full scalar isolation and observed hierarchical CKM are not claimed.',
 'inputs':'Nondegenerate singular-value targets and the zero cross-correlation coupling pattern are inputs. Mixing is a consequence of that declared potential, not a fit to observed CKM or a consequence of W33 alone.',
 'prior_owners':['analysis/w33_pass11281_family_sextet_selection.py','analysis/BT891_yukawa_texture_from_grading.md','analysis/BT919_mixing_cp_scorecard.md'],
 'primary_sources':['https://arxiv.org/abs/1103.2915','https://journals.aps.org/prd/abstract/10.1103/PhysRevD.35.1685']}
if __name__=='__main__':
 out=payload();(ROOT/'data/w33_pass11286_misaligned_family_vacua.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
