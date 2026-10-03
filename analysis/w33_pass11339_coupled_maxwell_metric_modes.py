"""Coupled Maxwell/metric torus-vector block of the actual warped wall saddle."""
from pathlib import Path
import json
import numpy as np
import sympy as s
from scipy.linalg import eigh
from numpy.polynomial.legendre import Legendre,leggauss
ROOT=Path(__file__).resolve().parents[1]

def spectrum(m,parity,n,quadrature=400):
 z,w=leggauss(quadrature);b=1+(z+1)/8;w/=8;F=8/(15*b)-b*b/30-1/(2*b*b);U=[];D=[]
 for k in range(n):
  P=Legendre.basis(k);a=m/2;x=b-1;u=x**a*P(z);d=x**a*P.deriv()(z)*8+(a*x**(a-1)*P(z) if a else 0)
  if parity=='odd':d=d*(1.25-b)-u;u=u*(1.25-b)
  U.append(u);D.append(d)
 U=np.array(U).T;D=np.array(D).T;M=U.T@(w[:,None]*U);K=D.T@((w*F)[:,None]*D)-2*U.T@((w/b**4)[:,None]*U)
 if m:K+=U.T@((w*m*m*.04/F)[:,None]*U)
 if m==0 and parity=='even':a=U.T@(w/b**4);K+=2*np.outer(a,a)/sum(w/b**4)
 return eigh(K,M,eigvals_only=True)

def ricci_control():
 b,c,e=s.symbols('b c e',positive=True);F=s.Function('F')(b);v=s.Function('v')(b);g=s.Matrix([[1/F,0,0,0],[0,F/c**2+b*b*e*e*v*v,b*b*e*v,0],[0,b*b*e*v,b*b,0],[0,0,0,b*b]]);gi=g.inv();D=lambda x,k:s.diff(x,b) if k==0 else 0
 G=[[[s.simplify(sum(gi[k,l]*(D(g[l,j],i)+D(g[l,i],j)-D(g[i,j],l))/2 for l in range(4))) for j in range(4)] for i in range(4)] for k in range(4)];R=0
 for i in range(4):
  for j in range(4):
   if gi[i,j]==0:continue
   rij=sum(D(G[k][i][j],k)-D(G[k][i][k],j)+sum(G[k][i][j]*G[l][k][l]-G[l][i][k]*G[k][j][l] for l in range(4)) for k in range(4));R+=gi[i,j]*rij
 R=s.simplify(R);difference=s.simplify(R-R.subs(e,0));assert s.simplify(difference+b*b*c*c*e*e*s.diff(v,b)**2/2)==0
 return str(difference)

def payload():
 b=s.symbols('b',positive=True);F=s.Rational(8,15)/b-b*b/30-1/(2*b*b);chi=1/b-s.Rational(4,5);v=s.Rational(5,2)/b**4-s.Rational(8,3)/b**3+s.Rational(1,6)
 assert s.simplify(-s.diff(F*s.diff(chi,b),b)-2*chi/b**4)==0 and chi.subs(b,s.Rational(5,4))==0
 assert v.subs(b,1)==0 and s.diff(v,b).subs(b,s.Rational(5,4))==0 and s.simplify(s.diff(v,b)/5+2*chi/b**4)==0
 supersol=s.factor(-s.diff(F,b,2)/2-s.diff(F,b)**2/(4*F)-2/b**4+s.Rational(1,25)/F)
 expected=(b-1)*(b+2)*(b**4-15*b**2-10*b-60)/(15*b**3*(b**3+b**2+b-15));assert s.simplify(supersol-expected)==0
 rows=[]
 for m in range(4):
  for parity in ['even','odd']:
   a=spectrum(m,parity,18)[:5];c=spectrum(m,parity,28)[:5];assert min(c)>-1e-8 and max(abs(a-c))<1e-5;rows.append({'m':m,'parity':parity,'eigenvalues':c.tolist(),'basis18_to28_error':float(max(abs(a-c)))})
 return {'status':'PASS','scope':'Coupled torus-independent Maxwell Wilson fields and mixed metric connections on the actual11324 Euclidean saddle; this full vector block is nonnegative with four exact zero modes, not total wall stability.','reduction':'ds²=gbase+b² sum_i(dy_i+Vi)^2, A=alpha+chi_i(dy_i+Vi). For eachi, H_i=dVi=h_i volbase. Quadratic density after torus integration: .5|dchi|²+(kappa² b⁴/4)h²+q chi h. Exact Ricci calculation and direct Maxwell contraction give these coefficients.','Ricci_connection_control':ricci_control(),'Maxwell_control':'In coordinate A=alpha+chi dy and Vphi=v, 1/4 sqrtg F² at second order gives measurebase[.5 F chi_b²-q c v chi_b]; integrating by parts with regular poles gives q chi h, h=c v_b. The wall has no Maxwell2-form source and induced volume is independent ofVi.','global_connection_constraint':'Vi is a globally defined connection on the trivial torus bundle: integral_S2 Hi=0. Elimination gives h=-2q(chi-mu)/(kappa² b⁴), mu=integral chi/b⁴ / integral1/b⁴. Omitting this global constraint would incorrectly give a negative mass to constant Wilson lines.','Schur_form':'Q=integralbase[|dchi|²-2q²(chi-mu)²/(kappa² b⁴)]. For m0 even, add the rank-one weighted-mean term; for odd or m!=0, mu=0. Constants remain Wilson-line zero modes.','exact_odd_zero':{'chi':str(chi),'Vphi':str(v),'Maxwell_nonzero_variation':'dchi/db=-1/b², so this is not a small Maxwell gauge transformation.','junctions':'chi odd with chiwall0; Vphi even with Vphi_prime_wall0. Both Maxwell flux continuity and mixed pure-tension Israel equation hold. Vphi vanishes quadratically in proper radius at both poles.','count':2},'Wilson_constant_zero_count':2,'vector_block_zero_count':4,'combined_shape_vector_zero_count':6,'positivity_proof':{'odd_m0':'chi0=1/b-4/5 is positive in the cap interior, regular at the pole and Dirichlet at the wall. Its exact zero equation and ground-state factorization give nonnegative oddm0 form.','even_m0':'Bare Neumann/Dirichlet Sturm interlacing puts exactly one negative eigenvalue in the bare even sector and its second eigenvalue strictly above the odd zero. The positive rank-one mean correction cannot lower the second level. Since its constant function is exactly zero, the lowest corrected level is zero and there is no negative level.','m_nonzero':'Positive supersolution sqrtF has L_m1 sqrtF/sqrtF equal to the stored rational expression. Both numerator polynomial b^4-15b²-10b-60 and denominator polynomial b³+b²+b-15 are negative throughout[1,5/4], so the ratio is nonnegative. Its wall derivative is positive; the factorization boundary term is nonnegative for even and zero for odd. Every|m|>1 adds a positive angular term.','m1_supersolution_ratio':str(supersol)},'spectral_controls':rows,'boundary':['Eigenvalues are Euclidean fluctuation values; integrating out the positive connection block determines negative-mode count by Schur complement, not the complete original coupled kinetic eigenvalues.','This vector block and prior11329 shape block leave the breathing/base/wall-bending/four-form block open.','The new odd zero is an infinitesimal mode; nonlinear integrability and quantum lifting are not proved. Wilson lines can be lifted by added charged matter, which is not supplied here.'],'prior_owners':['analysis/w33_pass11324_quantized_warped_history_wall.py','analysis/w33_pass11329_wall_shape_fluctuations.py'],'primary_sources':['https://arxiv.org/abs/0802.3564']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11339_coupled_maxwell_metric_modes.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['combined_shape_vector_zero_count'])
