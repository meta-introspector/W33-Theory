"""Exact complex-isotropic Yukawa ray and full seven-quartic search."""
from pathlib import Path
import sys,json
import numpy as np
import sympy as s
from scipy.optimize import root
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11313_symmetric_adjoint_portal_closure as P
from w33_pass11325_complex_yukawa_phase_audit import beta

def flavor():
 x=s.Rational(36865,19203);w=s.Rational(147095,57609);t=s.Rational(218635,76812);k=s.Rational(73,3)
 assert s.Rational(41,5)*x+t+s.Rational(9,4)*w==k and s.Rational(27,10)*x+s.Rational(15,2)*w==k and 2*x+s.Rational(36,5)*t==k
 A=np.diag([0,np.sqrt(float(x)),np.sqrt(float(x)),np.sqrt(float(t))]).astype(complex);B=np.zeros((4,4),complex);B[0,1]=np.sqrt(float(w)/2);B[0,2]=1j*np.sqrt(float(w)/2);B-=B.T
 assert max(np.linalg.norm(v) for v in beta(A,B,float(k)))<1e-11
 return x,w,t,A,B

def quartic_data(A,B):
 d=json.loads((ROOT/'data/w33_pass11313_symmetric_adjoint_portal_closure.json').read_text());C=np.array(d['scalar_beta_tensor']);D=np.array(d['gauge_quartic_source']);rng=np.random.default_rng(113400);X=[];box=[]
 for _ in range(32):
  S,K=P.sample(rng);M=np.kron(A,S)+np.kron(B,K);W=M.conj().T@M;X.append(P.invariants(S,K));box.append(np.trace(W@W).real)
 co=np.linalg.lstsq(X,box,rcond=None)[0];err=float(max(abs(np.array(X)@co-box)));assert err<1e-9
 ga=np.trace(A.conj().T@A).real;gb=np.trace(B.conj().T@B).real;linear=np.array([120-4*ga]*2+[96-4*gb]*2+[108-2*ga-2*gb]*3)-16/3
 return C,D,co,linear,err

def stronger_barrier(x,t):
 d=json.loads((ROOT/'data/w33_pass11313_symmetric_adjoint_portal_closure.json').read_text());C=[s.Matrix([[s.Rational(str(v)) for v in row] for row in a]) for a in d['scalar_beta_tensor']];D=s.Matrix([s.Rational(str(v)) for v in d['gauge_quartic_source']]);Q=C[0]+C[1]/5;ix=[0,1,4,5,6];minors=[s.factor(Q.extract(ix[:i],ix[:i]).det()) for i in range(1,6)];assert all(v>0 for v in minors)
 v=s.Matrix([1,s.Rational(1,5)]);q=Q.extract([0,1],[0,1]);k=s.factor(1/(v.T*q.inv()*v)[0]);c=s.factor(D[0]+D[1]/5-s.Rational(4,5)*(2*x*x+t*t));w=s.Rational(344,3)-4*(2*x+t);mu=s.factor(c-w*w/(4*k));assert mu>0
 a2=(1-1/s.sqrt(10))/6;b2=(1+2/s.sqrt(10))/6;assert s.simplify(4*a2+2*b2-1)==0 and s.simplify(4*a2*a2+2*b2*b2-s.Rational(1,5))==0
 return {'background':'S=diag(a,-a,a,-a,b,-b,0,0,0,0), a²=(1-1/sqrt10)/6,b²=(1+2/sqrt10)/6; TrS²1,TrS4=1/5.','functional':'z=lambdaSdouble+lambdaSsingle/5','positive_principal_minors':[str(v) for v in minors],'Riccati_coefficient':str(k),'linear_coefficient':str(w),'constant':str(c),'exact_positive_global_minimum':str(mu),'minimum_numeric':float(mu),'theorem':'For this exact complex-isotropic Yukawa ray and unchanged8-Weyl gauge budget, dz/dtau>=124 z²-w z+c has a strictly positive minimum. Thus no simultaneous real7-quartic AF fixed ray exists on this named Yukawa branch. The old rS3/10 bound fails; the new rS1/5 physical background closes its apparent loophole.'}

def payload():
 x,w,t,A,B=flavor();C,D,box,linear,err=quartic_data(A,B);fun=lambda z:np.einsum('aij,i,j->a',C,z,z)+D-linear*z-box
 rng=np.random.default_rng(11340);sols=[]
 for _ in range(800):
  z=root(fun,rng.normal(size=7)*.5+np.array([.3,-.1,.2,.1,0,0,0]),tol=1e-10)
  if max(abs(fun(z.x)))<1e-7 and all(np.linalg.norm(z.x-v)>1e-5 for v in sols):sols.append(z.x)
 rows=[]
 for z in sols:
  bounds=P.bounded_control(z);h=1e-5;J=np.column_stack([(fun(z+h*e)-fun(z-h*e))/(2*h) for e in np.eye(7)]);rows.append({'quartic_ratios':z.tolist(),'residual':float(max(abs(fun(z)))),'sufficient_boundedness':list(map(float,bounds)),'sufficiently_bounded':bool(bounds[0]>0 and bounds[1]>0 and bounds[3]>0),'ratio_linearization_eigenvalues':np.linalg.eigvals(J).real.tolist()})
 ga=2*x+t;gb=2*w;ba=2*x*x+t*t;bb=2*w*w;sa=s.factor(36-s.Rational(6,5)*ba-(s.Rational(344,3)-4*ga)**2/(4*s.Rational(16580,159)));kb=s.factor(s.Rational(102,5)-s.Rational(4,5)*bb-(s.Rational(272,3)-4*gb)**2/(4*s.Rational(6571,62)));assert sa<0 and kb<0
 return {'status':'PASS','scope':'Exact allowed complex4-active Yukawa fixed ray, stronger physical-background scalarAF obstruction on this branch, and finite800-start full7-quartic control.','exact_ratios':{'a_squared':str(x),'b_squared':str(w),'t_squared':str(t)},'texture':'A=diag(0,a,a,t), B01=b/sqrt2, B02=i b/sqrt2, B10=-B01,B20=-B02. The coupled direction (1,i) is complex-isotropic, so BA†B=0 despite nonzeroA,B. Four active plus four spectator vector Weyls, unchanged total8.','beta_ratio_residual':float(max(np.linalg.norm(v) for v in beta(A,B,73/3))),'A_real':A.real.tolist(),'A_imag':A.imag.tolist(),'B_real':B.real.tolist(),'B_imag':B.imag.tolist(),'exact_S_lower_bound':str(sa),'exact_K_lower_bound':str(kb),'fermion_box_coefficients':box.tolist(),'box_projection_error':err,'quartic_roots':rows,'stronger_exact_barrier':stronger_barrier(x,t),'root_count':len(rows),'sufficiently_bounded_roots':sum(r['sufficiently_bounded'] for r in rows),'boundary':['A negative lower bound removes that previous no-go test; it does not prove the corresponding beta admits a bounded fixed point.','The stronger physical-background barrier proves absence of7-quartic AF fixed points on this exact branch independently of the finite root search. Other larger textures, all-trajectory running and full-inventory completion remain open.','The complex isotropic branch is outside the prior two-active and paired-isotropic restrictions; those certified barriers remain valid.'],'prior_owners':['analysis/w33_pass11335_larger_complex_flavor_search.py','analysis/w33_pass11313_symmetric_adjoint_portal_closure.py'],'primary_sources':['https://arxiv.org/abs/hep-ph/0211440']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11340_isotropic_flavor_portal.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['root_count'],d['sufficiently_bounded_roots'])
