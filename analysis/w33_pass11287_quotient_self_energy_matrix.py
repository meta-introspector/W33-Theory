#!/usr/bin/env python3
"""Canonical Kähler-quotient jets and full22-field nonanalytic one-loop self-energy.

Bounded polynomial jets use NumPy/SciPy/SymPy. No vector/UV
threshold or full renormalized pole spectrum is claimed by the bubble witness.
"""
from pathlib import Path
import sys,json,itertools,os
import numpy as np
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11270_full_condensate_moduli as Q

def geometry():
 p,N=Q.reference();H=Q.L.generators()[:78];hp=np.einsum('aij,j->ai',H,p)
 R=(hp.conj()@hp.T).real;Ri=np.linalg.pinv(R,rcond=1e-10)
 A=np.einsum('ia,kij,jb->kab',N.conj(),H,N)
 assert np.max(abs(np.einsum('aij,j->ai',H,p).conj()@N))<1e-12
 assert np.linalg.matrix_rank(R,tol=1e-10)==70
 return p,N,H,Ri,A

class Jet:
 """Degree-four holomorphic polynomial jets in eleven coordinates."""
 __array_priority__=1000
 dim=11;order=4
 @classmethod
 def setup(cls):
  terms=[]
  for d in range(5):
   for ix in itertools.combinations_with_replacement(range(11),d):
    a=tuple(ix.count(i) for i in range(11));terms.append(a)
  cls.terms=terms;cls.index={a:i for i,a in enumerate(terms)};cls.degree=np.array([sum(a) for a in terms]);cls.pairs={}
  for dx in range(5):
   for dy in range(5):
    aa=[];bb=[];cc=[]
    for i in np.where(cls.degree<=dx)[0]:
     for k in np.where(cls.degree<=min(dy,4-cls.degree[i]))[0]:
      a=tuple(x+y for x,y in zip(terms[i],terms[k]));aa.append(i);bb.append(k);cc.append(cls.index[a])
    cls.pairs[dx,dy]=tuple(np.array(t,int) for t in (aa,bb,cc))
 def __init__(self,v=0,deg=0):
  self.c=np.zeros(len(self.terms),complex);self.c[0]=v;self.deg=deg
 @classmethod
 def coeffs(cls,c,deg):
  z=cls();z.c=c;z.deg=deg;return z
 def __add__(self,v):
  if not isinstance(v,Jet):v=Jet(v)
  return Jet.coeffs(self.c+v.c,max(self.deg,v.deg))
 __radd__=__add__
 def __neg__(self):return Jet.coeffs(-self.c,self.deg)
 def __sub__(self,v):return self+(-v)
 def __rsub__(self,v):return -self+v
 def __mul__(self,v):
  if not isinstance(v,Jet):return Jet.coeffs(self.c*v,self.deg)
  a,b,c=self.pairs[self.deg,v.deg];w=self.c[a]*v.c[b]
  out=np.bincount(c,weights=w.real,minlength=len(self.c))+1j*np.bincount(c,weights=w.imag,minlength=len(self.c))
  return Jet.coeffs(out,min(4,self.deg+v.deg))
 __rmul__=__mul__
 def __truediv__(self,v):return self*(1/v)
 def __pow__(self,n):
  if int(n)==n and n>=0:
   out=Jet(1)
   for _ in range(int(n)):out=out*self
   return out
  c=self.c[0];u=self/c-1;out=Jet(1);power=Jet(1);coef=1
  for k in range(1,5):
   power=power*u;coef*= (n-k+1)/k;out+=coef*power
  return c**n*out

def holomorphic_W_jets(p,N):
 import sympy as s
 import math
 Jet.setup();V=[]
 for i in range(81):
  v=Jet(p[i],1)
  for a in range(11):
   ix=tuple(int(k==a) for k in range(11));v.c[Jet.index[ix]]=N[i,a]
  V.append(v)
 V=np.array(V,dtype=object).reshape(27,3);I=Q.D.I;tri,d,eps=I.tensors()
 X=[[[0 for _ in range(3)] for _ in range(3)] for _ in range(27)]
 for i,j,k,sg in tri:
  for a,b,c in itertools.permutations((i,j,k)):
   for u in range(3):
    for v in range(3):X[a][u][v]+=int(sg)*V[b,u]*V[c,v]
 print('polynomial quotient I6 jet',flush=True)
 perms=list(itertools.permutations(range(3)));raw6=0
 for i,j,k,sg in tri:
  for a,b,c in perms:
   for e,f,g in perms:raw6+=6*int(sg*eps[a,b,c]*eps[e,f,g])*X[i][a][e]*X[j][b][f]*X[k][c][g]
 u6=9*raw6/4
 keys,z,S,T,_=Q.D.cubic_circuit();Fs=[]
 for u,v,w in keys:
  val=0
  for i,j,k,sg in tri:
   for a,b,c in itertools.permutations((i,j,k)):val+=int(sg)*V[a,u]*V[b,v]*V[c,w]
  Fs.append(val)
 print('polynomial quotient I12/I18 jets',flush=True)
 Sraw,Traw=s.lambdify(z,[S,T],'numpy',cse=True)(*Fs)
 u12=(152*u6*u6-3645*Sraw)/32
 i18=(19683*Traw-48384/5*u6**3+13824/5*u6*u12)/1492992
 W=3/14*(-729*i18)**(-1/3)
 def tensor(n):
  out=np.zeros((11,)*n,complex)
  for ix in itertools.product(range(11),repeat=n):
   a=tuple(ix.count(i) for i in range(11));out[ix]=W.c[Jet.index[a]]*math.prod(math.factorial(t) for t in a)
  return out
 return W.c[0],*[tensor(n) for n in (1,2,3,4)],i18.c[0]

def jets():
 p,N,H,Ri,A=geometry();w0,w1,w2,w3,w4,i18=holomorphic_W_jets(p,N)
 B=np.column_stack([np.eye(11),1j*np.eye(11)])/np.sqrt(2)
 # Scalar Taylor coefficients: V=|W_i|_g²-18 ReW on the canonical quotient cone.
 f1=w2@B;f2=np.einsum('ijk,ja,kb->iab',w3,B,B)/2
 f3=np.einsum('ijkl,ja,kb,lc->iabc',w4,B,B,B)/6
 W1=w1@B;W2=np.einsum('ij,ia,jb->ab',w2,B,B)/2;W3=np.einsum('ijk,ia,jb,kc->abc',w3,B,B,B)/6
 grad=2*np.real(w1.conj()@f1)-18*W1.real
 mu=np.real(np.einsum('ia,kij,jb->kab',B.conj(),A,B));mu=(mu+mu.transpose(0,2,1))/2
 t=np.einsum('ab,bij->aij',Ri,mu)
 mh=np.einsum('aji,jk->aik',A,B.conj());ma=np.einsum('aij,jk->aik',A,B)
 row=np.einsum('i,aji->aj',w1.conj(),A)
 g2w=-np.einsum('aj,akl->jkl',row,t)-np.einsum('ab,i,aik,bjl->jkl',Ri,w1.conj(),mh,ma)
 corr2=-np.real(np.einsum('jkl,j->kl',g2w,w1))
 cross3=-2*np.real(np.einsum('jkl,jm->klm',g2w,f1))
 hp=np.einsum('aij,j->ai',H,p);hn=np.einsum('aij,jk->aik',H,N)
 L=np.einsum('ai,bij->abj',hp.conj(),hn);L=(L+L.transpose(1,0,2))/2
 R1=2*np.real(np.einsum('abi,ik->abk',L,B))
 tr=np.einsum('ab,i,bik->ak',Ri,w1.conj(),mh);ts=np.einsum('ab,j,bjk->ak',Ri,w1,ma)
 trs=np.einsum('ab,i,j,bji->a',Ri,w1.conj(),w1,A)
 lr=np.einsum('i,abi->ab',w1.conj(),L);ls=np.einsum('j,abj->ab',w1,L.conj())
 g3w=(np.einsum('a,abi,bjk->ijk',trs,R1,t)+np.einsum('ai,abj,bk->ijk',tr,R1,ts)
       +np.einsum('ai,ab,bjk->ijk',tr,ls,t)+np.einsum('ai,ab,bjk->ijk',ts,lr,t)).real
 hcoef=np.real(f1.conj().T@f1+2*np.einsum('i,iab->ab',w1.conj(),f2))-18*W2.real+corr2
 ccoef=2*np.real(np.einsum('ia,ibc->abc',f1.conj(),f2)+np.einsum('i,iabc->abc',w1.conj(),f3))-18*W3.real+cross3-g3w
 Hess=hcoef+hcoef.T;cubic=sum(ccoef.transpose(ix) for ix in itertools.permutations(range(3)))
 anti=np.einsum('ab,aki,blj,k->ijl',Ri,A,A,w1)+np.einsum('ab,akj,bli,k->ijl',Ri,A,A,w1)
 Y=np.concatenate([(w3+anti).transpose(2,0,1)/np.sqrt(2),1j*(w3-anti).transpose(2,0,1)/np.sqrt(2)],axis=0)
 geo={'K4':'-mu^T R0^+ mu/2','K5':'mu^T R0^+ R1 R0^+ mu/2','orbit_rank':70,'curvature_vertex_norm':float(np.linalg.norm(anti)),'I18_at_reference':float(i18.real),'method':'Bounded degree-four holomorphic polynomial jets, then analytic canonical Kähler-quotient contractions; no large compiled autodifferentiation graph.'}
 return Hess,cubic,w2,Y,grad,geo

def takagi(M):
 R=M.real;I=M.imag;w,U=np.linalg.eigh(np.block([[R,-I],[-I,-R]]));sel=np.where(w>1e-8)[0]
 V=U[:len(M),sel]+1j*U[len(M):,sel];m=w[sel]
 assert np.max(abs(V.conj().T@V-np.eye(len(M))))<1e-9
 assert np.max(abs(V.T@M@V-np.diag(m)))<1e-8
 return m,V

def cut_groups(H,C,M,Y,m=.01):
 lam,O=np.linalg.eigh(H);lam[abs(lam)<1e-8]=0.;assert min(lam)>-1e-8
 fm,V=takagi(M);g=m*m*np.einsum('ia,jb,kc,ijk->abc',O,O,O,C)
 y=m*np.einsum('ia,jb,kc,ijk->abc',O,V,V,Y)
 sm=m*m*lam;fm=m*fm;groups={}
 def add(kind,a,b,R,I=None):
  key=(kind,round(float(a),12),round(float(b),12))
  if key not in groups:groups[key]={'kind':kind,'a':float(a),'b':float(b),'R':np.zeros((22,22)),'I':np.zeros((22,22))}
  groups[key]['R']+=R
  if I is not None:groups[key]['I']+=I
 for b in range(22):
  for c in range(b,22):
   v=g[:,b,c];add('scalar',sm[b],sm[c],(1 if b==c else 2)*np.outer(v,v)/(32*np.pi))
 for b in range(11):
  for c in range(b,11):
   v=y[:,b,c];n=1 if b==c else 2;add('fermion',fm[b]**2,fm[c]**2,n*np.outer(v.real,v.real)/(16*np.pi),n*np.outer(v.imag,v.imag)/(16*np.pi))
 return sm,fm,O,list(groups.values())

def rho(group,t):
 a,b=group['a'],group['b'];th=(np.sqrt(a)+np.sqrt(b))**2
 if t<=th:return np.zeros((22,22))
 beta=np.sqrt(max(0,1-2*(a+b)/t+(a-b)**2/t**2))
 if group['kind']=='scalar':return beta*group['R']
 return beta*((t-th)*group['R']+(t-(np.sqrt(a)-np.sqrt(b))**2)*group['I'])

def dispersion(group,s,s0):
 a,b=group['a'],group['b'];th=(np.sqrt(a)+np.sqrt(b))**2;up=max(100*s,2*th+abs(s0))
 def beta(t):return 1. if a==b==0 else np.sqrt(max(0,1-2*(a+b)/t+(a-b)**2/t**2))
 def integ(poly):
  f=lambda t:beta(t)*poly(t)/(t-s0)**2
  v=quad(f,th,up,weight='cauchy',wvar=s,epsabs=1e-10,limit=300)[0] if th<s<up else quad(lambda t:f(t)/(t-s),th,up,epsabs=1e-10,limit=300)[0]
  v+=quad(lambda t:f(t)/(t-s),up,np.inf,epsabs=1e-10,limit=300)[0]
  return (s-s0)**2*v/np.pi
 if group['kind']=='scalar':return integ(lambda t:1)*group['R']
 return integ(lambda t:t-th)*group['R']+integ(lambda t:t-(np.sqrt(a)-np.sqrt(b))**2)*group['I']

def finite_slice_control(C):
 from scipy.linalg import expm
 p,N,H,*_=geometry();v=np.random.default_rng(8728).normal(size=22);v/=np.linalg.norm(v);z=(v[:11]+1j*v[11:])/np.sqrt(2)
 def value(h):
  q=p+h*N@z
  for _ in range(4):
   act=np.einsum('aij,j->ai',H,q);mu=np.real(np.einsum('i,ai->a',q.conj(),act));R=(act.conj()@act.T).real
   if max(abs(mu))<1e-13:break
   t=-.5*np.linalg.pinv(R,rcond=1e-10)@mu;q=expm(np.einsum('a,aij->ij',t,H))@q
  return Q.potential(q),float(max(abs(mu)))
 expected=float(np.einsum('ijk,i,j,k',C,v,v,v));rows=[]
 for h in (.004,.002,.001):
  values=[value(t*h) for t in (2,1,-1,-2)];got=(values[0][0]-2*values[1][0]+2*values[2][0]-values[3][0])/(2*h**3)
  rows.append({'step':h,'third_derivative':float(got),'error':float(got-expected),'moment_residual':max(x[1] for x in values)})
 extrapolated=(4*rows[2]['third_derivative']-rows[1]['third_derivative'])/3
 assert abs(extrapolated-expected)<2e-5 and abs(rows[2]['error'])<abs(rows[1]['error'])/3
 return {'seed':8728,'expected':expected,'rows':rows,'Richardson_error':float(extrapolated-expected),'method':'Original 81-field potential evaluated after complexified-gauge Newton projection onto the actual D-flat sheet; independent of polynomial jet contraction.'}

def payload():
 H,C,M,Y,grad,geo=jets();expected=np.array([0.]*8+[162/49]*10+[486/49]*2+[486/7,648/7])
 ev=np.linalg.eigvalsh(H);err=float(max(abs(ev-expected)));print('scalar profile error',err,flush=True);assert err<1e-7 and max(abs(grad))<1e-8
 sm,fm,O,groups=cut_groups(H,C,M,Y);s=sm[-1];s0=-.0001
 spectral=sum((rho(c,s) for c in groups),start=np.zeros((22,22)));real=sum((dispersion(c,s,s0) for c in groups),start=np.zeros((22,22)))
 assert min(np.linalg.eigvalsh(spectral))>-1e-11
 gold=sum((rho(c,s) for c in groups if c['kind']=='scalar' and max(c['a'],c['b'])<1e-10),start=np.zeros((22,22)))
 goldtarget=s*s/(8*np.pi);gold_err=abs(gold[-1,-1]-goldtarget);assert gold_err<1e-10
 return {'status':'PASS','result_scope':'PASS_FULL_22_FIELD_NONANALYTIC_SELF_ENERGY_IN_CANONICAL_KAHLER_QUOTIENT',
 'independent_D_flat_cubic_control':finite_slice_control(C),'geometry':geo,'gradient_max':float(max(abs(grad))),'scalar_profile_error':err,'Hessian':H.tolist(),'cubic_tensor':C.tolist(),'Weyl_mass_real':M.real.tolist(),'Weyl_mass_imag':M.imag.tolist(),'Yukawa_real':Y.real.tolist(),'Yukawa_imag':Y.imag.tolist(),
 'scalar_mass_squared':sm.tolist(),'fermion_masses':fm.tolist(),'probe_momentum_squared':float(s),'subtraction_point':s0,
 'self_energy_real':real.tolist(),'self_energy_absorptive':spectral.tolist(),'absorptive_eigenvalues':np.linalg.eigvalsh(spectral).tolist(),'radial_width_over_mass':float(spectral[-1,-1]/s),'Goldstone_cut_Ward_error':float(gold_err),
 'matrix_method':'All22 external fields,22 internal scalars and11 internal Weyl modes. Full covariant cubic potential and curvature-corrected Weyl vertices. Twice-subtracted spectral dispersion; matrix absorptive part is PSD by channel outer products.',
 'quotient_jet_scope':'Canonical Kähler quotient through degree5, obtained by minimizing complexified gauge norm around the actual70 orbit. This suffices for potential cubic and Weyl linear vertices at the regular vacuum; arbitrary noncanonical Kähler terms are not fixed.',
 'counterterm_boundary':'Four-scalar/curvature tadpoles and polynomial local terms are absorbed into declared matching conditions. Goldstone Ward mass/tadpole matching is required before extracting a complete physical pole spectrum; this reports the full nonanalytic bubble matrix at one common probe, not22 solved poles.',
 'boundaries':['Complex128 polynomial-jet coefficients and analytic contractions, not directed intervals.','Heavy vector multiplets and UV matching omitted; gauge coupling is not selected.','A full momentum-dependent matrix calculation is stronger than the earlier radial-only fixed-slice cuts, but does not determine observed masses.'],
 'prior_owners':['analysis/w33_pass11270_full_condensate_moduli.py','analysis/w33_pass11282_all_modulus_radial_self_energy.py'],
 'primary_sources':['https://arxiv.org/abs/hep-th/0203081','https://arxiv.org/abs/hep-ph/0502168','https://arxiv.org/html/2505.07931v1']}
if __name__=='__main__':
 out=payload();(ROOT/'data/w33_pass11287_quotient_self_energy_matrix.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'],out['radial_width_over_mass'])
