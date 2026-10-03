"""Finite larger-complex Yukawa fixed-ray search; no completeness claim."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11325_complex_yukawa_phase_audit import beta

def bases(n):
 S=[];K=[]
 for i in range(n):
  for j in range(i,n):
   a=np.zeros((n,n));a[i,j]=a[j,i]=1 if i==j else 1/np.sqrt(2);S.append(a)
   if i!=j:
    a=np.zeros((n,n));a[i,j]=1/np.sqrt(2);a[j,i]=-1/np.sqrt(2);K.append(a)
 return np.array(S),np.array(K)
def unpack(x,S,K):
 ns=len(S);nk=len(K);z=x[:ns+nk]+1j*x[ns+nk:];return np.einsum('a,aij->ij',z[:ns],S),np.einsum('a,aij->ij',z[ns:],K)
def pack(A,B,S,K):
 z=np.r_[np.einsum('aij,ij->a',S,A),np.einsum('aij,ij->a',K,B)];return np.r_[z.real,z.imag]
def barrier(A):
 H=A.conj().T@A;gamma=np.trace(H).real;box=np.trace(H@H).real;return float(36-1.2*box-(120-4*gamma-16/3)**2/(4*(16580/159)))
def payload():
 rng=np.random.default_rng(11335);rows=[]
 for n,starts in [(3,10),(4,10),(6,4),(8,4)]:
  S,K=bases(n)
  for ix in range(starts):
   x=rng.normal(size=2*n*n)/np.sqrt(n)*1.2
   def fun(x):
    A,B=unpack(x,S,K);fa,fb=beta(A,B,73/3);return pack(fa,fb,S,K)
   z=least_squares(fun,x,gtol=2e-9,ftol=2e-10,xtol=2e-10,max_nfev=180);A,B=unpack(z.x,S,K);res=float(np.linalg.norm(fun(z.x)))
   if res<1e-6 and np.linalg.norm(z.x)>1e-4:
    sv=np.linalg.svd(A,compute_uv=False);bw=np.linalg.svd(B,compute_uv=False);rows.append({'n':n,'start':ix,'residual':res,'A_singular_values':sv.tolist(),'B_singular_values':bw.tolist(),'S_portal_barrier':barrier(A),'K_portal_barrier':float(102/5-.8*np.sum(bw**4)-(96-4*np.sum(bw**2)-16/3)**2/(4*(6571/62))),'coordinates':z.x.tolist()})
 assert rows
 return {'status':'PASS','scope':'Finite seeded stationary fixed-matrix search, not exhaustive classification or exclusion of running-frame relative rays. Every accepted witness replays betaA=betaB=0.','rows':rows,'minimum_sampled_barrier':min(r['S_portal_barrier'] for r in rows),'both_barriers_nonpositive_count':sum(r['S_portal_barrier']<=0 and r['K_portal_barrier']<=0 for r in rows),'search_starts':{'3':10,'4':10,'6':4,'8':4},'additional_bound':'ReTr(A† betaA)+k gamma >= gamma²+(31/5)Tr(A†A)²+(7/2)Tr(A†A B†B). This follows ReTr(A†B)² >= -Tr(A†A B†B). The relaxed bound alone admits negative portal lower bounds and is not a no-go theorem.','prior_owners':['analysis/w33_pass11325_complex_yukawa_phase_audit.py','analysis/w33_pass11313_symmetric_adjoint_portal_closure.py'],'primary_sources':['https://arxiv.org/abs/hep-ph/0211440']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11335_larger_complex_flavor_search.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],len(d['rows']),d['minimum_sampled_barrier'])
