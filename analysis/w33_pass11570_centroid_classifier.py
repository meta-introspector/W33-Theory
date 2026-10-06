import sys,itertools,numpy as np,sympy as sp
from math import lcm
sys.path.insert(0,'analysis')
import w33_pass10950_clock_albert_lorentz_spinor as P
J=P.build_clock_albert(); n=J['n']; prod=J['prod']; den=lcm(*[x.denominator for x in prod.flat])
Pi=np.array([[[int(x*den) for x in prod[i,j]] for j in range(n)] for i in range(n)],dtype=np.int64)
Ls=np.stack([Pi[i].T for i in range(n)])
comms=[Ls[i]@Ls[j]-Ls[j]@Ls[i] for i,j in itertools.combinations(range(n),2)]
basisD=[]; rr=[]
for M in comms:
 v=M.flatten().tolist()
 if P.rank(rr+[v])>len(rr):basisD.append(M);rr.append(v)
 if len(basisD)==52:break
H=[0,9,26,10,11,1,2,24,25]; C=[i for i in range(27) if i not in H]
cons=[[int(D[c,h]) for D in basisD] for h in H for c in C]
ns=P.nullspace(cons)
def comb(cv):
 q=lcm(*[x.denominator for x in cv]); return sum((int(x*q)*M for x,M in zip(cv,basisD)),np.zeros((27,27),dtype=np.int64))
Hstab=[comb(v) for v in ns]
cons2=cons+[[int(D[a,26]) for D in basisD] for a in range(27)]
SM=[comb(v) for v in P.nullspace(cons2)]
def basis_from_span(mats):
 out=[]; rr=[]
 for M in mats:
  v=M.flatten().tolist()
  if P.rank(rr+[v])>len(rr):out.append(M);rr.append(v)
 return out
def derived(base):
 return basis_from_span([A@B-B@A for i,A in enumerate(base) for B in base[i+1:]])
def coords_setup(base):
 F=sp.Matrix.hstack(*[sp.Matrix(M.reshape(-1,1).tolist()) for M in base])
 # choose pivot rows via pivots of transpose
 piv=F.T.rref()[1]
 assert len(piv)==len(base)
 Pm=F[list(piv),:]
 Pinv=Pm.inv()
 return F,list(piv),Pinv
def adjoints(base):
 F,piv,Pinv=coords_setup(base); ads=[]
 for A in base:
  cols=[]
  for B in base:
   V=sp.Matrix((A@B-B@A).reshape(-1,1).tolist())
   c=Pinv*V[piv,:]
   cols.append(c)
  ads.append(sp.Matrix.hstack(*cols))
 return ads
def centroid_info(base):
 ads=adjoints(base); d=len(base)
 vars=sp.symbols('x0:'+str(d*d)); X=sp.Matrix(d,d,vars)
 eq=[]
 for A in ads:
  M=X*A-A*X
  eq.extend(list(M))
 Mcoef=sp.linear_eq_to_matrix(eq,vars)[0]
 ns=Mcoef.nullspace()
 print('centroid dim',len(ns))
 mats=[sp.Matrix(d,d,v) for v in ns]
 Id=sp.eye(d)
 # pick non-scalar
 C=next(M for M in mats if not (M-M[0,0]*Id).is_zero_matrix)
 print('charpoly',sp.factor(C.charpoly().as_expr()))
 print('eigs',C.eigenvals())
 return len(ns),C.eigenvals()
print('Hstab'); centroid_info(Hstab)
Dsm=derived(SM)
print('Dsm'); centroid_info(Dsm)
