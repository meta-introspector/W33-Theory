import sys,itertools,numpy as np,sympy as sp
from math import lcm
from fractions import Fraction as Fr
sys.path.insert(0,'analysis')
import w33_pass10950_clock_albert_lorentz_spinor as P
J=P.build_clock_albert(); n=J['n']; prodF=J['prod']; den=lcm(*[x.denominator for x in prodF.flat])
Pi=np.array([[[int(x*den) for x in prodF[i,j]] for j in range(n)] for i in range(n)],dtype=np.int64)
Ls=np.stack([Pi[i].T for i in range(n)])
I=np.eye(n,dtype=np.int64); idx,G=J['idx'],J['G']; evec=sum(I[idx[g]] for g in G); cvec=I[26]
def Lv(x):return np.tensordot(np.asarray(x,dtype=np.int64),Ls,axes=1)
def mulv(x,y):return np.asarray(x,dtype=np.int64)@np.tensordot(Pi,np.asarray(y,dtype=np.int64),axes=([1],[0]))
trL=np.array([np.trace(Ls[i]) for i in range(n)])
def trform(x,y):return Fr(int(np.dot(trL,mulv(x,y))),9*den*den)
# F4
comms=[Ls[i]@Ls[j]-Ls[j]@Ls[i] for i,j in itertools.combinations(range(n),2)]
basisD=[];rr=[]
for M in comms:
 v=M.flatten().tolist()
 if P.rank(rr+[v])>len(rr):basisD.append(M);rr.append(v)
 if len(basisD)==52:break
H=[0,9,26,10,11,1,2,24,25]; C=[i for i in range(27) if i not in H]
cons=[[int(D[c,h]) for D in basisD] for h in H for c in C]
cons2=cons+[[int(D[a,26]) for D in basisD] for a in range(27)]
def comb(cv,base=basisD):
 q=lcm(*[x.denominator for x in cv]);return sum((int(x*q)*M for x,M in zip(cv,base)),np.zeros((27,27),dtype=np.int64))
SM=[comb(v) for v in P.nullspace(cons2)]
def span_basis(mats):
 out=[];r=[]
 for M in mats:
  v=M.flatten().tolist()
  if P.rank(r+[v])>len(r):out.append(M);r.append(v)
 return out
Dsm=span_basis([A@B-B@A for i,A in enumerate(SM) for B in SM[i+1:]])
# coordinates on algebra
def coords_setup(base):
 F=sp.Matrix.hstack(*[sp.Matrix(M.reshape(-1,1).tolist()) for M in base])
 piv=list(F.T.rref()[1]); Pinv=F[piv,:].inv(); return F,piv,Pinv
def adjoints(base):
 F,piv,Pinv=coords_setup(base);ads=[]
 for A in base:
  cols=[]
  for B in base:
   V=sp.Matrix((A@B-B@A).reshape(-1,1).tolist());cols.append(Pinv*V[piv,:])
  ads.append(sp.Matrix.hstack(*cols))
 return ads
ads=adjoints(Dsm);d0=len(Dsm); xs=sp.symbols('x0:'+str(d0*d0));X=sp.Matrix(d0,d0,xs);eq=[]
for A in ads:eq.extend(list(X*A-A*X))
Mcoef=sp.linear_eq_to_matrix(eq,xs)[0];cent=[sp.Matrix(d0,d0,v) for v in Mcoef.nullspace()]
Cen=next(M for M in cent if not (M-M[0,0]*sp.eye(d0)).is_zero_matrix)
print('cent eig',Cen.eigenvals())
# ideal bases from eigenspaces; columns coeffs in Dsm
ideals=[]
for lam,mult in Cen.eigenvals().items():
 vs=(Cen-lam*sp.eye(d0)).nullspace()
 mats=[]
 for v in vs:
  q=lcm(*[Fr(sp.Rational(x).p,sp.Rational(x).q).denominator for x in v])
  mats.append(sum((int(sp.Rational(x)*q)*M for x,M in zip(v,Dsm)),np.zeros((27,27),dtype=np.int64)))
 ideals.append((int(mult),lam,mats))
# spatial 9 basis from pass10950
Lc=Lv(cvec)
A0=P.nullspace(Lc.tolist())
A0i=[np.array([int(x*lcm(*[y.denominator for y in v])) for x in v],dtype=np.int64) for v in A0]
uvec=evec-cvec;uu=trform(uvec,uvec);spatial=[]
for v in A0i:
 coef=trform(v,uvec)/uu;w=[Fr(int(a))-coef*int(b) for a,b in zip(v,uvec)];q=lcm(*[x.denominator for x in w]);spatial.append(np.array([int(x*q) for x in w],dtype=np.int64))
sp9=[]
for s in spatial:
 if P.rank([x.tolist() for x in sp9+[s]])>len(sp9):sp9.append(s)
S=sp.Matrix(np.stack(sp9).T.tolist());Spinv=(S.T*S).inv()*S.T
print('sp9',len(sp9))
def rep9(M):
 R=Spinv*sp.Matrix((M@np.stack(sp9).T).tolist()); assert S*R==sp.Matrix((M@np.stack(sp9).T).tolist());return R
for mult,lam,mats in ideals:
 reps=[rep9(M) for M in mats]
 eq=sp.Matrix.vstack(*reps)
 fix=eq.nullspace()
 print('ideal',mult,'lam',lam,'fixed dim',len(fix))
 # common image rank and traces
 print(' ranks reps span',sp.Matrix.hstack(*[R.reshape(81,1) for R in reps]).rank())
 if fix: print('fix vecs',[list(v.T) for v in fix])
# center generator of SM
# solve center coefficients
zs=sp.symbols('z0:'+str(len(SM)));Z=sum((zs[i]*sp.Matrix(SM[i].tolist()) for i in range(len(SM))),sp.zeros(27))
eqz=[]
for B in SM:eqz.extend(list(Z*sp.Matrix(B.tolist())-sp.Matrix(B.tolist())*Z))
Az=sp.linear_eq_to_matrix(eqz,zs)[0];zns=Az.nullspace();print('SM center',len(zns))
ZM=sum((sp.Rational(x)*sp.Matrix(M.tolist()) for x,M in zip(zns[0],SM)),sp.zeros(27))
Rz=rep9(np.array(ZM.tolist(),dtype=object))
print('center rep eig',Rz.eigenvals())
