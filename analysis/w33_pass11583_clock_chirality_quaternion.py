import itertools,json,sys
import numpy as np,sympy as sp
from fractions import Fraction as Fr
from math import lcm
from scipy.linalg import expm
sys.path.insert(0,'analysis')
import w33_pass10950_clock_albert_lorentz_spinor as P
# exact Cl10 chirality
gd=json.load(open('data/w33_pass10961_albert_clifford9_gammas.json',encoding='utf-8'))
g9=[sp.Matrix([[sp.Rational(x) for x in row] for row in M]) for M in gd['gamma9']]
I16=sp.eye(16);Z16=sp.zeros(16)
g=[sp.Matrix.vstack(sp.Matrix.hstack(Z16,x),sp.Matrix.hstack(x,Z16)) for x in g9]+[sp.diag(I16,-I16)]
Prod=sp.eye(32)
for x in g: Prod=Prod*x
Chi=sp.I*Prod
assert Chi*Chi==sp.eye(32) and Chi.eigenvals()=={-1:16,1:16}
for i,j in itertools.combinations(range(10),2):
 assert g[i]*g[j]*Chi==Chi*g[i]*g[j]
print('Spin10 chirality canonical',Chi.eigenvals())

# reconstruct deterministic Pass10956 half-turn in the same Peirce16 basis
J=P.build_clock_albert();n=J['n'];prod=J['prod'];den=lcm(*[x.denominator for x in prod.flat])
Pi=np.array([[[int(x*den) for x in prod[i,j]] for j in range(n)] for i in range(n)],dtype=np.int64)
Ls=np.stack([Pi[i].T for i in range(n)])
def Lv(x): return np.tensordot(np.asarray(x,dtype=np.int64),Ls,axes=1)
def mulv(x,y): return np.asarray(x,dtype=np.int64)@np.tensordot(Pi,np.asarray(y,dtype=np.int64),axes=([1],[0]))
I27=np.eye(n,dtype=np.int64);idx,G=J['idx'],J['G'];evec=sum(I27[idx[z]] for z in G);cvec=I27[idx[G[0]]];e2,e3=I27[idx[G[1]]],I27[idx[G[2]]];s0=e2-e3
trL=np.array([np.trace(Ls[i]) for i in range(n)])
Lc=Lv(cvec)
Ahalf=P.nullspace((2*Lc-den*I27).tolist())
Hi=[np.array([int(x*lcm(*[y.denominator for y in v])) for x in v],dtype=np.int64) for v in Ahalf]
Harr=np.stack(Hi).T;Hmat=sp.Matrix(Harr.tolist());Hpinv=(Hmat.T*Hmat).inv()*Hmat.T
def rep16(M):
 img=sp.Matrix((M@Harr).tolist());R=Hpinv*img;assert Hmat*R==img;return R
A0=P.nullspace(Lc.tolist())
A0i=[np.array([int(x*lcm(*[y.denominator for y in v])) for x in v],dtype=np.int64) for v in A0]
uvec=evec-cvec
def trform(x,y): return Fr(int(np.dot(trL,mulv(x,y))),9*den*den)
uu=trform(uvec,uvec);spatial=[]
for v in A0i:
 coef=trform(v,uvec)/uu;w=[Fr(int(a))-coef*int(b) for a,b in zip(v,uvec)];q=lcm(*[x.denominator for x in w]);spatial.append(np.array([int(x*q) for x in w],dtype=np.int64))
sp9=[]
for s in spatial:
 if P.rank([x.tolist() for x in sp9+[s]])>len(sp9):sp9.append(s)
ss=trform(s0,s0);z=None
for v in sp9:
 coef=trform(v,s0)/ss;w=[Fr(int(a))-coef*int(b) for a,b in zip(v,s0)];q=lcm(*[x.denominator for x in w]);vv=np.array([int(x*q) for x in w],dtype=np.int64)
 if any(vv) and trform(vv,vv)!=0:z=vv;break
Drot=Lv(s0)@Lv(z)-Lv(z)@Lv(s0);Rrot=rep16(Drot)
A0mat=np.stack(A0i).T.astype(float)
R0=np.linalg.lstsq(A0mat,Drot.astype(float)@A0mat,rcond=None)[0]
freq=sorted(abs(x.imag) for x in np.linalg.eigvals(R0) if abs(x.imag)>1e-9)[0]
U=expm(np.pi/freq*np.array(Rrot.tolist(),float))
print('U2+I',np.linalg.norm(U@U+np.eye(16)))
# doubled clock as Pass10961
zeta=np.exp(1j*np.pi/4);GA=np.block([[zeta**-1*U,np.zeros((16,16))],[np.zeros((16,16)),zeta*U]])
Ci=np.array(Chi.evalf(),complex);Gi=np.linalg.inv(GA)
C1=GA@Ci@Gi;C2=np.linalg.matrix_power(GA,2)@Ci@np.linalg.matrix_power(Gi,2)
print('G Chi plus/minus errors',np.linalg.norm(C1-Ci),np.linalg.norm(C1+Ci))
print('G2 Chi plus/minus errors',np.linalg.norm(C2-Ci),np.linalg.norm(C2+Ci))
print('C1 involution',np.linalg.norm(C1@C1-np.eye(32)),'anticomm C1/Chi',np.linalg.norm(C1@Ci+Ci@C1))
K=Ci@C1
print('K square +I',np.linalg.norm(K@K+np.eye(32)))
print('G4 Chi + error',np.linalg.norm(np.linalg.matrix_power(GA,4)@Ci@np.linalg.matrix_power(Gi,4)-Ci))
# projectors and clock image overlap
Pp=(np.eye(32)+Ci)/2;Pm=(np.eye(32)-Ci)/2
print('G maps P+ to P-?',np.linalg.norm(GA@Pp@Gi-Pm))
print('G2 maps P+ to P-?',np.linalg.norm(np.linalg.matrix_power(GA,2)@Pp@np.linalg.matrix_power(Gi,2)-Pm))
print('PASS')
