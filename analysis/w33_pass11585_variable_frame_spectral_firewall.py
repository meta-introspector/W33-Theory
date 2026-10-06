import itertools,math,sys
import numpy as np
from scipy.linalg import expm
sys.path.insert(0,'analysis')
import w33_pass11557_pin_equivariant_event_dirac as P57
sites=list(itertools.product(range(3),repeat=3));idx={x:i for i,x in enumerate(sites)};pairs=[(0,1),(0,2),(1,2)]
vals=[0.,1.,-1.]
def shift(x,i,s):
 y=list(x);y[i]=(y[i]+s)%3;return tuple(y)
def raw_frame(x,eps):
 a,b,c=[vals[t] for t in x]
 return np.diag([1+eps*(a+b),1+eps*(b+c),1+eps*(c+a)])
def frames(eps):
 E={x:raw_frame(x,eps) for x in sites}
 vol=sum(np.linalg.det(M) for M in E.values());scale=(27/vol)**(1/3)
 return {x:scale*M for x,M in E.items()}
def delta(fun,x,i): return (fun(shift(x,i,1))-fun(shift(x,i,-1)))/2
def cartan_A(E):
 unk=[(i,a,b) for i in range(3) for a,b in pairs];rows=[]
 for a in range(3):
  for i,j in pairs:
   row=[]
   for ii,aa,bb in unk:
    v=0.
    if ii==i:
     if a==aa:v+=E[bb,j]
     if a==bb:v-=E[aa,j]
    if ii==j:
     if a==aa:v-=E[bb,i]
     if a==bb:v+=E[aa,i]
    row.append(v)
   rows.append(row)
 return np.array(rows),unk
def connection(E):
 out={}
 for x in sites:
  A,unk=cartan_A(E[x]);C=[]
  for a in range(3):
   for i,j in pairs:
    C.append(delta(lambda y:E[y][a,j],x,i)-delta(lambda y:E[y][a,i],x,j))
  w=np.linalg.solve(A,-np.array(C));Om=[np.zeros((3,3)) for _ in range(3)]
  for coeff,(i,a,b) in zip(w,unk):Om[i][a,b]=coeff;Om[i][b,a]=-coeff
  out[x]=Om
 return out
def curvature(E,om):
 F={};Rs=[]
 for x in sites:
  for i,j in pairs:
   D=(om[shift(x,i,1)][j]-om[shift(x,i,-1)][j])/2-(om[shift(x,j,1)][i]-om[shift(x,j,-1)][i])/2
   F[x,i,j]=D+om[x][i]@om[x][j]-om[x][j]@om[x][i]
 for x in sites:
  Ei=np.linalg.inv(E[x]);R=0.
  for i,j in pairs:
   for a in range(3):
    for b in range(3):R+=2*Ei[i,a]*Ei[j,b]*F[x,i,j][a,b]
  Rs.append(R)
 return F,np.array(Rs)
# shifts/differences on sites
T=[]
for i in range(3):
 M=np.zeros((27,27))
 for x in sites:M[idx[x],idx[shift(x,i,1)]]=1
 T.append(M)
D=[(M-M.T)/2 for M in T]
W=sum(np.eye(27)-(M+M.T)/2 for M in T)
gam,grade=P57.small_clifford();g=[np.array(x.evalf(),complex) for x in gam];gr=np.array(grade.evalf(),complex)
def operator(E,om,r=.25):
 Q=np.zeros((108,108),complex)
 for a in range(3):
  for i in range(3):
   valsE=np.array([E[x][a,i] for x in sites]);Em=np.diag(valsE)
   Q+=1j*np.kron(g[a],(Em@D[i]+D[i]@Em)/2)
 # local spin connection term
 # local spin-connection block
 for p,x in enumerate(sites):
  S=np.zeros((4,4),complex)
  for i in range(3):
   O=sum((.5*om[x][i][b,c]*(g[b]@g[c]) for b,c in pairs),np.zeros((4,4),complex))
   for a in range(3):S+=1j*g[a]*E[x][a,i]*O
  for s1 in range(4):
   for s2 in range(4):Q[s1*27+p,s2*27+p]+=S[s1,s2]
 Q=(Q+Q.conj().T)/2
 Q+=r*np.kron(gr,W)
 return (Q+Q.conj().T)/2
rows=[]
flat=None
for eps in [0.,.02,.04,.06,.08,.10]:
 E=frames(eps);om=connection(E);F,R=curvature(E,om);Q=operator(E,om)
 ev=np.linalg.eigvalsh(Q);vint=sum(np.linalg.det(E[x]) for x in sites);Sint=sum(np.linalg.det(E[x])*R[idx[x]] for x in sites)
 hs=[float(np.exp(-t*ev*ev).sum()) for t in [.02,.05,.1,.2]]
 if flat is None:flat=hs
 row=dict(eps=eps,volume=vint,intR=float(Sint),trQ2=float(np.sum(ev*ev)),heat=hs,delta=[h-f for h,f in zip(hs,flat)])
 rows.append(row);print(row,flush=True)
# ratios in weak field; check whether deltaH/intR stabilizes
for j,t in enumerate([.02,.05,.1,.2]):
 rr=[r['delta'][j]/(r['eps']**2) for r in rows[1:]]
 X=np.array([[r['eps']**2,r['eps']**4] for r in rows[1:]])
 y=np.array([r['delta'][j] for r in rows[1:]])
 co=np.linalg.lstsq(X,y,rcond=None)[0];res=np.linalg.norm(X@co-y)
 print('t',t,'delta/eps2',rr,'quad_coeff',co.tolist(),'fit_res',res,flush=True)
import runpy, sympy as sp
ex=runpy.run_path('analysis/w33_pass11571_discrete_curvature.py')
weighted=sp.factor(sum(ex['frame'](x).det()*ex['Rs'][ex['idx'][x]] for x in ex['sites']))
unweighted=sp.factor(sum(ex['Rs']))
print('exact_unweighted_sumR',unweighted,'exact_weighted_intR',weighted)
assert weighted==0 and str(unweighted)=='-360725/209088'
print('PASS')
