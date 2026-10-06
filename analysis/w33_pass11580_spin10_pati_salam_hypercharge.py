import json,runpy,itertools
import sympy as sp, numpy as np
from fractions import Fraction as Fr
from math import lcm

ns=runpy.run_path('analysis/w33_pass11574_event_plane_from_sm.py')
P=ns['P']; SM=ns['SM']; ideals=ns['ideals']; ZM=sp.Matrix(ns['ZM'])
I=ns['I']; Ls=ns['Ls']; den=ns['den']
# exact Peirce16 basis
Lc=np.tensordot(I[26],Ls,axes=1)
Ahalf=P.nullspace((2*Lc-den*I).tolist())
Hi=[np.array([int(x*lcm(*[y.denominator for y in v])) for x in v],dtype=np.int64) for v in Ahalf]
H=sp.Matrix(np.stack(Hi).T.tolist()); piv=list(H.T.rref()[1]); Hleft=H[piv,:].inv()
def rep16(M):
    M=sp.Matrix(M.tolist()) if not isinstance(M,sp.MatrixBase) else M
    img=M*H
    return Hleft*img[piv,:]
su3_27=next(mats for mult,lam,mats in ideals if mult==8)
su3=[rep16(M) for M in su3_27]
B=rep16(ZM)/2 # eigenvalues i * q_BL, q_BL=+-1,+-3
print('B eig',B.eigenvals())

# exact gamma9 from committed certificate
gdata=json.load(open('data/w33_pass10961_albert_clifford9_gammas.json',encoding='utf-8'))
g9=[sp.Matrix([[sp.Rational(x) for x in row] for row in M]) for M in gdata['gamma9']]
Z16=sp.zeros(16); I16=sp.eye(16)
g10=[sp.Matrix.vstack(sp.Matrix.hstack(Z16,g),sp.Matrix.hstack(g,Z16)) for g in g9]
G10=sp.diag(I16,-I16); g10.append(G10)
biv=[g10[i]*g10[j] for i,j in itertools.combinations(range(10),2)]
assert len(biv)==45
# embed su3+B into 32 as same action on both Albert halves
emb=[sp.diag(R,R) for R in su3+[B]]
# verify every embedded generator actually belongs to the same spin(10) bivector algebra
BF=sp.Matrix.hstack(*[X.reshape(1024,1) for X in biv]); bp=list(BF.T.rref()[1]); Bleft=BF[bp,:].inv()
for C in emb:
    coeff=Bleft*C.reshape(1024,1)[bp,:]
    assert BF*coeff==C.reshape(1024,1)
print('embedded su3+B lies in spin10 bivectors')
# solve centralizer in spin10 bivectors
cols=[]
for X in biv:
    cols.append(sp.Matrix.vstack(*[(X*C-C*X).reshape(1024,1) for C in emb]))
A=sp.Matrix.hstack(*cols)
ker=A.nullspace()
print('centralizer dim',len(ker))
def combo(v):
    return sum((v[i]*biv[i] for i in range(45)),sp.zeros(32))
cent=[combo(v) for v in ker]
# closure & algebra centroid
F=sp.Matrix.hstack(*[M.reshape(1024,1) for M in cent]); rp=list(F.T.rref()[1]); left=F[rp,:].inv()
def coords(M):return left*M.reshape(1024,1)[rp,:]
ads=[]
for X in cent:
    ads.append(sp.Matrix.hstack(*[coords(X*Y-Y*X) for Y in cent]))
xx=sp.symbols('x0:'+str(len(cent)**2)); XX=sp.Matrix(len(cent),len(cent),xx); eq=[]
for a in ads:eq.extend(list(XX*a-a*XX))
Cns=sp.linear_eq_to_matrix(eq,xx)[0].nullspace()
print('centroid dim',len(Cns))
Cmats=[sp.Matrix(len(cent),len(cent),v) for v in Cns]
C0=next(M for M in Cmats if not (M-M[0,0]*sp.eye(len(cent))).is_zero_matrix)
print('centroid eig',C0.eigenvals())
ideals4=[]
for lam,mult in C0.eigenvals().items():
    mats=[]
    for v in (C0-lam*sp.eye(len(cent))).nullspace():
        mats.append(sum((v[i]*cent[i] for i in range(len(cent))),sp.zeros(32)))
    ideals4.append((int(mult),lam,mats))
print('ideal dims',[x[0] for x in ideals4])

# chirality of Cl10: i*product has square 1
prod=sp.eye(32)
for g in g10:prod=prod*g
Chi=sp.I*prod
assert Chi*Chi==sp.eye(32)
assert all(Chi*X-X*Chi==sp.zeros(32) for X in biv)
print('chi square',Chi*Chi==sp.eye(32),'chi eig',Chi.eigenvals())
W=(Chi-sp.eye(32)).nullspace(); Wm=sp.Matrix.hstack(*W)
wp=list(Wm.T.rref()[1]); Wleft=Wm[wp,:].inv()
def repW(M):
    img=M*Wm
    assert Wm*(Wleft*img[wp,:])==img
    return Wleft*img[wp,:]
Bw=repW(sp.diag(B,B))
print('W B eig',Bw.eigenvals())

# classify the two su2 ideals on Weyl16 via fixed spaces and B spectra
rows=[]
for ii,(mult,lam,mats) in enumerate(ideals4):
    if mult!=3:
        print('central ideal',ii,'dim',mult)
        continue
    reps=[repW(M) for M in mats]
    fixed=sp.Matrix.vstack(*reps).nullspace()
    Fb=sp.Matrix.hstack(*fixed)
    fp=list(Fb.T.rref()[1]); Fl=Fb[fp,:].inv()
    Bfix=Fl*(Bw*Fb)[fp,:]
    # active quotient inferred via whole minus fixed spectra; inspect a generic generator
    Hc=reps[0]
    eig=Hc.eigenvals()
    rows.append((ii,reps,fixed,Bfix.eigenvals(),eig))
    print('ideal',ii,'fixed',len(fixed),'Bfix',Bfix.eigenvals(),'H0eig',eig)

# pick SU2_R as ideal whose fixed sector contains qBL +1^6,-3^2 (left doublets)
ridx=next(ii for ii,reps,fixed,bfix,eig in rows if bfix.get(sp.I,0)==6 and bfix.get(-3*sp.I,0)==2)
three_ids=[ii for ii,reps,fixed,bfix,eig in rows]
lidx=next(ii for ii in three_ids if ii!=ridx)
print('L idx',lidx,'R idx',ridx)
Rreps=rows[ridx][1]
# choose a generator with two nonzero imaginary eigenvalues +/- i*c and zero^8
Hr=next(M for M in Rreps if len(M.eigenvals())>=3)
print('Hr eig',Hr.eigenvals())
# exact c^2 from -nonzero eigen^2
nonzero=[lam for lam in Hr.eigenvals() if lam!=0]
c2=sp.simplify(-(nonzero[0]**2)); c=sp.sqrt(c2)
R=sp.simplify(Hr/c)
print('R normalized eig',R.eigenvals())
print('comm BR',Bw*R-R*Bw==sp.zeros(16))
# simultaneous spectrum by R eigenspaces; report B then y=q+3r
spec=[]
for rlam,rmult in R.eigenvals().items():
    V=(R-rlam*sp.eye(16)).nullspace(); Vmat=sp.Matrix.hstack(*V); vp=list(Vmat.T.rref()[1]); Vl=Vmat[vp,:].inv()
    Br=Vl*(Bw*Vmat)[vp,:]
    re=sp.simplify(rlam/sp.I)
    for blam,bmult in Br.eigenvals().items():
        q=sp.simplify(blam/sp.I); y=sp.simplify(q+3*re)
        spec.append((re,q,y,int(bmult)))
print('joint r,q,y,mult',spec)
# expected one-family hypercharge integer spectrum
from collections import Counter
cnt=Counter()
for r,q,y,m in spec:cnt[str(y)]+=m
print('Y spectrum',dict(cnt))
assert cnt==Counter({'1':6,'-3':2,'2':3,'-4':3,'6':1,'0':1})
