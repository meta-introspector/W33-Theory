import sympy as s,itertools
from fractions import Fraction
pairs=[(0,1),(0,2),(1,2)]
sites=list(itertools.product(range(3),repeat=3)); idx={x:i for i,x in enumerate(sites)}
vals=[0,1,-1]
eps=s.Rational(1,10)
def frame(x):
 a,b,c=[vals[t] for t in x]
 return s.diag(1+eps*(a+b),1+eps*(b+c),1+eps*(c+a))
def shift(x,i,sgn):
 y=list(x);y[i]=(y[i]+sgn)%3;return tuple(y)
def delta_field(fun,x,i):
 return (fun(shift(x,i,1))-fun(shift(x,i,-1)))/2
def buildA(E):
 unk=[(i,a,b) for i in range(3) for a,b in pairs];rows=[]
 for a in range(3):
  for i,j in pairs:
   row=[]
   for ii,aa,bb in unk:
    v=0
    if ii==i:
     if a==aa:v+=E[bb,j]
     if a==bb:v-=E[aa,j]
    if ii==j:
     if a==aa:v-=E[bb,i]
     if a==bb:v+=E[aa,i]
    row.append(v)
   rows.append(row)
 return s.Matrix(rows),unk
omegas={}
torsmax=0
for x in sites:
 E=frame(x); A,unk=buildA(E); C=[]
 for a in range(3):
  for i,j in pairs:
   C.append(delta_field(lambda y:frame(y)[a,j],x,i)-delta_field(lambda y:frame(y)[a,i],x,j))
 w=A.LUsolve(-s.Matrix(C))
 Om=[s.zeros(3) for _ in range(3)]
 for coeff,(i,a,b) in zip(w,unk):Om[i][a,b]=coeff;Om[i][b,a]=-coeff
 omegas[x]=Om
 assert A*w+s.Matrix(C)==s.zeros(9,1)
# curvature
F={}
for x in sites:
 for i,j in pairs:
  D=(omegas[shift(x,i,1)][j]-omegas[shift(x,i,-1)][j])/2-(omegas[shift(x,j,1)][i]-omegas[shift(x,j,-1)][i])/2
  F[x,i,j]=s.simplify(D+omegas[x][i]*omegas[x][j]-omegas[x][j]*omegas[x][i])
# scalar contraction R=2 sum_{i<j,a,b} Einv[i,a] Einv[j,b] F_ij[a,b]
Rs=[]
for x in sites:
 Ei=frame(x).inv(); R=0
 for i,j in pairs:
  for a in range(3):
   for b in range(3):R+=2*Ei[i,a]*Ei[j,b]*F[x,i,j][a,b]
 Rs.append(s.factor(R))
print('nonzero F sites',sum(any(v!=0 for v in F[x,i,j]) for x in sites for i,j in pairs),'of',81)
print('R distinct',sorted(set(map(str,Rs)))[:20],'count',len(set(Rs)))
print('sumR',s.factor(sum(Rs)))
print('max numerator digits',max(len(str(abs(s.fraction(v)[0]))) for M in F.values() for v in M))
print('sampleR',[(sites[i],Rs[i]) for i in range(6)])
