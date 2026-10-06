import numpy as np,sympy as sp
from math import lcm
from fractions import Fraction as Fr
exec(open('analysis/w33_pass11570_executable_sm_intersection.py',encoding='utf-8').read())
# Peirce 16 for c=e26
Lc=np.tensordot(I[26],Ls,axes=1)
Ahalf=P.nullspace((2*Lc-den*I).tolist())
Hi=[np.array([int(x*lcm(*[y.denominator for y in v])) for x in v],dtype=np.int64) for v in Ahalf]
Hmat=sp.Matrix(np.stack(Hi).T.tolist()); Hpinv=(Hmat.T*Hmat).inv()*Hmat.T
def rep16(M):
 R=Hpinv*sp.Matrix((M@np.stack(Hi).T).tolist()); assert Hmat*R==sp.Matrix((M@np.stack(Hi).T).tolist()); return R
reps=[rep16(M) for M in SM]
# commutant
x=sp.symbols('x0:'+str(16*16));X=sp.Matrix(16,16,x);eq=[]
for R in reps:eq.extend(list(X*R-R*X))
A=sp.linear_eq_to_matrix(eq,x)[0]; ns=A.nullspace(); print('commutant dim',len(ns))
for k,v in enumerate(ns):
 M=sp.Matrix(16,16,v);print('comm basis',k,'trace',sp.simplify(M.trace()),'char',sp.factor(M.charpoly().as_expr()))
# center of SM algebra
z=sp.symbols('z0:'+str(len(SM)));Z=sum((z[i]*sp.Matrix(SM[i].tolist()) for i in range(len(SM))),sp.zeros(27));eqz=[]
for B in SM:eqz.extend(list(Z*sp.Matrix(B.tolist())-sp.Matrix(B.tolist())*Z))
zv=sp.linear_eq_to_matrix(eqz,z)[0].nullspace()[0]
Z27=sum((sp.Rational(a)*sp.Matrix(B.tolist()) for a,B in zip(zv,SM)),sp.zeros(27));Rz=Hpinv*(Z27*Hmat)
print('center eig16',Rz.eigenvals());print('trQ',sp.simplify(Rz.trace()),'trQ3',sp.simplify((Rz**3).trace()))
# perturbative cubic gauge anomaly tensor
mx=0;bad=0
for a in range(len(reps)):
 for b in range(a,len(reps)):
  for c in range(b,len(reps)):
   val=sp.simplify(sp.trace(reps[a]*(reps[b]*reps[c]+reps[c]*reps[b])))
   if val!=0: bad+=1; print('nonzero',a,b,c,val); raise SystemExit
print('cubic anomaly entries nonzero',bad)
