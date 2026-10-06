import sys,itertools,numpy as np,sympy as sp
from math import lcm
from fractions import Fraction as Fr
sys.path.insert(0,'analysis')
import w33_pass10950_clock_albert_lorentz_spinor as P
J=P.build_clock_albert(); n=J['n']; prod=J['prod']; den=lcm(*[x.denominator for x in prod.flat])
Pi=np.array([[[int(x*den) for x in prod[i,j]] for j in range(n)] for i in range(n)],dtype=np.int64)
Ls=np.stack([Pi[i].T for i in range(n)])
I=np.eye(n,dtype=np.int64)
# exact F4 basis
comms=[Ls[i]@Ls[j]-Ls[j]@Ls[i] for i,j in itertools.combinations(range(n),2)]
basisD=[]; rows=[]
for M in comms:
 v=M.flatten().tolist()
 if P.rank(rows+[v])>len(rows):
  basisD.append(M);rows.append(v)
 if len(basisD)==52:break
print('F4',len(basisD))
H=[0,9,26,10,11,1,2,24,25]; Hset=set(H); C=[i for i in range(n) if i not in Hset]
# closure
for i in H:
 for j in H:
  nz=set(np.flatnonzero(Pi[i,j]))
  assert nz<=Hset,(i,j,nz)
print('H closure ok')
# constraints preserve H
cons=[]
for h in H:
 for c in C:
  cons.append([int(D[c,h]) for D in basisD])
ns=P.nullspace(cons)
print('H-stab dim',len(ns))
def comb(cv,base=basisD):
 q=lcm(*[x.denominator for x in cv]) if cv else 1
 return sum((int(x*q)*M for x,M in zip(cv,base)),np.zeros((27,27),dtype=np.int64))
Hstab=[comb(v) for v in ns]
print('Hstab rank',P.rank([M.flatten().tolist() for M in Hstab]))
# intersection fix primitive idempotent c=e26
cons2=list(cons)
for a in range(27): cons2.append([int(D[a,26]) for D in basisD])
ns2=P.nullspace(cons2)
print('intersection dim',len(ns2))
SM=[comb(v) for v in ns2]
# helper exact coordinates and Lie invariants
def basis_from_span(mats):
 out=[]; rr=[]
 for M in mats:
  v=M.flatten().tolist()
  if P.rank(rr+[v])>len(rr):
   out.append(M);rr.append(v)
 return out
def derived(base):
 br=[A@B-B@A for i,A in enumerate(base) for B in base[i+1:]]
 return basis_from_span(br)
def center_dim(base):
 cons=[]
 for B in base:
  # coeffs x_a, require [sum x_a A_a,B]=0
  for q in range(27):
   for r in range(27):
    cons.append([int((A@B-B@A)[q,r]) for A in base])
 return len(P.nullspace(cons))
for name,base in [('Hstab',Hstab),('SM',SM)]:
 der=derived(base)
 print(name,'dim',len(base),'center',center_dim(base),'derived',len(der))
# intersection action on H and complement ranks/traces
# generic trace form inertia
for name,base in [('Hstab',Hstab),('SM',SM)]:
 Gm=np.array([[np.trace(A@B) for B in base] for A in base],dtype=object)
 print(name,'trace inertia',P.inertia_int(Gm))
# derived SM center/trace
Dsm=derived(SM)
print('Dsm center',center_dim(Dsm),'dim',len(Dsm))
