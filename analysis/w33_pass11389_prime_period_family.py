#!/usr/bin/env python3
"""Independent small-prime witnesses for the selected-parabolic spatial cover.

The orbit argument proves fixed harmonic rankq for every finite field. This
producer checks the stronger integral period/root-lattice law only at prime
q=2,3,5; it does not implement extension fields or certify allq periods.
"""
from itertools import product,combinations
import numpy as np,sympy as s
from sympy.matrices.normalforms import hermite_normal_form

def prime_period_witness(q):
 assert s.isprime(q)
 def canon(v):
  for x in v:
   if x%q:
    k=pow(int(x)%q,-1,q);return tuple(int(y)*k%q for y in v)
 pts=sorted({canon(x)for x in product(range(q),repeat=4)if any(x)});idx={p:i for i,p in enumerate(pts)};n=len(pts)
 def sy(x,y):return (x[0]*y[2]-x[2]*y[0]+x[1]*y[3]-x[3]*y[1])%q
 lines=set()
 for i,j in combinations(range(n),2):
  if sy(pts[i],pts[j])==0:
   lines.add(tuple(sorted({idx[canon(tuple((a*x+b*y)%q for x,y in zip(pts[i],pts[j])))]for a,b in product(range(q),repeat=2)if a or b})))
 lines=sorted(lines);li={frozenset(l):i for i,l in enumerate(lines)}
 edges=[(p,n+l)for l,ps in enumerate(lines)for p in ps];ei={e:i for i,e in enumerate(edges)}
 gs=[]
 for a,b,c in product(range(q),repeat=3):
  m=np.eye(4,dtype=int);m[:2,2:]=[[a,b],[b,c]]
  pp=[idx[canon(m@np.array(p))]for p in pts];gs.append(pp+[n+li[frozenset(pp[p]for p in l)]for l in lines])
 seen=set();orbs=[]
 for i,(a,b)in enumerate(edges):
  if i not in seen:
   o=sorted({ei[(g[a],g[b])]for g in gs});orbs.append(o);seen.update(o)
 D=s.zeros(2*n,len(edges));O=s.zeros(len(edges),len(orbs))
 for i,(a,b)in enumerate(edges):D[a,i]=-1;D[b,i]=1
 for j,o in enumerate(orbs):
  for i in o:O[i,j]=1
 B=(D*O).nullspace();C=O*s.Matrix.hstack(*B)
 for j in range(C.cols):C[:,j]*=s.ilcm(*[x.q for x in C[:,j]])
 phi=s.zeros(2*n,q);adj=[[]for _ in range(2*n)]
 for i,(a,b)in enumerate(edges):adj[a].append((b,i,1));adj[b].append((a,i,-1))
 seen={0};queue=[0]
 for a in queue:
  for b,i,sg in adj[a]:
   if b not in seen:seen.add(b);queue.append(b);phi[b,:]=phi[a,:]+sg*C[i,:]
 P=C-D.T*phi;L=hermite_normal_form(P.T.applyfunc(s.Integer));S=q**3+q*q+q+1
 assert len(B)==q and len(orbs)==4*(q+1)
 assert C.T*C==q**3*S*(s.eye(q)+s.ones(q))
 assert L==hermite_normal_form(S*(s.eye(q)+s.ones(q)))
 U=(s.eye(q)+s.ones(q)).inv()*(L/S)
 assert all(x.q==1 for x in U) and abs(U.det())==1
 return {'q':q,'point_count':n,'flags':len(edges),'edge_orbits':len(orbs),'harmonic_rank':len(B),'gram_law':C.T*C==q**3*S*(s.eye(q)+s.ones(q)),'period_law':L==hermite_normal_form(S*(s.eye(q)+s.ones(q))),'root_lattice':'A'+str(q),'cartesian_quadratic_coefficient':str(s.Rational(q**3,2*S*S)),'raw_period_lattice':[[str(x)for x in L.row(i)]for i in range(q)]}
