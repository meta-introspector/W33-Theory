#!/usr/bin/env python3
"""Pass 11067: an explicit nontrivial H27 2-cocycle twists split K81 into chamber U81."""
from __future__ import annotations
import argparse,itertools,json
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11067_explicit_h27_central_cocycle_weld.json"

def E(i,j):
    M=sp.zeros(4);M[i,j]=1;return M
I=sp.eye(4)
X0=E(0,1)-E(3,2);X1=E(1,3);X2=E(0,3)+E(1,2);X3=E(0,2)
def xr(X,t):return I+t*X
def umat(a,b,c,d):return xr(X0,a)*xr(X1,b)*xr(X2,c)*xr(X3,d)

def hmul(x,y):
    a,b,c=x;A,B,C=y
    return ((a+A)%3,(b+B)%3,(c+C-A*b)%3)

def kappa(x,y):
    a,b,c=x;A,B,C=y
    return (A*A*b-2*A*c)%3

def rank_mod(A,p=3):
    M=np.array(A,dtype=np.int64)%p;n,m=M.shape;r=0
    for col in range(m):
        z=next((i for i in range(r,n) if M[i,col]%p),None)
        if z is None:continue
        M[[r,z]]=M[[z,r]]
        iv=pow(int(M[r,col]),-1,p);M[r]=(M[r]*iv)%p
        for i in range(r+1,n):
            q=int(M[i,col]%p)
            if q:M[i]=(M[i]-q*M[r])%p
        r+=1
        if r==n:break
    return r

def payload():
    a,b,c,d,A,B,C,D=sp.symbols("a b c d A B C D")
    lhs=sp.expand(umat(a,b,c,d)*umat(A,B,C,D))
    rhs=sp.expand(umat(a+A,b+B,c+C-A*b,d+D+A*A*b-2*A*c))
    assert sp.simplify(lhs-rhs)==sp.zeros(4)

    H=list(itertools.product(range(3),repeat=3))
    hi={x:i for i,x in enumerate(H)}
    for x,y,z in itertools.product(H,repeat=3):
        left=(kappa(x,y)+kappa(hmul(x,y),z))%3
        right=(kappa(y,z)+kappa(x,hmul(y,z)))%3
        assert left==right

    M=[];rhsb=[]
    for x,y in itertools.product(H,repeat=2):
        row=[0]*27
        row[hi[x]]=(row[hi[x]]+1)%3
        row[hi[y]]=(row[hi[y]]+1)%3
        row[hi[hmul(x,y)]]=(row[hi[hmul(x,y)]]-1)%3
        M.append(row);rhsb.append(kappa(x,y))
    M=np.array(M,dtype=np.int64);rhsb=np.array(rhsb,dtype=np.int64).reshape(-1,1)
    r=rank_mod(M);ra=rank_mod(np.hstack([M,rhsb]))
    assert (r,ra)==(25,26)

    old=json.loads((ROOT/"data/w33_pass11065_split_nonsplit_order81_bridge.json").read_text())
    assert old["U81"]["quotient_by_center"]["is_H27"]

    return {
      "schema":"w33.pass11067.explicit-h27-central-cocycle-weld.v1",
      "status":"PASS_U81_IS_THE_EXPLICIT_NONTRIVIAL_CENTRAL_COCYCLE_TWIST_OF_SPLIT_K81_ON_THE_SAME_81_LABELS",
      "headline":"In C2 root coordinates u(a,b,c,d)=x0(a)x1(b)x2(c)x3(d), multiplication is u(a,b,c,d)u(A,B,C,D)=u(a+A,b+B,c+C-Ab,d+D+A^2 b-2Ac). The first three coordinates are exactly the frozen H27 law, while the fourth is a central C3 coordinate twisted by kappa((a,b,c),(A,B,C))=A^2 b-2Ac.",
      "exact_group_law":{
        "normal_form":"u(a,b,c,d)=x0(a)x1(b)x2(c)x3(d)",
        "product":"(a,b,c,d)*(A,B,C,D)=(a+A,b+B,c+C-Ab,d+D+A^2 b-2Ac) mod 3",
        "quotient_law":"(a,b,c)*(A,B,C)=(a+A,b+B,c+C-Ab) mod 3",
        "quotient":"H27",
        "central_cocycle":"kappa(x,y)=A^2 b-2 A c mod 3"
      },
      "cohomology_certificate":{
        "cocycle_identity_checked_triples":19683,
        "coboundary_linear_system_shape":[729,27],
        "coefficient_rank_mod3":25,
        "augmented_rank_mod3":26,
        "is_coboundary":False
      },
      "same_carrier_bridge":{
        "K81_product":"(h,d)*(h',D)=(hh',d+D)",
        "U81_product":"(h,d)*(h',D)=(hh',d+D+kappa(h,h'))",
        "underlying_set":"H27 x F3, 81 labels"
      },
      "boundary":"The cocycle is an exact finite-group weld. Its polynomial degree does not by itself make it a physical interaction Hamiltonian or elapsed-time law.",
      "parents":["data/w33_pass11065_split_nonsplit_order81_bridge.json","data/PART_W33_PASS5120_U81_STATE_PROGRAM_TRANSPORT.json"],
      "checks":{"symbolic_matrix_product_identity":True,"cocycle_identity_exhaustive":True,"noncoboundary_rank_witness":True,"same_H27_quotient":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args()
    p=payload();text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"rank":[25,26]},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
