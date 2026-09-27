#!/usr/bin/env python3
"""Pass 11039: exact rank of the 4x9 Bell-relative response tomography map."""
from __future__ import annotations
import argparse, itertools, json
from collections import Counter
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11039_bell_relative_tomography_rank.json"
Q=3


def canon(v):
    v=tuple(int(x)%Q for x in v)
    for x in v:
        if x:
            z=pow(x,-1,Q)
            return tuple(z*y%Q for y in v)
    raise ValueError


def omega(x,y):
    return (x[0]*y[1]-x[1]*y[0]-x[2]*y[3]+x[3]*y[2])%Q


def geometry():
    pts=sorted({canon(v) for v in itertools.product(range(Q),repeat=4) if any(v)})
    lines=set()
    for i,a in enumerate(pts):
        for b in pts[i+1:]:
            if omega(a,b): continue
            L=frozenset(canon(tuple((r*x+s*y)%Q for x,y in zip(a,b)))
                        for r,s in itertools.product(range(Q),repeat=2) if r or s)
            if len(L)==4: lines.add(L)
    return pts,sorted(lines,key=lambda L:sorted(L))
def rank_mod(M,p):
    A=[[int(M[i,j])%p for j in range(M.cols)] for i in range(M.rows)]
    r=0
    for c in range(M.cols):
        z=next((i for i in range(r,M.rows) if A[i][c]),None)
        if z is None: continue
        A[r],A[z]=A[z],A[r]
        iv=pow(A[r][c],-1,p); A[r]=[x*iv%p for x in A[r]]
        for i in range(M.rows):
            if i!=r and A[i][c]:
                q=A[i][c]
                A[i]=[(x-q*y)%p for x,y in zip(A[i],A[r])]
        r+=1
    return r


def eig_hist(M):
    out={}
    for x,m in M.eigenvals().items():
        out[str(x)]=int(m)
    return dict(sorted(out.items(),key=lambda kv:float(sp.N(sp.sympify(kv[0]))),reverse=True))


def payload():
    pts,lines=geometry()
    B=frozenset(canon((u[0],u[1],u[0],u[1]))
                for u in itertools.product(range(Q),repeat=2) if any(u))
    off=[p for p in pts if p not in B]
    trans=[L for L in lines if not(L&B)]
    assert len(off)==36 and len(trans)==27

    sectors={b:sorted(p for p in off if omega(p,b)==0) for b in sorted(B)}
    assert all(len(v)==9 for v in sectors.values())
    order=[p for b in sorted(B) for p in sectors[b]]
    R=sp.Matrix([[int(p in L) for L in trans] for p in order])
    assert set(map(int,list(R*sp.ones(27,1))))=={3}
    assert set(map(int,list(R.T*sp.ones(36,1))))=={4}
    rq=R.rank()
    assert rq==21
    mod_ranks={str(p):rank_mod(R,p) for p in (2,3,5,7,11)}
    assert set(mod_ranks.values())=={21}
    spectrum=eig_hist(R.T*R)
    assert spectrum=={"12":"1","6":"12","3":"8","0":"6"} or spectrum=={"12":1,"6":12,"3":8,"0":6}

    sector_ranks={}
    increments=[]
    prev=0
    for k in range(1,5):
        hist=Counter()
        for ss in itertools.combinations(range(4),k):
            rows=[j for s in ss for j in range(9*s,9*s+9)]
            hist[R[rows,:].rank()]+=1
        sector_ranks[str(k)]={str(a):b for a,b in sorted(hist.items())}
        unique=next(iter(hist)) if len(hist)==1 else None
        assert unique is not None
        increments.append(unique-prev); prev=unique
    assert increments==[9,6,4,2]

    G=R*R.T
    cross={}
    for i in range(4):
        for j in range(i+1,4):
            blk=G[9*i:9*i+9,9*j:9*j+9]
            vals=Counter(map(int,blk))
            rowsums=set(map(int,list(blk*sp.ones(9,1))))
            assert vals==Counter({0:54,1:27}) and rowsums=={3}
            assert blk.rank()==3
            cross[f"{i}-{j}"]={"rank":3,"ones":27,"zeros":54,"row_sum":3}
    checks={
      "36_offBell_by_27_transverse_incidence":R.shape==(36,27),
      "row_weight_3":set(map(int,list(R*sp.ones(27,1))))=={3},
      "column_weight_4":set(map(int,list(R.T*sp.ones(36,1))))=={4},
      "rational_rank_21":rq==21,
      "rank_21_mod_small_primes":set(mod_ranks.values())=={21},
      "Gram_spectrum_12_6_3_0":spectrum=={"12":1,"6":12,"3":8,"0":6},
      "sector_rank_increments_9_6_4_2":increments==[9,6,4,2],
      "all_cross_sector_blocks_rank3":all(x["rank"]==3 for x in cross.values()),
    }
    assert all(checks.values())
    return {
      "schema":"w33.pass11039.bell-relative-tomography-rank.v1",
      "status":"PASS",
      "headline":(
        "The proposed 4 x 9 Bell-relative response layout is exact but highly "
        "redundant. Incidence of the 36 off-Bell probes with the 27 transverse "
        "Lagrangian contexts is a 36x27 matrix of row weight 3 and column weight 4 "
        "with rank exactly 21 (over Q and over F2,F3,F5,F7,F11). Its Gram spectrum "
        "is 12^1, 6^12, 3^8, 0^6."
      ),
      "incidence":{
        "rows":"36 off-Bell W33 points ordered as four Bell-relative sectors of 9",
        "columns":"27 Lagrangian lines transverse to the Bell line",
        "shape":[36,27],
        "row_weight":3,
        "column_weight":4,
        "rank_Q":rq,
        "right_kernel_dimension":27-rq,
        "left_kernel_dimension":36-rq,
        "small_prime_ranks":mod_ranks,
        "Gram_RtR_spectrum":spectrum,
      },
      "sector_growth":{
        "rank_by_number_of_complete_9_point_sheets":sector_ranks,
        "successive_rank_increments":increments,
        "reading":(
          "one sheet carries rank 9; a second adds 6; a third adds 4; the fourth "
          "adds only 2. Thus the 36-response atlas has 15 exact linear redundancies "
          "and cannot reconstruct an arbitrary 27-parameter history vector without "
          "six additional gauge/normalization conditions."
        ),
      },
      "cross_sector_blocks":cross,
      "physical_reading":(
        "The 4x9 layout is best treated as an overcomplete response frame, not 36 "
        "independent observables. Its nonzero singular sectors have multiplicities "
        "1,12,8, suggestively close to the repo's finite-light-cone decomposition, "
        "but this packet does not identify those eigenspaces with physical memory "
        "orders without an explicit intertwiner."
      ),
      "boundary":(
        "This is incidence tomography on W33 operator directions. Pass 11038 shows "
        "that the off-Bell rays are not individually CPTP channels, so a laboratory "
        "protocol still needs a CP measurement/dilation map whose reconstructed "
        "linear responses realize this incidence frame."
      ),
      "checks":checks,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT)
    a=ap.parse_args(); p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text: raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(text)
    print(json.dumps({"status":p["status"],"rank":p["incidence"]["rank_Q"],
      "kernel":p["incidence"]["right_kernel_dimension"],
      "increments":p["sector_growth"]["successive_rank_increments"]},sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
