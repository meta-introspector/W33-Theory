#!/usr/bin/env python3
"""Pass 11041: identify the diagonal-weld 12D intersection relative to 1+6+12+8."""
from __future__ import annotations
import argparse, itertools, json, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_e6_cubic_fourier54_alignment as ALIGN

OUT=ROOT/"data/w33_pass11041_phase_weld_pairwise12_dressing.json"


def rref(M,p):
    A=[[int(x)%p for x in row] for row in M]
    nr=len(A); nc=len(A[0]) if nr else 0; r=0; piv=[]
    for c in range(nc):
        z=next((i for i in range(r,nr) if A[i][c]),None)
        if z is None: continue
        A[r],A[z]=A[z],A[r]
        iv=pow(A[r][c],-1,p)
        A[r]=[x*iv%p for x in A[r]]
        for i in range(nr):
            if i!=r and A[i][c]:
                q=A[i][c]
                A[i]=[(x-q*y)%p for x,y in zip(A[i],A[r])]
        piv.append(c); r+=1
        if r==nr: break
    return A,piv


def rank(M,p):
    return len(rref(M,p)[1])


def nullspace(M,p):
    R,piv=rref(M,p)
    nc=len(M[0]) if M else 0
    free=[j for j in range(nc) if j not in piv]
    out=[]
    for f in free:
        x=[0]*nc; x[f]=1
        for rr,c in reversed(list(enumerate(piv))):
            x[c]=(-sum(R[rr][j]*x[j] for j in free))%p
        out.append(x)
    return out
def build_spaces(p,e6_to_h,records,background):
    roots=[x for x in range(2,p) if (x*x+x+1)%p==0]
    assert len(roots)==2
    omega=min(roots)

    def mm(A,B):
        return [[sum(A[i][k]*B[k][j] for k in range(len(B)))%p
                 for j in range(len(B[0]))] for i in range(len(A))]
    def mpow(A,e):
        R=[[1,0,0],[0,1,0],[0,0,1]]; B=A
        while e:
            if e&1: R=mm(R,B)
            B=mm(B,B); e//=2
        return R
    X=[[0,0,1],[1,0,0],[0,1,0]]
    Z=[[1,0,0],[0,omega,0],[0,0,pow(omega,2,p)]]
    def rho(h):
        a,b,c=h
        M=mm(mpow(Z,a),mpow(X,b)); s=pow(omega,c,p)
        return [[s*x%p for x in row] for row in M]

    labels=[]; cols=[]
    for t,r,i in itertools.product(range(3),repeat=3):
        labels.append((t,r,i)); col=[]
        for eid in range(27):
            M=rho(e6_to_h[eid])
            for phase in range(3):
                col.append(M[i][r]*pow(omega,(t*phase)%3,p)%p)
        cols.append(col)
    S=[[cols[j][i] for j in range(27)] for i in range(81)]
    assert rank(S,p)==27

    D=[[0]*81 for _ in range(81)]
    for u,x,o,c in records:
        D[o][x]=(D[o][x]+c*background[u])%p
    assert rank(D,p)==66

    # y is in the S1-coordinate intersection iff S*y lies in Im(D).
    # Equivalently every left-null vector of D annihilates S*y.
    DT=[list(row) for row in zip(*D)]
    left=nullspace(DT,p)
    assert len(left)==15
    C=[[sum(l[a]*S[a][j] for a in range(81))%p for j in range(27)]
       for l in left]
    K=nullspace(C,p)
    assert len(K)==12
    return labels,K
def sector_packet(labels,K,p):
    sectors={k:[j for j,lab in enumerate(labels)
                if sum(x!=0 for x in lab)==k] for k in range(4)}
    assert {k:len(v) for k,v in sectors.items()}=={0:1,1:6,2:12,3:8}

    proj={}
    intersections={}
    for k,inds in sectors.items():
        P=[[v[j] for j in inds] for v in K]
        proj[k]=rank(P,p)
        E=[[1 if j==idx else 0 for j in range(27)] for idx in inds]
        intersections[k]=len(K)+len(E)-rank(K+E,p)
    assert proj=={0:1,1:6,2:12,3:8}
    assert intersections=={0:0,1:0,2:0,3:0}

    # Since dim K = dim C2 = 12 and projection to C2 has rank 12,
    # K is the graph of unique linear maps C2 -> C0,C1,C3.
    return {
      "sector_dimensions":{"C0":1,"C1":6,"C2":12,"C3":8},
      "projection_ranks":{f"C{k}":proj[k] for k in range(4)},
      "pure_sector_intersection_dimensions":{f"C{k}":intersections[k] for k in range(4)},
      "pairwise_projection_isomorphism":proj[2]==12,
      "graph_component_ranks":{"C2_to_C0":proj[0],"C2_to_C1":proj[1],"C2_to_C3":proj[3]},
    }


def payload():
    bridge=json.loads((ROOT/"data/w33_e6id_current_h27_gauge_bridge.json").read_text())
    e6_to_h={int(i):tuple(map(int,h))
              for i,h in bridge["maps"]["e6id_to_current_H27_address"].items()}
    records=ALIGN.ordered_records()
    coords=[(e6_to_h[eid],phase) for eid in range(27) for phase in range(3)]
    backgrounds={
      "center_plus_external":[1+((h[2]+phase)%3) for h,phase in coords],
      "center_minus_external":[1+((h[2]-phase)%3) for h,phase in coords],
    }
    cert={}
    for name,bg in backgrounds.items():
        cert[name]={}
        for p in (103,109):
            labels,K=build_spaces(p,e6_to_h,records,bg)
            cert[name][str(p)]=sector_packet(labels,K,p)

    # Require the same structural verdict for both diagonal slopes and both split primes.
    for name in cert:
        for p,row in cert[name].items():
            assert row["projection_ranks"]=={"C0":1,"C1":6,"C2":12,"C3":8}
            assert row["pure_sector_intersection_dimensions"]=={
                "C0":0,"C1":0,"C2":0,"C3":0}
            assert row["pairwise_projection_isomorphism"]

    checks={
      "both_diagonal_welds_tested":set(cert)=={
          "center_plus_external","center_minus_external"},
      "split_primes_103_109":all(set(v)=={"103","109"} for v in cert.values()),
      "intersection_dimension_12":True,
      "correlation_dimensions_1_6_12_8":True,
      "full_projection_to_pairwise_C2":True,
      "zero_literal_intersection_with_pure_C2":True,
      "full_projection_to_C0_C1_C3":True,
      "same_verdict_both_weld_orientations":True,
    }
    assert all(checks.values())
    return {
      "schema":"w33.pass11041.phase-weld-pairwise12-dressing.v1",
      "status":"PASS",
      "headline":(
        "The diagonal phase weld's 12-dimensional S1 intersection is not literally "
        "the pure pairwise-correlation C2 sector. In the frozen S1 coefficient "
        "factorization (t,r,i) in F3^3, with dimensions 1+6+12+8 by correlation "
        "order, the weld intersection has zero intersection with every pure Ck. "
        "However its projection onto C2 has full rank 12, so it is canonically a "
        "graph over the entire pairwise 12, dressed by lower- and higher-order pieces."
      ),
      "canonical_factorization":{
        "S1_dimension":27,
        "basis_labels":"(external character t, regular multiplicity r, Schrodinger component i)",
        "local_index_set":"F3 x F3 x F3",
        "reference_value":0,
        "correlation_order_dimensions":{"C0":1,"C1":6,"C2":12,"C3":8},
      },
      "certificates":cert,
      "refined_phase_weld_mechanism":{
        "literal_statement_that_intersection_equals_C2":False,
        "pairwise_core_complete":True,
        "pairwise_projection_dimension":12,
        "dressing_component_ranks":{"scalar_C0":1,"single_C1":6,"triple_C3":8},
        "description":(
          "K12 is transverse to C2 and pi_C2:K12->C2 is an isomorphism. "
          "Thus K12={x+A0x+A1x+A3x : x in C2}. The diagonal weld activates "
          "all pairwise coordinates, but only as a dressed subspace coupling "
          "simultaneously to scalar, one-body, and three-body sectors."
        ),
      },
      "importance":(
        "This resolves the numerical 12 coincidence left open by the phase-weld "
        "paper and the uploaded temporal note. The missing pairwise 12 is present "
        "as a complete quotient coordinate, not as an isolated invariant block. "
        "That is a stronger mechanism for mediation between the 6- and 8-dimensional "
        "sectors, but it also forbids interpreting the weld as 'purely pairwise'."
      ),
      "boundary":(
        "The 1+6+12+8 split depends on the frozen tensor-factor basis and reference "
        "value 0 in each of the three S1 coefficient factors. The split-prime "
        "certificates prove the stated transversality over characteristic zero only "
        "to the extent that the relevant nonzero minors lift from Z[omega]; a fully "
        "symbolic Q(omega) graph matrix is still worth freezing."
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
        if not a.output.exists() or a.output.read_text()!=text:
            raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(text)
    row=p["certificates"]["center_plus_external"]["103"]
    print(json.dumps({"status":p["status"],"projections":row["projection_ranks"],
      "intersections":row["pure_sector_intersection_dimensions"]},sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
