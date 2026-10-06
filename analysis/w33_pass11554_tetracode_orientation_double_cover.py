#!/usr/bin/env python3
"""Pass 11554: the tetracode orientation lift is a non-split double cover.

Pass 10946 uses oriented representatives of P1(F3).  GL(2,3) acts on those
representatives, while projectivization gives PGL(2,3) ~= S4.  This producer
builds the central C2 cocycle from an explicit section and proves by F2 linear
algebra that it is not a coboundary: no globally S4-equivariant sign choice of
the four projective clock representatives exists.
"""
from __future__ import annotations
import itertools, json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/PART_W33_PASS11554_TETRACODE_ORIENTATION_DOUBLE_COVER.json"
P=3
R=((1,0),(0,1),(1,1),(1,2))

def mod(x): return x%P
def det(A): return mod(A[0]*A[3]-A[1]*A[2])
def mv(A,v): return (mod(A[0]*v[0]+A[1]*v[1]),mod(A[2]*v[0]+A[3]*v[1]))
def projective(v):
    s=pow(next(x for x in v if x),-1,3)
    return tuple(mod(s*x) for x in v)
def perm(A): return tuple(R.index(projective(mv(A,r))) for r in R)
def mul(A,B):
    return tuple(mod(sum(A[2*i+k]*B[2*k+j] for k in range(2)))
                 for i in range(2) for j in range(2))
def neg(A): return tuple(mod(-x) for x in A)
def compose(p,q): return tuple(p[q[i]] for i in range(4))
def parity(p): return (-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))

def rank2(A,ncols):
    A=[list(map(lambda x:x&1,row[:ncols])) for row in A]; r=0
    for c in range(ncols):
        q=next((i for i in range(r,len(A)) if A[i][c]),None)
        if q is None: continue
        A[r],A[q]=A[q],A[r]
        for i in range(len(A)):
            if i!=r and A[i][c]: A[i]=[x^y for x,y in zip(A[i],A[r])]
        r+=1
    return r

def main():
    G=[A for A in itertools.product(range(3),repeat=4) if det(A)]
    assert len(G)==48
    perms=sorted(set(perm(A) for A in G)); assert len(perms)==24
    fibers={p:[A for A in G if perm(A)==p] for p in perms}
    assert {len(v) for v in fibers.values()}=={2}
    assert all(v[1]==neg(v[0]) or v[0]==neg(v[1]) for v in fibers.values())

    section={p:min(fibers[p]) for p in perms}
    pi={p:i for i,p in enumerate(perms)}
    cocycle=[]; equations=[]
    for p in perms:
        for q in perms:
            r=compose(p,q)
            M=mul(section[p],section[q]); S=section[r]
            bit=0 if M==S else 1
            assert M==S or M==neg(S)
            cocycle.append({"p":list(p),"q":list(q),"c":bit})
            row=[0]*25
            row[pi[p]]^=1; row[pi[q]]^=1; row[pi[r]]^=1; row[24]=bit
            equations.append(row)
    r=rank2([x[:24] for x in equations],24)
    ra=rank2(equations,25)
    assert ra==r+1
    nonsplit=True

    assert all((1 if det(A)==1 else -1)==parity(perm(A)) for A in G)
    minusI=(2,0,0,2)
    assert perm(minusI)==(0,1,2,3)

    p10946=json.loads((ROOT/"data/w33_pass10946_clock_code_cone_objectwise.json").read_text())
    assert p10946["symmetry"]["GL2_order"]==48
    assert p10946["symmetry"]["same_projective_permutation_has_two_global_sign_lifts"] is True

    out={
      "schema":"w33.pass11554.tetracode_orientation_double_cover.v1",
      "status":"PASS_NON_SPLIT_GL23_TO_S4_ORIENTATION_EXTENSION","pass":11554,
      "extension":{"exact_sequence":"1 -> {+I,-I} -> GL(2,3) -> PGL(2,3)~=S4 -> 1",
        "GL2_order":48,"projective_image_order":24,"kernel_order":2,
        "splits":False,"section_cocycle_is_coboundary":False,
        "F2_coboundary_rank":r,"F2_augmented_rank":ra},
      "orientation_obstruction":{
        "statement":"There is no choice of one GL(2,3) lift for every S4 projective clock permutation that is closed under multiplication.",
        "central_minus_I":"acts trivially on P1(F3) but flips all oriented representatives simultaneously",
        "Pass10946_meaning":"the explicit tetracode evaluation signs are a genuine double-cover lift, not removable by an S4-equivariant gauge choice"},
      "determinant_character":{"law":"det(A) equals the sign of the induced S4 permutation",
        "SL2_preimage":"det=+1 maps onto A4"},
      "two_C2_firewall":"The central -I here sends v->-v and is killed by the Veronese map vv^T. It is distinct from the later nonsquare scaling S->-S that exchanges the two temporal cone/Weil sheets in Pass11550.",
      "cocycle_ones":sum(x["c"] for x in cocycle),
      "boundary":"Exact finite central-extension obstruction. It does not identify GL(2,3) with the binary octahedral group (they are distinct double covers of S4), nor turn the lift into a physical spin structure without further geometry."
    }
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":out["status"],"ranks":[r,ra],"cocycle_ones":out["cocycle_ones"]},indent=2))
if __name__=="__main__": main()
