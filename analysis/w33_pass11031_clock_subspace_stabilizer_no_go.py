#!/usr/bin/env python3
"""Pass 11031: exact stabilizer of the naive 3D clock augmentation."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass10975_e6_cubic_reynolds_clock_weld as P75
import w33_pass11021_minimal_signed_clock_carrier as P21
import w33_pass11022_binary_octahedral_clock_decomposition as P22

OUT=ROOT/"data/w33_pass11031_clock_subspace_stabilizer_no_go.json"
IDENT=(1,0,0,1)


def rank(vs):
    return sp.Matrix.hstack(*[sp.Matrix(v) for v in vs]).rank()


def order(g):
    h=IDENT
    for n in range(1,49):
        h=P22.mm(h,g)
        if h==IDENT:
            return n
    raise AssertionError("order > 48")
def payload():
    e6,dirs,rows,perms,section=P21.split_section()
    p_to_m={r[3]:tuple(map(int,r[0])) for r in rows}
    by_m={p_to_m[p]:(p,s) for p,s in section}
    group=tuple(sorted(by_m))
    assert len(group)==48

    fibres=[]
    indicators=[]
    for d in dirs:
        f=sorted(i for i,h in e6.items()
                 if h[:2]!=(0,0) and P75.norm_dir(h[:2])==d)
        assert len(f)==6
        v=[0]*27
        for i in f:
            v[i]=1
        fibres.append(f)
        indicators.append(tuple(v))

    aug=[tuple(a-b for a,b in zip(indicators[i],indicators[3]))
         for i in range(3)]
    assert rank(aug)==3
    B=sp.Matrix.hstack(*map(sp.Matrix,aug))

    stabilizer=[]
    induced=[]
    escape_ranks={}
    for g in group:
        images=[P21.act_vec(by_m[g],v) for v in aug]
        rr=rank(aug+images)
        escape_ranks[str(g)]=rr
        if rr==3:
            stabilizer.append(g)
            cols=[]
            for w in images:
                sol=next(iter(sp.linsolve((B,sp.Matrix(w)))))
                cols.append(sp.Matrix(sol))
            induced.append(sp.Matrix.hstack(*cols))
    assert stabilizer==[IDENT]
    nontrivial_escape=[escape_ranks[str(g)] for g in group if g!=IDENT]
    histogram={str(r):nontrivial_escape.count(r) for r in sorted(set(nontrivial_escape))}
    central=(2,0,0,2)
    assert escape_ranks[str(central)]==6
    full_closure=P21.orbit_span_rank(section,aug)[0]
    assert full_closure==24

    checks={
        "augmentation_rank_3": rank(aug)==3,
        "setwise_stabilizer_order_1": len(stabilizer)==1,
        "only_stabilizer_is_identity": stabilizer==[IDENT],
        "central_minusI_already_expands_to_rank6":
            escape_ranks[str(central)]==6,
        "full_signed_orbit_span_rank24": full_closure==24,
    }
    assert all(checks.values())
    return {
        "schema":"w33.pass11031.clock-subspace-stabilizer-no-go.v1",
        "status":"PASS",
        "headline":(
            "In the frozen signed H27 embedding the naive 3D clock augmentation "
            "has trivial setwise stabilizer inside the exact GL2(3) action. "
            "Every nonidentity signed symmetry sends at least one clock generator "
            "outside the 3D subspace; even central -I expands its joint span to 6, "
            "and full closure is 24."
        ),
        "carrier":{
            "exact_group":"GL2(3)",
            "group_order":len(group),
            "clock_augmentation_dimension":3,
            "coordinate_carrier_dimension":24,
        },
        "stabilizer":{
            "order":len(stabilizer),
            "elements":[list(g) for g in stabilizer],
            "induced_matrices":[list(map(int,list(M))) for M in induced],
            "nonidentity_joint_span_rank_histogram":histogram,
        },
        "symmetry_breaking_consequence":(
            "There is no nontrivial subgroup of the current exact signed GL2(3) "
            "whose every element preserves this particular coarse 3D clock "
            "augmentation: any such subgroup lies in the setwise stabilizer, "
            "which is trivial. Exact linear dynamics selecting this embedded "
            "3D image must therefore completely break the signed action, or use "
            "a quotient/readout map rather than an invariant-subspace mechanism."
        ),
        "boundary":(
            "This is basis-independent for the fixed 3D subspace but embedding-"
            "dependent: it does not rule out a different three-dimensional field "
            "inside V24, a nonlinear order parameter, or a different signed gauge "
            "combined with a correspondingly transformed physical readout."
        ),
        "checks":checks,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT)
    a=ap.parse_args()
    p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:
            raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(text)
    print(json.dumps({
        "status":p["status"],
        "stabilizer":p["stabilizer"]["order"],
        "escape_hist":p["stabilizer"]["nonidentity_joint_span_rank_histogram"],
    },sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
