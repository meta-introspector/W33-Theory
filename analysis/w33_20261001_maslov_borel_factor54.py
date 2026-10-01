#!/usr/bin/env python3
"""Explain the remaining factor 2 in the 256-pair Maslov chirality.

Pass 11222 proved 27 | chi but left chi/27 = +/-2 unexplained.
This producer identifies the effective relation stabilizer with the W33 flag
Borel U81:C2 and then resolves the defect already on its Sylow-3 core U81.
"""
from __future__ import annotations
import argparse, itertools, json, sys
from collections import Counter
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11222_why_54 as W
import w33_pass11209_maslov_chirality as M
import w33_pass11084_pgsp_u81_four_sheet_cover as N84

OUT=ROOT/"data/w33_20261001_maslov_borel_factor54.json"

def mkey(A): return tuple(int(x)%3 for x in A.ravel())
def pkey(A): return min(mkey(A),mkey((-A)%3))
def eye(n): return np.eye(n,dtype=np.int64)%3
def mm(A,B): return (A@B)%3
def closure(gens, projective=False):
    I=eye(gens[0].shape[0])
    key=pkey if projective else mkey
    found={key(I):I}; frontier=[I]
    while frontier:
        a=frontier.pop()
        for g in gens:
            b=mm(a,g); k=key(b)
            if k not in found:
                found[k]=b; frontier.append(b)
    return found

def porder(A):
    I=pkey(eye(A.shape[0])); x=eye(A.shape[0])
    for n in range(1,100):
        x=mm(x,A)
        if pkey(x)==I: return n
    raise AssertionError("projective order exceeded bound")

def projective_center(elements):
    vals=list(elements.values())
    out=set()
    for a in vals:
        if all(pkey(mm(a,b))==pkey(mm(b,a)) for b in vals):
            out.add(pkey(a))
    return out

SOURCE_GENS=[
np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],
          [0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,1,1]],dtype=np.int64),
np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],
          [0,0,1,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]],dtype=np.int64),
np.array([[1,0,0,0,0,0],[1,1,0,0,0,0],[0,0,1,0,0,0],
          [0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]],dtype=np.int64),
2*eye(6),
np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0,0,1,0],
          [0,0,0,0,0,1],[0,0,1,0,0,0],[0,0,0,1,0,0]],dtype=np.int64),
np.array([[0,0,1,0,0,0],[0,0,0,1,0,0],[1,0,0,0,0,0],
          [0,1,0,0,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]],dtype=np.int64),
]
TARGET_GENS=[
np.array([[1,0,2,2],[0,1,2,2],[0,0,1,0],[0,0,0,1]],dtype=np.int64),
np.array([[1,0,2,1],[0,1,1,2],[0,0,1,0],[0,0,0,1]],dtype=np.int64),
np.array([[1,0,0,0],[0,1,0,2],[0,0,1,0],[0,0,0,1]],dtype=np.int64),
eye(4),
np.array([[1,0,0,1],[0,2,2,0],[0,0,1,0],[0,0,0,2]],dtype=np.int64),
np.array([[1,1,0,1],[0,2,2,0],[0,0,1,0],[0,0,1,2]],dtype=np.int64),
]

def target_borel():
    rootgens=[(eye(4)+A)%3 for A in N84.X]
    U=closure(rootgens)
    B=closure(rootgens+[N84.T])
    assert len(U)==81 and len(B)==162
    return U,B
def verify_projective_borel_isomorphism(raw):
    rawkeys={mkey(g) for g in raw}
    assert all(mkey(g) in rawkeys for g in SOURCE_GENS)
    generated=closure(SOURCE_GENS)
    assert len(generated)==324 and set(generated)==rawkeys
    target_U,target_B=target_borel()
    mapping={pkey(eye(6)):mkey(eye(4))}
    reps={mkey(eye(4)):eye(4)}
    frontier=[(eye(6),eye(4))]
    while frontier:
        a,b=frontier.pop()
        for sg,tg in zip(SOURCE_GENS,TARGET_GENS):
            aa=mm(a,sg); bb=mm(b,tg)
            ka,kb=pkey(aa),mkey(bb)
            if ka in mapping:
                assert mapping[ka]==kb
            else:
                mapping[ka]=kb; reps[kb]=bb; frontier.append((aa,bb))
    assert len(mapping)==162 and len(set(mapping.values()))==162
    assert set(mapping.values())==set(target_B)
    return target_U,target_B,{
        "domain_projective_order":len(mapping),
        "image_order":len(set(mapping.values())),
        "bijective":True,
        "gap_smallgroup_projective":[162,10],
        "gap_smallgroup_raw":[324,68],
        "gap_structure_projective":"((C3 x C3 x C3) : C3) : C2",
        "gap_structure_raw":"C2 x (((C3 x C3 x C3) : C3) : C2)",
    }
def action_tables(gs,P0,PF):
    k0={W.lkey(L):i for i,L in enumerate(P0)}
    kF={W.lkey(L):i for i,L in enumerate(PF)}
    return (
      [[k0[W.lkey((g@L)%3)] for L in P0] for g in gs],
      [[kF[W.lkey((g@L)%3)] for L in PF] for g in gs])

def odd_triples(P0,PF):
    i,m=(a.ravel() for a in np.meshgrid(np.arange(64),np.arange(64),indexing="ij"))
    dim=(6-M.rank_batch(np.concatenate([P0[i],PF[m]],2))).reshape(64,64)
    prof0=[tuple(sorted(r)) for r in dim.tolist()]
    pairs=[(a,b) for a in range(64) for b in range(64) if prof0[a]>prof0[b]]
    i1=np.repeat([a for a,_ in pairs],64)
    i2=np.repeat([b for _,b in pairs],64)
    mmidx=np.tile(np.arange(64),len(pairs))
    cls=M.kashiwara(P0[i1],P0[i2],PF[mmidx])
    return {(a,b,c):int(k) for a,b,c,k in
            zip(i1.tolist(),i2.tolist(),mmidx.tolist(),cls.tolist()) if int(k)%2}

def orbit(t,A0,AF):
    return frozenset((A0[j][t[0]],A0[j][t[1]],AF[j][t[2]])
                     for j in range(len(A0)))
def orbit_partition(cl,A0,AF):
    seen=set(); out=[]
    for t,c in cl.items():
        if t in seen: continue
        o=orbit(t,A0,AF); seen|=o
        assert {cl[x] for x in o}=={c}
        out.append((o,c))
    return out

def census(orbits):
    return Counter((len(o),c) for o,c in orbits)

def imbalance(C):
    sizes=sorted({s for s,_ in C})
    return {str(s):C.get((s,1),0)-C.get((s,3),0)
            for s in sizes if C.get((s,1),0)!=C.get((s,3),0)}

def analyse_relation(B):
    raw=W.relation_stabilizer(B)
    assert len(raw)==324
    proj={}
    for g in raw: proj.setdefault(pkey(g),g)
    assert len(proj)==162
    vals=list(proj.values())
    U=[g for g in vals if porder(g) in (1,3,9)]
    assert len(U)==81
    center=projective_center({pkey(g):g for g in U})
    assert len(center)==3
    P0=M.product_lagrangians(eye(6))
    PF=M.product_lagrangians(B%3)
    B0,BF=action_tables(vals,P0,PF)
    U0,UF=action_tables(U,P0,PF)
    cl=odd_triples(P0,PF)
    borb=orbit_partition(cl,B0,BF)
    uorb=orbit_partition(cl,U0,UF)
    bc,uc=census(borb),census(uorb)
    uid={}
    for j,(o,_c) in enumerate(uorb):
        for x in o: uid[x]=j
    split=Counter()
    for o,c in borb:
        parts={uid[x] for x in o}
        split[(len(o),c,len(parts),tuple(sorted(len(uorb[j][0]) for j in parts)))]+=1

    ukeys=[pkey(g) for g in U]
    stab_types=Counter()
    for o,c in uorb:
        if len(o)!=27: continue
        t=next(iter(o))
        st={ukeys[j] for j in range(len(U))
            if (U0[j][t[0]],U0[j][t[1]],UF[j][t[2]])==t}
        assert len(st)==3
        stab_types[(c,"center" if st==center else "noncentral")]+=1
    chi=sum(len(o)*(1 if c==1 else -1) for o,c in uorb)
    return {
      "raw_stabilizer_order":len(raw),
      "effective_projective_order":len(proj),
      "scalar_kernel":["+I","-I"],
      "sylow3_order":len(U),
      "sylow3_center_order":len(center),
      "B_census":{f"{s}:{c}":n for (s,c),n in sorted(bc.items())},
      "U_census":{f"{s}:{c}":n for (s,c),n in sorted(uc.items())},
      "B_imbalance":imbalance(bc),
      "U_imbalance":imbalance(uc),
      "B_to_U_splitting":{str(k):v for k,v in sorted(split.items())},
      "size27_stabilizer_types":{str(k):v for k,v in sorted(stab_types.items())},
      "chi_from_U_orbits":chi,
      "factor54_identity":"|U81|-|U81/C3_noncentral| = 81-27 = 54",
    }
def payload():
    reps=json.loads((ROOT/"data/w33_pass11212_one_chirality.json").read_text())["representatives"]
    prior=json.loads((ROOT/"data/w33_pass11222_why_54.json").read_text())
    B6=np.array(reps["6"]["B"],dtype=np.int64)
    raw6=W.relation_stabilizer(B6)
    _u,_b,iso=verify_projective_borel_isomorphism(raw6)
    results={}
    for label in ("6","7"):
        row=analyse_relation(np.array(reps[label]["B"],dtype=np.int64))
        row["twist"]=prior["orbitals"][label]["twist"]
        row["prior_chi"]=prior["orbitals"][label]["chi"]
        assert row["chi_from_U_orbits"]==row["prior_chi"]
        results[label]=row
    assert results["6"]["U_imbalance"]=={"27":1,"81":-1}
    assert results["7"]["U_imbalance"]=={"27":-1,"81":1}
    assert results["6"]["chi_from_U_orbits"]==-54
    assert results["7"]["chi_from_U_orbits"]==54
    assert all("noncentral" in k for r in results.values()
               for k in r["size27_stabilizer_types"])
    return {
      "schema":"w33.20261001.maslov-borel-factor54.v1",
      "status":"PASS_MASLOV_54_IS_THE_U81_FREE_VERSUS_C3_STABILIZED_ORBIT_DEFECT",
      "headline":(
        "The unexplained factor two in the 256-pair Maslov chirality is resolved. "
        "The raw order-324 relation stabilizer has an ineffective scalar +/-I; its "
        "projective quotient is the W33 flag Borel U81:C2. Restricting the odd "
        "Maslov triples to U81 leaves exactly one unmatched 27-orbit and one "
        "oppositely signed unmatched free 81-orbit. Hence chi=+/- (81-27)=+/-54. "
        "The 27-orbit stabilizer is a noncentral C3. The residual Borel C2 only "
        "fuses U81 orbits and is not the source of chirality."
      ),
      "effective_stabilizer_isomorphism":iso,
      "relations":results,
      "gap_independent_audit":{
        "raw_idgroup":[324,68],"projective_idgroup":[162,10],
        "projective_borel_idgroup":[162,10],
        "raw_direct_product":True,
      },
      "correction_firewall":{
        "false_order324_identification":
          "Gamma is not N_PGSp(U81): their centers, derived groups and order censuses differ.",
        "correct_identification":
          "Gamma/{+/-I} is isomorphic to N_PSp(U81)=U81:C2.",
        "factor2_not_deck_bit":
          "The chirality defect already exists on U81; the Borel C2 only fuses pairs of U81 orbits.",
        "size27_not_center_quotient":
          "Its stabilizer is a noncentral C3, not Z(U81).",
      },
      "boundary":(
        "This is an exact finite-group/orbit explanation of the integer 54 in this "
        "Maslov chirality. It does not make 54 a universal physical constant, a "
        "time rate, a coupling, or a continuum prediction."
      ),
      "parents":[
        "data/w33_pass11222_why_54.json",
        "data/w33_pass11212_one_chirality.json",
        "data/w33_pass11084_pgsp_u81_four_sheet_cover.json",
        "data/w33_pass11070_dual_parabolic_flag_diamond.json",
      ],
    }
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT)
    a=ap.parse_args()
    p=payload(); text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:
            raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(text,encoding="utf-8")
    print(json.dumps({
      "status":p["status"],
      "chi6":p["relations"]["6"]["chi_from_U_orbits"],
      "chi7":p["relations"]["7"]["chi_from_U_orbits"],
      "idgroup":p["effective_stabilizer_isomorphism"]["gap_smallgroup_projective"],
    },sort_keys=True))
if __name__=="__main__": raise SystemExit(main())
