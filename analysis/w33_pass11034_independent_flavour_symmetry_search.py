#!/usr/bin/env python3
"""Pass 11034: exhaustive GL2(3) three-generation flavour/parity search."""
from __future__ import annotations
import argparse, itertools, json, sys
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11022_binary_octahedral_clock_decomposition as P22

OUT=ROOT/"data/w33_pass11034_independent_flavour_symmetry_search.json"


def invariant_mult(a,b,c,sizes):
    v=sp.simplify(sum(n*x*y*z for n,x,y,z in zip(sizes,a,b,c))/48)
    assert v.is_Integer
    return int(v)


def payload():
    rows=P22.exact_character_table()
    chars=dict(rows); sizes=P22.CLASS_SIZES
    names=[n for n,_ in rows]
    dims={n:int(chars[n][0]) for n in names}

    # All semisimple 3-dimensional complex characters of GL2(3).
    types=[]
    for mult in itertools.product(range(4),repeat=len(names)):
        if sum(mult[i]*dims[n] for i,n in enumerate(names))!=3:
            continue
        ch=tuple(sp.simplify(sum(mult[i]*chars[n][j] for i,n in enumerate(names)))
                 for j in range(8))
        label="+".join((str(mult[i])+"*" if mult[i]>1 else "")+n
                       for i,n in enumerate(names) if mult[i])
        types.append((label,ch,mult))
    assert len(types)==12
    ones=[("trivial",chars["trivial"]),("det",chars["det"])]
    N=len(types)
    tri={(i,j,k):invariant_mult(types[i][1],types[j][1],types[k][1],sizes)
         for i in range(N) for j in range(N) for k in range(N)}

    # Full-rank Yukawa pairing A x B x h requires B*h ~= A^*.
    pair={}
    for a in range(N):
      for h in range(2):
        pair[a,h]=[
          b for b in range(N)
          if all(sp.simplify(types[b][1][j]*ones[h][1][j]
                             -sp.conjugate(types[a][1][j]))==0
                 for j in range(8))
        ]
        assert len(pair[a,h])==1

    fullrank=0; safe=[]
    for Q in range(N):
      for Hu,Hd in itertools.product(range(2),repeat=2):
        for U in pair[Q,Hu]:
          for D in pair[Q,Hd]:
            for L in range(N):
              for E in pair[L,Hd]:
                fullrank+=1
                r=(tri[U,D,D],tri[Q,L,D],tri[L,L,E])
                if r==(0,0,0):
                    safe.append((Q,U,D,L,E,Hu,Hd))
    assert fullrank==576 and len(safe)==4

    safe_rows=[]
    for Q,U,D,L,E,Hu,Hd in safe:
        safe_rows.append({
          "Q":types[Q][0],"U":types[U][0],"D":types[D][0],
          "L":types[L][0],"E":types[E][0],
          "Hu":ones[Hu][0],"Hd":ones[Hd][0],
          "RPV_invariant_multiplicities":[0,0,0],
        })
    # Irreducible triplet-only family assignment: exact no-go.
    trip=[i for i,x in enumerate(types)
          if x[0] in ("standard3","standard3_twist")]
    natural=0; natural_safe=0; natural_patterns=set()
    for Q in trip:
      for Hu,Hd in itertools.product(range(2),repeat=2):
        U=pair[Q,Hu][0]; D=pair[Q,Hd][0]
        if U not in trip or D not in trip: continue
        for L in trip:
          E=pair[L,Hd][0]
          if E not in trip: continue
          natural+=1
          r=(tri[U,D,D],tri[Q,L,D],tri[L,L,E])
          natural_patterns.add(r)
          natural_safe+=r==(0,0,0)
    assert natural==16 and natural_safe==0 and natural_patterns=={(1,1,1)}

    # Every safe full-rank assignment factors through the 1D determinant quotient.
    abelian_labels={"3*trivial","3*det"}
    all_safe_abelian=all(
      all(row[f] in abelian_labels for f in ("Q","U","D","L","E"))
      for row in safe_rows)
    assert all_safe_abelian
    central_det=chars["det"][2]  # class 2 is central -I in frozen order
    assert central_det==1

    checks={
      "twelve_semisimple_generation_reps":len(types)==12,
      "fullrank_yukawa_assignments_576":fullrank==576,
      "irreducible_triplet_assignments_16":natural==16,
      "irreducible_triplets_always_allow_all_three_RPV":
          natural_safe==0 and natural_patterns=={(1,1,1)},
      "exactly_four_fullrank_RPV_free_assignments":len(safe)==4,
      "all_four_factor_through_determinant_Z2":all_safe_abelian,
      "determinant_Z2_is_independent_of_central_minusI":central_det==1,
    }
    assert all(checks.values())
    return {
      "schema":"w33.pass11034.independent-flavour-symmetry-search.v1",
      "status":"PASS",
      "headline":(
        "An exhaustive representation-ring search finds a sharp boundary for "
        "GL2(3) as an independent proton-protection symmetry. If each three-"
        "generation family is an irreducible triplet, all 16 full-rank Yukawa "
        "assignments also allow UDD, QLD and LLE once each. Across all 12 "
        "semisimple 3D family representations there are 576 full-rank Yukawa "
        "assignments and exactly four RPV-free ones; every survivor collapses "
        "to the generation-blind determinant Z2 quotient."
      ),
      "search_space":{
        "three_dimensional_semisimple_rep_types":[x[0] for x in types],
        "Higgs_reps":["trivial","det"],
        "forced_FI_singlet":"trivial, so its VEV preserves the candidate symmetry",
        "fullrank_yukawa_assignments":fullrank,
      },
      "irreducible_triplet_firewall":{
        "assignments":natural,
        "RPV_free":natural_safe,
        "RPV_multiplicity_pattern":[1,1,1],
        "reading":"nonabelian triplet flavour alone cannot forbid the regenerated RPV operators",
      },
      "surviving_independent_parities":safe_rows,
      "survivor_structure":{
        "count":len(safe_rows),
        "factor":"det: GL2(3) -> F3^* ~= Z2",
        "det_on_central_minusI":"+1",
        "independent_of_clock_central_parity":True,
        "generation_mixing":False,
      },
      "interpretation":(
        "The search did not find a genuinely nonabelian GL2(3) flavour rescue. "
        "It did uncover a narrower algebraic escape: a determinant Z2 that is "
        "independent of the clock's central -I and can remain unbroken when the "
        "forced FI singlet is trivial. Because it acts generation-blindly it is "
        "an extra parity, not a texture mechanism."
      ),
      "boundary":(
        "This is a representation-ring selection-rule search, not yet an "
        "embedding into the committed heterotic field labels. The decisive next "
        "test is whether any of the four determinant charge patterns is realized "
        "by an exact space-group/geometric automorphism of the 23/28 candidate "
        "orbifold vacua and remains unbroken by every required singlet VEV."
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
    print(json.dumps({"status":p["status"],
      "fullrank":p["search_space"]["fullrank_yukawa_assignments"],
      "triplet_safe":p["irreducible_triplet_firewall"]["RPV_free"],
      "det_parities":p["survivor_structure"]["count"]},sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
