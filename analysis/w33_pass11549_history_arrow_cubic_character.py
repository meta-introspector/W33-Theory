#!/usr/bin/env python3
"""Pass 11549: closed-form history arrow in Hamming coordinates.

Pass 11547 identifies the 27 histories with F_3^3 Hamming addresses and null
steps with the eight full-support vectors d in {+/-1}^3.  Pass 11548 identifies
the repo PSp linear group with the even-sign subgroup ker(sigma).

Define the cubic sign character on an oriented null step
    chi(d) = product_i sign(d_i),  d_i in {1,2} = {+1,-1}.
Then chi(-d)=-chi(d).

This verifier proves:
  * the old 108-edge invariant cycle is exactly -chi(H(s_j-s_i)) in the
    repository's canonical i<j edge ordering;
  * equivalently, reversing the arbitrary global sign, the arrow cycle is the
    local directed-edge rule chi(d);
  * every null edge belongs to exactly one temporal triangle;
  * the 36 triangles are exactly 4 projective null directions x 9 affine
    lines, with the four oriented directions selected by chi=+1;
  * orienting every affine line along its chi=+1 direction and summing the 36
    triangle boundaries gives the chi edge cycle;
  * repo PSp=ker(sigma) preserves chi; the outer sigma=-1 coset flips it.

No thermodynamic/cosmological arrow or continuum causal structure is inferred.
"""
from __future__ import annotations

import itertools
import json
from collections import Counter
from pathlib import Path

P=3
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data"/"PART_W33_PASS11549_HISTORY_ARROW_CUBIC_CHARACTER.json"
V=list(itertools.product(range(3),repeat=3))
VID={v:i for i,v in enumerate(V)}


def mod(x): return x%P
def add(u,v): return tuple(mod(u[i]+v[i]) for i in range(3))
def sub(u,v): return tuple(mod(u[i]-v[i]) for i in range(3))


def qsym(v):
    a,b,c=v
    return mod(a*c-b*b)


def Hsym(v):
    a,b,c=v
    return (mod(2*a+b+c),mod(2*a+2*b+c),mod(2*a+2*c))


HINV={Hsym(v):v for v in V}
assert len(HINV)==27


def hamming_chi(d):
    """Map 1 -> +1, 2=-1 -> -1 and multiply the three coordinates."""
    assert all(x in (1,2) for x in d)
    out=1
    for x in d:
        out*=1 if x==1 else -1
    return out


def null_edge(i,j):
    return i!=j and qsym(sub(V[j],V[i]))==0


EDGES=[(i,j) for i in range(27) for j in range(i+1,27) if null_edge(i,j)]
EIDX={e:i for i,e in enumerate(EDGES)}
TRIS=[t for t in itertools.combinations(range(27),3)
      if all(null_edge(i,j) for i,j in itertools.combinations(t,2))]
TIDX={t:i for i,t in enumerate(TRIS)}


def edge_add(vec,i,j,coef=1):
    if i<j: vec[EIDX[(i,j)]]+=coef
    else: vec[EIDX[(j,i)]]-=coef


def sorted_triangle_boundary(t):
    a,b,c=t
    v=[0]*len(EDGES)
    edge_add(v,b,c,+1)
    edge_add(v,a,c,-1)
    edge_add(v,a,b,+1)
    return v


def signed_permutation_group():
    out=[]
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1,2),repeat=3):
            M=[[0]*3 for _ in range(3)]
            for i,j in enumerate(perm): M[i][j]=signs[i]
            out.append(tuple(tuple(r) for r in M))
    assert len(set(out))==48
    return sorted(set(out))


def mv(M,v):
    return tuple(mod(sum(M[i][j]*v[j] for j in range(3))) for i in range(3))


def sign_product(M):
    out=1
    for row in M:
        nz=[x for x in row if x]
        assert len(nz)==1
        out*=1 if nz[0]==1 else -1
    return out


def main():
    assert len(EDGES)==108 and len(TRIS)==36

    # Every null edge is contained in exactly one triangle.
    edge_triangle_count=Counter()
    for t in TRIS:
        for e in itertools.combinations(t,2):
            edge_triangle_count[tuple(sorted(e))]+=1
    assert set(edge_triangle_count.values())=={1}
    assert len(edge_triangle_count)==108

    # Closed local formula in canonical i<j edge ordering.
    chi_cycle=[]
    for i,j in EDGES:
        d=Hsym(sub(V[j],V[i]))
        assert all(d)
        chi_cycle.append(hamming_chi(d))
    assert Counter(chi_cycle)==Counter({-1:72,+1:36})

    # Compare exactly with the old frozen PSp-invariant cycle.
    old=json.loads(
        (ROOT/"data"/"w33_20260924_history_invariant_cycle_orientation.json").read_text()
    )
    old_cycle=old["invariant_cycle"]["edge_coefficients"]
    assert len(old_cycle)==108
    assert old_cycle==[-x for x in chi_cycle]

    # The old all-+1 sorted triangle chain gives that same old global sign.
    old_from_triangles=[0]*108
    for t in TRIS:
        b=sorted_triangle_boundary(t)
        for k,x in enumerate(b): old_from_triangles[k]+=x
    assert old_from_triangles==old_cycle

    # Four oriented body diagonals chi=+1.
    oriented_dirs=[d for d in V if all(d) and hamming_chi(d)==1]
    assert len(oriented_dirs)==4
    assert all(hamming_chi(tuple(mod(-x) for x in d))==-1 for d in oriented_dirs)

    # Their affine cosets are exactly the 36 temporal triangles: 9 per direction.
    line_records=[]
    seen=set()
    chi_triangle_cycle=[0]*108
    per_direction=Counter()
    for dh in oriented_dirs:
        d=HINV[dh]
        direction_lines=set()
        for x in V:
            line=tuple(sorted((
                VID[x],
                VID[add(x,d)],
                VID[add(x,add(d,d))],
            )))
            direction_lines.add(line)
        assert len(direction_lines)==9
        per_direction[str(dh)]=9
        for line in sorted(direction_lines):
            assert line in TIDX
            assert line not in seen
            seen.add(line)

            # Start at the smallest-index vertex only to make a deterministic
            # cyclic representative; cyclic rotation does not change boundary.
            x0=V[line[0]]
            x1=add(x0,d)
            x2=add(x1,d)
            assert add(x2,d)==x0

            edge_add(chi_triangle_cycle,VID[x0],VID[x1],+1)
            edge_add(chi_triangle_cycle,VID[x1],VID[x2],+1)
            edge_add(chi_triangle_cycle,VID[x2],VID[x0],+1)
            line_records.append({
                "hamming_direction":list(dh),
                "sym2_direction":list(d),
                "triangle":list(line),
            })

    assert len(seen)==36 and seen==set(TRIS)
    assert chi_triangle_cycle==chi_cycle

    # The local character transforms exactly by the Pass-11548 sigma character.
    O=signed_permutation_group()
    PSp=[M for M in O if sign_product(M)==1]
    outer=[M for M in O if sign_product(M)==-1]
    assert len(PSp)==len(outer)==24
    null_steps=[d for d in V if all(d)]
    for M in PSp:
        assert all(hamming_chi(mv(M,d))==hamming_chi(d) for d in null_steps)
    for M in outer:
        assert all(hamming_chi(mv(M,d))==-hamming_chi(d) for d in null_steps)

    # Translation invariance upgrades 24/24 to the affine 648/648 split.
    correction=json.loads(
        (ROOT/"data"/"PART_W33_PASS11548_HISTORY_ORIENTATION_SQUARE_CORRECTION.json").read_text()
    )
    assert correction["corrected_orientation_dictionary"]["global_history_outer_bit"]["character"].startswith("sigma")

    out={
      "schema":"w33.pass11549.history_arrow_cubic_character.v1",
      "status":"PASS_CLOSED_FORM_HISTORY_ARROW_CUBIC_CHARACTER",
      "pass":11549,
      "hamming_arrow":{
        "null_steps":"d in {+/-1}^3",
        "local_character":"chi(d)=product_i sign(d_i)",
        "antisymmetry":"chi(-d)=-chi(d)",
        "canonical_edge_formula":"for repo edge i<j, old coefficient = -chi(Hsym(s_j-s_i))",
        "exact_old_edge_matches":108,
        "global_sign_note":"the opposite overall sign is the same one-dimensional orientation line",
      },
      "triangle_foliation":{
        "temporal_triangles":36,
        "edges":108,
        "each_edge_in_exactly_one_triangle":True,
        "oriented_projective_directions":4,
        "affine_lines_per_direction":9,
        "factorization":"36 = 4 * 9",
        "oriented_hamming_directions":[list(d) for d in oriented_dirs],
        "direction_line_counts":dict(sorted(per_direction.items())),
        "all_sorted_triangle_coefficients_in_old_chain":"+1",
        "chi_oriented_triangle_boundaries_equal_chi_edge_cycle":True,
      },
      "symmetry":{
        "repo_PSp_linear":"ker(sigma), order 24",
        "repo_PSp_affine_order":648,
        "PSp_action_on_chi":"preserves pointwise as a directed-edge rule",
        "outer_linear_coset_order":24,
        "outer_affine_coset_size":648,
        "outer_action_on_chi":"negates",
        "closed_form_reason":"chi(Md)=sigma(M) chi(d)",
      },
      "old_certificate_weld":{
        "source":"data/w33_20260924_history_invariant_cycle_orientation.json",
        "old_support_size":old["invariant_cycle"]["support_size"],
        "old_PSp_fixed":old["invariant_cycle"]["PSp_fixed_pointwise"],
        "old_outer_flips":old["PGSp_outer_action"]["outer_coset_flips_orientation"],
        "replaces_group_average_with_local_formula":True,
      },
      "interpretation":"The finite history orientation is the chirality of a body-diagonal Hamming step: whether an even or odd number of the three trits flips sign. Four positive body diagonals orient four families of nine affine temporal triangles.",
      "boundary":"Exact finite directed-edge/triangle theorem. It simplifies the repository's finite history-orientation certificate; it does not establish a thermodynamic arrow, cosmological initial condition, continuum causal orientation, or observed CPT dynamics.",
    }
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
      "status":out["status"],
      "edge_matches":108,
      "triangles":"4*9=36",
      "PSp":"preserves chi",
      "outer":"flips chi",
    },indent=2,sort_keys=True))


if __name__=="__main__":
    main()
