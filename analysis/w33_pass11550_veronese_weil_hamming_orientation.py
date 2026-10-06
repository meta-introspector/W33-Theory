#!/usr/bin/env python3
"""Pass 11550: Veronese sheet = Weil-phase sheet = Hamming arrow sheet.

Inputs already frozen in the repo:
  * Pass 10946: four P1(F3) clock rays map by nu(v)=vv^T to four rank-one
    temporal null matrices.
  * 2026-09-24 fixed9 certificate: the eight oriented rank-one dual null modes
    split 4+4 under PSp, with normalized Weil phases -i and +i.
  * Passes 11547-11549: Hamming coordinates and cubic null-step character
    chi(z)=sign(z1)sign(z2)sign(z3).

This verifier proves objectwise:
  * the four Veronese matrices vv^T are exactly the old -i Weil orbit;
  * their negatives are exactly the +i orbit;
  * Hamming chi=-1 on the Veronese sheet and +1 on its negative;
  * hence Weil_phase = i*chi on all eight oriented rank-one null modes;
  * nonsquare scaling -1 swaps the two sheets, flips chi, and conjugates
    the Weil phase;
  * p(z)=z1*z2*z3 is a degree-3 polynomial invariant of the actual repo
    PSp linear group W(D3)=even signed permutations and transforms by sigma
    under the full O(3,3) signed-permutation group.

The Veronese map kills v -> -v, so this cone-sheet orientation does NOT remove
Pass10946's explicit oriented-representative choice for tetracode evaluation.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

P=3
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data"/"PART_W33_PASS11550_VERONESE_WEIL_HAMMING_ORIENTATION.json"
V=list(itertools.product(range(3),repeat=3))


def mod(x): return x%P
def neg(v): return tuple(mod(-x) for x in v)


def Hsym(v):
    a,b,c=v
    return (mod(2*a+b+c),mod(2*a+2*b+c),mod(2*a+2*c))


def chi(z):
    assert all(x in (1,2) for x in z)
    out=1
    for x in z: out*=1 if x==1 else -1
    return out


def cubic(z):
    return mod(z[0]*z[1]*z[2])


def veronese(v):
    x,y=v
    return (mod(x*x),mod(x*y),mod(y*y))


def mv(M,v):
    return tuple(mod(sum(M[i][j]*v[j] for j in range(3))) for i in range(3))


def signed_permutation_group():
    out=[]
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1,2),repeat=3):
            M=[[0]*3 for _ in range(3)]
            for i,j in enumerate(perm): M[i][j]=signs[i]
            out.append(tuple(tuple(r) for r in M))
    return sorted(set(out))


def sign_product(M):
    out=1
    for row in M:
        nz=[x for x in row if x]
        assert len(nz)==1
        out*=1 if nz[0]==1 else -1
    return out


def matrix_tuple(M):
    return (int(M[0][0]),int(M[0][1]),int(M[1][1]))


def main():
    p10946=json.loads(
        (ROOT/"data"/"w33_pass10946_clock_code_cone_objectwise.json").read_text()
    )
    fixed9=json.loads(
        (ROOT/"data"/"w33_20260924_fixed9_chiral_null_fourplusfour.json").read_text()
    )
    p11548=json.loads(
        (ROOT/"data"/"PART_W33_PASS11548_HISTORY_ORIENTATION_SQUARE_CORRECTION.json").read_text()
    )
    p11549=json.loads(
        (ROOT/"data"/"PART_W33_PASS11549_HISTORY_ARROW_CUBIC_CHARACTER.json").read_text()
    )
    assert p10946["status"]=="PASS_CLOCK_CODE_CONE_OBJECTWISE_THEOREM"
    assert fixed9["status"]=="PASS_FIXED9_IS_VACUUM_PLUS_TWO_ORIENTED_NULL_FOURS"
    assert p11548["status"]=="PASS_OBJECTWISE_CORRECTION_OF_PASS11540_CHARACTER_ASSIGNMENT"
    assert p11549["status"]=="PASS_CLOSED_FORM_HISTORY_ARROW_CUBIC_CHARACTER"

    projective_rays=[tuple(x) for x in p10946["clock_line"]["projective_rays"]]
    oriented_eval=[tuple(x) for x in p10946["clock_line"]["oriented_evaluation_representatives"]]
    cone_rows=p10946["cone_dictionary"]["coordinate_to_existing_geometry"]

    # Veronese is independent of the +/- representative of a projective ray.
    for ray,orep,row in zip(projective_rays,oriented_eval,cone_rows):
        assert veronese(ray)==veronese(orep)
        assert tuple(row["symmetric_square_null_ray"])==veronese(ray)

    ver_sheet={veronese(r) for r in projective_rays}
    assert len(ver_sheet)==4
    neg_sheet={neg(s) for s in ver_sheet}
    assert len(neg_sheet)==4 and ver_sheet.isdisjoint(neg_sheet)

    # Hamming cubic sheet labels.
    ver_h={s:Hsym(s) for s in ver_sheet}
    neg_h={s:Hsym(s) for s in neg_sheet}
    assert {chi(z) for z in ver_h.values()}=={-1}
    assert {chi(z) for z in neg_h.values()}=={+1}

    # The old Weil-phase rows store the actual symmetric matrices objectwise.
    minus_i=set()
    plus_i=set()
    row_table=[]
    for row in fixed9["rows"]:
        s=matrix_tuple(row["dual_symmetric_matrix"])
        z=Hsym(s)
        ch=chi(z)
        phase=row["phase"]
        expected="-i" if ch==-1 else "i"
        assert phase==expected
        (minus_i if phase=="-i" else plus_i).add(s)
        row_table.append({
            "symmetric_matrix_coords":list(s),
            "hamming_coords":list(z),
            "chi":ch,
            "old_weil_phase":phase,
            "formula":"i*chi",
        })

    assert minus_i==ver_sheet
    assert plus_i==neg_sheet
    assert neg_sheet=={neg(s) for s in minus_i}

    # Outer nonsquare scaling -1 swaps every objectwise row with its conjugate.
    phase_by_s={matrix_tuple(r["dual_symmetric_matrix"]):r["phase"] for r in fixed9["rows"]}
    conjugate={"i":"-i","-i":"i"}
    for s,phase in phase_by_s.items():
        assert phase_by_s[neg(s)]==conjugate[phase]
        assert chi(Hsym(neg(s)))==-chi(Hsym(s))

    # W(D3) invariant-theory statement, verified directly over all F3^3.
    O=signed_permutation_group()
    PSp=[M for M in O if sign_product(M)==1]
    outer=[M for M in O if sign_product(M)==-1]
    assert len(PSp)==len(outer)==24
    for M in PSp:
        assert all(cubic(mv(M,z))==cubic(z) for z in V)
    for M in outer:
        assert all(cubic(mv(M,z))==mod(-cubic(z)) for z in V)

    # On the null cube, cubic value 1/2 is exactly chi +1/-1.
    for z in V:
        if all(z):
            assert cubic(z)==(1 if chi(z)==1 else 2)

    out={
      "schema":"w33.pass11550.veronese_weil_hamming_orientation.v1",
      "status":"PASS_VERONESE_WEIL_HAMMING_ORIENTATION_WELD",
      "pass":11550,
      "veronese_sheet":{
        "map":"nu(x,y)=(x^2,x*y,y^2)",
        "projective_rays":[list(x) for x in projective_rays],
        "rank_one_matrices":[list(x) for x in sorted(ver_sheet)],
        "hamming_images":[list(ver_h[x]) for x in sorted(ver_sheet)],
        "chi_value":-1,
        "old_weil_phase":"-i",
      },
      "negative_veronese_sheet":{
        "rank_one_nonsquare_scaled_matrices":[list(x) for x in sorted(neg_sheet)],
        "hamming_images":[list(neg_h[x]) for x in sorted(neg_sheet)],
        "chi_value":+1,
        "old_weil_phase":"+i",
      },
      "closed_formula":{
        "all_eight_rows":row_table,
        "Weil_phase_equals":"i * chi",
        "objectwise_matches":8,
        "outer_minus_one":"s -> -s; chi -> -chi; +i <-> -i",
      },
      "cubic_invariant":{
        "polynomial":"p(z)=z1*z2*z3 over F3",
        "repo_PSp_linear":"W(D3)=ker(sigma)",
        "PSp_law":"p(Mz)=p(z)",
        "full_O_law":"p(Mz)=sigma(M)*p(z)",
        "null_cube_relation":"p=1 corresponds chi=+1; p=2=-1 corresponds chi=-1",
        "standard_invariant_theory_note":"type D_n Weyl invariants have degrees 2,4,...,2n-2,n; for D3 the degree-3 invariant is the coordinate product",
      },
      "tetracode_orientation_firewall":{
        "veronese_kills_projective_sign":True,
        "statement":"nu(v)=nu(-v), so the cone sheet selected by the rank-one Veronese lift does not determine the oriented P1 representatives used to freeze the tetracode evaluation signs in Pass10946",
        "what_is_new":"the cone/Weil/time-orientation sheet is now canonical objectwise and is compatible with, but distinct from, the code-coordinate orientation lift",
      },
      "interpretation":"The four projective clock/null directions acquire two oriented rank-one sheets. The square/Veronese sheet is exactly the -i spectral sheet; nonsquare scaling gives its +i conjugate. The finite history arrow chi and the old normalized Weil phase are the same Z2 datum in real-sign and phase language.",
      "boundary":"Exact finite F3 coding/quadratic/Weil-phase weld. It does not turn the sign sheet into thermodynamic time, derive continuum CPT, or select a physical tetracode sign convention.",
    }
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":out["status"],
        "Veronese":"chi=-1 <-> -i",
        "negative":"chi=+1 <-> +i",
        "formula":"Weil=i*chi",
        "rows":8,
    },indent=2,sort_keys=True))


if __name__=="__main__":
    main()
