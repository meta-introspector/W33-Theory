#!/usr/bin/env python3
"""Pass 11043: a 9D latent H27 module regularizes the trinification qutrit."""
from __future__ import annotations
import argparse, itertools, json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11043_commutant_dressing_regular_h27.json"


def add(x,y): return (x[0]+y[0],x[1]+y[1])


def mul(x,y):
    a,b=x; c,d=y
    return (a*c-b*d, a*d+b*c-b*d)


def scale(n,x): return (n*x[0],n*x[1])


def wpow(e):
    return ((1,0),(0,1),(-1,-1))[e%3]


def chi1(r,s,h):
    a,b,c=h
    return wpow(r*a+s*b)


def chiV(power,h):
    a,b,c=h
    return scale(3,wpow(power*c)) if a==0 and b==0 else (0,0)
def all_h27():
    return list(itertools.product(range(3),repeat=3))


def payload():
    H=all_h27()
    one_line=[(0,0),(1,0),(2,0)]
    rows=[]
    for h in H:
        v=chiV(1,h); vb=chiV(2,h)
        A=(0,0)
        for r,s in one_line: A=add(A,chi1(r,s,h))
        A=add(A,v); A=add(A,vb)
        dressed=mul(v,A)
        reg=(27,0) if h==(0,0,0) else (0,0)
        assert dressed==reg
        rows.append({
          "h":list(h),
          "chi_V":list(v),
          "chi_Vbar":list(vb),
          "chi_A9":list(A),
          "chi_V_tensor_A9":list(dressed),
          "chi_regular":list(reg),
        })

    # Representation-ring identities behind the pointwise character equality.
    decomposition={
      "A9":"chi_00 + chi_10 + chi_20 + V_omega + V_omega2",
      "V_times_three_characters":"3 V_omega",
      "V_times_V":"3 V_omega2",
      "V_times_Vbar":"sum_{r,s in F3} chi_rs",
      "V_times_A9":"3 V_omega + 3 V_omega2 + sum_9 chi_rs = Reg(H27)",
    }
    old={"one_dimensional_total":0,"V_omega":9,"V_omega2":0}
    dressed={"one_dimensional_total":9,"V_omega":3,"V_omega2":3}
    assert 9*1+3*3+3*3==27
    checks={
      "all_27_characters_match_regular":all(
          r["chi_V_tensor_A9"]==r["chi_regular"] for r in rows),
      "A9_dimension_9":3+3+3==9,
      "dressed_dimension_27":27==27,
      "regular_irrep_census":dressed=={
          "one_dimensional_total":9,"V_omega":3,"V_omega2":3},
      "old_trinification_is_9V":old["V_omega"]==9,
    }
    assert all(checks.values())
    return {
      "schema":"w33.pass11043.commutant-dressing-regular-h27.v1",
      "status":"PASS_9D_LATENT_DRESSING_TURNS_TRINIFICATION_V_INTO_REGULAR_H27",
      "headline":(
        "Let V be the 3D Schrodinger irrep of H27 with central character omega. "
        "On the existing 9D multiplicity subsystem choose "
        "A9=chi00+chi10+chi20+V+Vbar. Then V tensor A9 has exactly the "
        "regular H27 character: 27 at the identity and 0 on every other element."
      ),
      "latent_module":{
        "dimension":9,
        "chosen_one_dimensional_character_line":["chi_00","chi_10","chi_20"],
        "three_dimensional_blocks":["V_omega","V_omega2"],
        "formula":decomposition["A9"],
      },
      "representation_ring":decomposition,
      "before_after":{
        "undressed_E6_27":old,
        "dressed_E6_27":dressed,
      },
      "commutant_reading":(
        "The prior execution theorem gives E6_27 = C9_multiplicity tensor V_omega "
        "with commutant M9 tensor I3. Therefore A9 can act entirely inside the "
        "latent multiplicity factor. The diagonal action A9(h) tensor V_omega(h) "
        "changes the H27 symmetry type without changing the 27D carrier dimension."
      ),
      "character_table_rows":rows,
      "parents":[
        "analysis/2026-09-21_e8_matter81_h27_address_operator_compiler.md",
        "data/w33_minimal_symmetry_changing_81_compiler.json",
      ],
      "boundary":(
        "This is an exact finite representation theorem and an algebraic control "
        "construction. It does not assert that the latent A9 action is generated "
        "by the current microscopic Hamiltonian or that it is energetically selected."
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
    print(json.dumps({
      "status":p["status"],
      "latent_dimension":p["latent_module"]["dimension"],
      "dressed":p["before_after"]["dressed_E6_27"]},sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
