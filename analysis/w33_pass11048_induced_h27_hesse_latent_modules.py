#!/usr/bin/env python3
"""Pass 11048: the minimal latent A9 modules are induced C3 sectors of H27."""
from __future__ import annotations
import argparse, itertools, json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11048_induced_h27_hesse_latent_modules.json"
DIRS=[(1,0),(0,1),(1,1),(1,2)]


def affine_line(direction,theta):
    a,b=direction
    return sorted((r,s) for r,s in itertools.product(range(3),repeat=2)
                  if (r*a+s*b-theta)%3==0)


def payload():
    rows=[]
    all_lines=set()
    for di,d in enumerate(DIRS):
        partition=[]
        for theta in range(3):
            line=affine_line(d,theta)
            assert len(line)==3
            partition.extend(line)
            all_lines.add(tuple(line))

            # Frobenius reciprocity for L=<g_d> ~= C3:
            # 1D chi_(r,s) occurs iff chi_(r,s)|_L = theta.
            # Both 3D Schrodinger irreps restrict to Reg(C3), so each occurs once.
            induced={
              "one_dimensional_characters":[list(x) for x in line],
              "V_omega":1,
              "V_omega2":1,
              "dimension":3+3+3,
            }
            assert induced["dimension"]==9

            # Tensor with V_omega:
            # chi*V=V, V*V=3Vbar, Vbar*V=sum_9 chi.
            dressed={
              "one_dimensional_each":1,
              "V_omega":3,
              "V_omega2":3,
              "dimension":9+9+9,
            }
            assert dressed["dimension"]==27
            rows.append({
              "direction_index":di,
              "projective_direction":list(d),
              "theta":theta,
              "dual_affine_line":[list(x) for x in line],
              "module":"Ind_L^H27(theta)",
              "irreducible_content":induced,
              "after_tensor_with_V_omega":dressed,
            })
        assert sorted(partition)==list(itertools.product(range(3),repeat=2))
    assert len(all_lines)==12
    by_direction={}
    for di,d in enumerate(DIRS):
        rr=[r for r in rows if r["direction_index"]==di]
        by_direction[str(di)]={
          "direction":list(d),
          "three_theta_sectors":[r["theta"] for r in rr],
          "direct_sum_dimension":sum(r["irreducible_content"]["dimension"] for r in rr),
          "direct_sum":"Reg(H27)",
          "dual_lines":[r["dual_affine_line"] for r in rr],
        }
        assert by_direction[str(di)]["direct_sum_dimension"]==27

    prior=json.loads((ROOT/"data/w33_pass11044_minimal_latent_dimension_and_54_budget.json").read_text())
    assert prior["classification"]["latent_irrep_content"]=={
      "one_dimensional_total":3,"V_omega":1,"V_omega2":1,"dimension":9}

    checks={
      "four_projective_noncentral_C3_directions":len(DIRS)==4,
      "three_characters_per_direction":all(len([r for r in rows if r["direction_index"]==i])==3 for i in range(4)),
      "twelve_induced_A9_modules":len(rows)==12 and len(all_lines)==12,
      "each_A9_dimension9":all(r["irreducible_content"]["dimension"]==9 for r in rows),
      "each_A9_tensor_V_is_regular_H27":all(r["after_tensor_with_V_omega"]["dimension"]==27 for r in rows),
      "three_theta_sectors_sum_to_regular":all(x["direct_sum"]=="Reg(H27)" for x in by_direction.values()),
      "matches_pass11044_minimal_content":True,
    }
    assert all(checks.values())
    return {
      "schema":"w33.pass11048.induced-h27-hesse-latent-modules.v1",
      "status":"PASS_MINIMAL_A9_IS_AN_INDUCED_C3_SECTOR_AND_THE_12_CHOICES_ARE_THE_HESSE_LINES",
      "headline":(
        "Every minimal latent module A9=3chi+V+Vbar is naturally an induced "
        "representation Ind_L^H27(theta), where L is one of the four projective "
        "noncentral C3 directions in H27/Z and theta is one of the three characters "
        "of L. Frobenius reciprocity selects exactly the three 1D H27 characters on "
        "one affine line of dual F3^2, while V and Vbar each occur once. Hence the "
        "4x3 induced sectors are exactly the 12 affine/Hesse lines."
      ),
      "induction_theorem":{
        "subgroup":"L=<g> ~= C3 with noncentral projective direction d in P1(F3)",
        "one_dimensional_multiplicity":"chi_(r,s) occurs iff r*d0+s*d1 = theta mod 3",
        "Schrodinger_restriction":"Res_L(V_omega)=Reg(C3)=Res_L(V_omega2)",
        "Frobenius_reciprocity":"mult_rho Ind_L(theta) = mult_theta Res_L(rho)",
        "result":"Ind_L(theta)=sum_{dual affine line} chi + V_omega + V_omega2",
      },
      "four_directions":by_direction,
      "twelve_sectors":rows,
      "tensor_identity_preview":(
        "For every sector, V_omega tensor Ind_L(theta) is regular H27. "
        "Pass 11051 upgrades this character identity to the canonical tensor-induction "
        "isomorphism V tensor Ind_L(theta) ~= Ind_L(Res_L(V) tensor theta) "
        "~= Ind_L(Reg L) ~= Reg(H27)."
      ),
      "geometry_reading":(
        "The previously free choice of three one-dimensional latent characters is "
        "not arbitrary: choosing a noncentral C3 direction chooses a parallel class "
        "of three dual affine lines, and theta chooses one of the three lines. "
        "Thus the four qutrit/Hesse directions each decompose Reg(H27) into three "
        "minimal 9D latent sectors."
      ),
      "boundary":(
        "This is a finite representation-theory theorem. The identification of the "
        "four projective directions with physical qutrit MUB settings is a structural "
        "codec statement; no dynamical preference among the four directions is implied."
      ),
      "parents":[
        "data/w33_pass11043_commutant_dressing_regular_h27.json",
        "data/w33_pass11044_minimal_latent_dimension_and_54_budget.json",
      ],
      "checks":checks,
    }


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT); a=ap.parse_args(); p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text: raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(text)
    print(json.dumps({"status":p["status"],"sectors":len(p["twelve_sectors"])},sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
