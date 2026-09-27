#!/usr/bin/env python3
"""Pass 11038: W33 spread census and complete-positivity instrument firewall."""
from __future__ import annotations
import argparse, itertools, json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11038_spread_instrument_cp_firewall.json"
Q=3


def canon(v):
    v=tuple(int(x)%Q for x in v)
    for x in v:
        if x:
            inv=pow(x,-1,Q)
            return tuple((inv*y)%Q for y in v)
    raise ValueError("zero")


def omega(x,y):
    return (x[0]*y[1]-x[1]*y[0]-x[2]*y[3]+x[3]*y[2])%Q


def geometry():
    pts=sorted({canon(v) for v in itertools.product(range(Q),repeat=4) if any(v)})
    pidx={p:i for i,p in enumerate(pts)}
    lines=set()
    for i,a in enumerate(pts):
        for b in pts[i+1:]:
            if omega(a,b): continue
            L=frozenset(canon(tuple((r*x+s*y)%Q for x,y in zip(a,b)))
                        for r,s in itertools.product(range(Q),repeat=2) if r or s)
            if len(L)==4: lines.add(L)
    lines=sorted(lines,key=lambda L:sorted(pidx[p] for p in L))
    assert len(pts)==len(lines)==40
    return pts,pidx,lines
def spreads(pts,pidx,lines):
    masks=[]
    through=[[] for _ in pts]
    for li,L in enumerate(lines):
        m=0
        for p in L:
            j=pidx[p]; m|=1<<j; through[j].append(li)
        masks.append(m)
    full=(1<<len(pts))-1
    out=[]
    def rec(mask,chosen):
        if mask==full:
            out.append(tuple(sorted(chosen))); return
        uncovered=[j for j in range(len(pts)) if not(mask>>j&1)]
        j=min(uncovered,key=lambda x:sum((masks[li]&mask)==0 for li in through[x]))
        for li in through[j]:
            if masks[li]&mask: continue
            rec(mask|masks[li],chosen+[li])
    rec(0,[])
    return sorted(set(out))


def payload():
    pts,pidx,lines=geometry()
    B=frozenset(canon((u[0],u[1],u[0],u[1]))
                for u in itertools.product(range(Q),repeat=2) if any(u))
    bi=lines.index(B)
    S=spreads(pts,pidx,lines)
    assert len(S)==36
    withB=[s for s in S if bi in s]
    assert len(withB)==9

    pair_intersections=Counter(
        len(set(a)&set(b)) for i,a in enumerate(S) for b in S[i+1:])
    assert pair_intersections==Counter({1:360,4:270})
    # The four Bell rays are exactly the projective two-sided Weyl basis maps
    # Phi_(u,u), hence the only basis rays whose rank-one Choi matrices are PSD.
    # Every spread covers all 40 points, so every spread contains exactly four
    # CP basis rays and 36 non-CP basis rays.
    spread_cp_profiles=Counter(
        tuple(sorted(len(lines[li]&B) for li in s)) for s in S)
    expected_profiles=Counter({
        (0,0,0,0,0,0,0,0,0,4):9,
        (0,0,0,0,0,0,1,1,1,1):27,
    })
    assert spread_cp_profiles==expected_profiles

    sectors={b:[p for p in pts if p not in B and omega(p,b)==0] for b in B}
    assert sorted(map(len,sectors.values()))==[9,9,9,9]
    transverse=[i for i,L in enumerate(lines) if not(L&B)]
    assert len(transverse)==27
    assert all(
        all(len(lines[li]&set(sectors[b]))==1 for b in B)
        for s in withB for li in s if li!=bi
    )
    incidence=Counter(li for s in withB for li in s if li!=bi)
    assert len(incidence)==27 and set(incidence.values())=={3}

    checks={
      "W33_40_points_40_lines":len(pts)==len(lines)==40,
      "spread_count_36":len(S)==36,
      "Bell_containing_spreads_9":len(withB)==9,
      "spread_pair_intersections_1_or_4":
          pair_intersections==Counter({1:360,4:270}),
      "each_spread_has_4_CP_basis_rays":all(
          sum(len(lines[li]&B) for li in s)==4 for s in S),
      "Bell_spread_profile_one_CP_context":spread_cp_profiles[
          (0,0,0,0,0,0,0,0,0,4)]==9,
      "nonBell_spread_profile_four_single_CP_contexts":spread_cp_profiles[
          (0,0,0,0,0,0,1,1,1,1)]==27,
      "Bell_relative_sectors_4_by_9":sorted(map(len,sectors.values()))==[9]*4,
      "Bell_spread_transverse_lines_hit_each_sector_once":True,
      "each_transverse_line_occurs_in_three_Bell_spreads":set(incidence.values())=={3},
    }
    assert all(checks.values())
    return {
      "schema":"w33.pass11038.spread-instrument-cp-firewall.v1",
      "status":"PASS",
      "headline":(
        "W33 has exactly 36 line spreads; 9 contain the Bell line. The spread "
        "intersection census is {1:360,4:270}. But a spread is not directly a "
        "ten-context family of physical quantum channels: among the 40 projective "
        "two-sided Weyl basis rays only the four Bell-diagonal rays are completely "
        "positive. Every spread therefore contains 4 CP basis rays and 36 non-CP "
        "basis rays."
      ),
      "spread_census":{
        "spreads":len(S),
        "lines_per_spread":10,
        "Bell_line_index":bi,
        "Bell_containing_spreads":len(withB),
        "not_containing_Bell":len(S)-len(withB),
        "pair_intersection_histogram":{str(k):v for k,v in sorted(pair_intersections.items())},
        "Bell_spread_example":list(withB[0]),
      },
      "complete_positivity":{
        "basis_map":"Phi_(u,v)(A)=P_u A P_v^dagger",
        "rank_one_Choi":"|P_u>><<P_v|",
        "PSD_condition":"u=v up to the fixed Pauli phase convention",
        "CP_projective_points":4,
        "CP_locus":"the Bell line",
        "non_CP_projective_points":36,
        "spread_CP_context_profiles":{
          "9_Bell_spreads":"one line with 4 CP rays; nine lines with 0",
          "27_other_spreads":"four lines with 1 CP ray; six lines with 0",
        },
      },
      "Bell_relative_geometry":{
        "off_Bell_points":36,
        "sectors":4,
        "points_per_sector":9,
        "transverse_lines":27,
        "Bell_spread_rule":"each of the other nine lines hits every Bell-relative sector once",
        "transverse_line_Bell_spread_incidence":3,
      },
      "physical_firewall":(
        "The uploaded temporal architecture is safe if a W33 point is called an "
        "operator/intervention direction. It is too strong to call all 40 basis "
        "rays physical channels: complete positivity selects only the Bell line. "
        "A ten-setting laboratory instrument built from a spread therefore needs "
        "an explicit CP completion (positive linear combinations, ancilla dilation, "
        "or a different measurement map) for the 36 off-Bell directions."
      ),
      "boundary":(
        "This does not say off-Bell W33 directions are unobservable; they form an "
        "exact operator basis and can enter tomography by linear reconstruction. "
        "It says only that the raw rank-one superoperators are not individually "
        "implementable CPTP interventions."
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
    print(json.dumps({"status":p["status"],"spreads":p["spread_census"]["spreads"],
      "Bell_spreads":p["spread_census"]["Bell_containing_spreads"],
      "CP_points":p["complete_positivity"]["CP_projective_points"]},sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
