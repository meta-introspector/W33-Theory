#!/usr/bin/env python3
"""Pass 11084: PGSp extends the oriented-chamber double cover to a U81:V4 four-sheet cover."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11084_pgsp_u81_four_sheet_cover.json"
J=np.array([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]],dtype=np.int64)%3
I=np.eye(4,dtype=np.int64)%3
T=np.diag([1,2,1,2]).astype(np.int64)%3
S=np.diag([1,2,2,1]).astype(np.int64)%3
def E(i,j):
    M=np.zeros((4,4),dtype=np.int64);M[i,j]=1;return M
X=[(E(0,1)-E(3,2))%3,E(1,3)%3,(E(0,3)+E(1,2))%3,E(0,2)%3]
def signs(M):
    out=[]
    Minv=M
    for A in X:
        C=(M@A@Minv)%3
        out.append(next(z for z in (1,2) if np.array_equal(C,(z*A)%3)))
    return out
def payload():
    borel=json.loads((ROOT/"data/PART_W33_PASS4519_FLAG_BOREL_SYLOW3_NORMALIZER.json").read_text())
    internal=json.loads((ROOT/"data/w33_pass11072_internal_borel_deck_involution.json").read_text())
    outer=json.loads((ROOT/"data/w33_pass11073_outer_similitude_tick_reversal.json").read_text())
    assert borel["group_order"]==25920 and borel["normalizer"]["order"]==162 and borel["sylow3"]["order"]==81
    assert borel["flag"]["index"]==160
    assert np.array_equal((T.T@J@T)%3,J) and np.array_equal((S.T@J@S)%3,(2*J)%3)
    assert np.array_equal((T@T)%3,I) and np.array_equal((S@S)%3,I)
    assert np.array_equal((T@S)%3,(S@T)%3)
    assert signs(T)==[2,1,2,1] and signs(S)==[2,2,1,2]
    assert internal["checks"]["T_is_symplectic"] and outer["checks"]["S_is_minus_symplectic_similitude"]
    psp=25920;pgsp=2*psp;nsyl=160;normalizer=pgsp//nsyl
    assert normalizer==324 and normalizer//81==4 and pgsp//81==640
    return {
      "schema":"w33.pass11084.pgsp-u81-four-sheet-cover.v1",
      "status":"PASS_PGSP_NORMALIZER_OF_U81_IS_U81_SEMIDIRECT_V4_AND_G_OVER_U81_IS_A_FOUR_SHEET_FLAG_COVER",
      "headline":"All Sylow-3 subgroups of PGSp(4,3) lie in its index-two PSp subgroup, so the 160 chamber U81 subgroups are unchanged. Hence N_PGSp(U81) has order 51840/160=324. The internal symplectic deck involution T and the outer anti-symplectic similitude S both normalize the C2 positive-root group, commute, and generate a V4 complementary to the 3-group U81. Therefore N_PGSp(U81)/U81 ~= V4 and PGSp/U81 has 640 elements, four over each of 160 flags.",
      "orders":{"PSp":psp,"PGSp":pgsp,"U81":81,"sylow3_subgroups":nsyl,
                "normalizer_PSp":162,"normalizer_PGSp":normalizer,
                "normalizer_quotient_order":4,"PGSp_over_U81":640,"flags":160},
      "deck_generators":{"T":{"matrix":"diag(1,-1,1,-1)","location":"PSp","root_signs":["-","+","-","+"]},
                         "S":{"matrix":"diag(1,-1,-1,1)","location":"PGSp-PSp","root_signs":["-","-","+","-"]},
                         "relations":["T^2=1","S^2=1","TS=ST"],"quotient":"V4"},
      "sylow_argument":"A 3-subgroup maps trivially to PGSp/PSp~=C2, so every PGSp Sylow-3 subgroup is contained in PSp. Thus the PSp census of 160 Sylow-3 subgroups is also the PGSp census.",
      "boundary":"The four sheets are exact finite-group orientation data. Their physical interpretation as time/chirality choices requires a dynamical selection principle.",
      "parents":["data/PART_W33_PASS4519_FLAG_BOREL_SYLOW3_NORMALIZER.json","data/w33_pass11072_internal_borel_deck_involution.json","data/w33_pass11073_outer_similitude_tick_reversal.json"],
      "checks":{"PGSp_order51840":True,"same_160_sylow3":True,"normalizer_order324":True,
                "T_S_commuting_involutions":True,"root_group_normalized":True,"quotient_V4":True,"four_sheet_cover640_to160":True}
    }
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or json.loads(a.output.read_text())!=p:raise SystemExit("certificate drift")
    else:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"cover":"640 -> 160","deck":"V4"},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
