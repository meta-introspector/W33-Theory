#!/usr/bin/env python3
"""Pass 11066: K81 -> U81 replaces a degree-1 external C3 by a degree-3 highest-root clock."""
from __future__ import annotations
import argparse,itertools,json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11066_jennings_clock_promotion.json"

class Span:
    def __init__(self,p,n):self.p=p;self.n=n;self.rows={}
    def add(self,v):
        p=self.p;v=[int(x)%p for x in v]
        while True:
            k=next((i for i,x in enumerate(v) if x),None)
            if k is None:return False
            if k not in self.rows:break
            f=v[k];r=self.rows[k];v=[(v[i]-f*r[i])%p for i in range(self.n)]
        z=pow(v[k],-1,p);v=[z*x%p for x in v]
        for j,r in list(self.rows.items()):
            if r[k]:
                f=r[k];self.rows[j]=[(r[i]-f*v[i])%p for i in range(self.n)]
        self.rows[k]=v;return True
    @property
    def rank(self):return len(self.rows)
    def basis(self):return list(self.rows.values())

def kmul(x,y):
    (a,b,c),p=x;(d,e,f),q=y
    return (((a+d)%3,(b+e)%3,(c+f-d*b)%3),(p+q)%3)

def augmentation_dims():
    G=[((a,b,c),p) for a,b,c,p in itertools.product(range(3),repeat=4)]
    idx={g:i for i,g in enumerate(G)}
    gens=[((1,0,0),0),((0,1,0),0),((0,0,0),1)]
    perms=[[idx[kmul(a,g)] for a in G] for g in gens]
    n=81;cur=Span(3,n)
    for perm in perms:
        for i,j in enumerate(perm):
            v=[0]*n;v[j]=1;v[i]=(v[i]-1)%3;cur.add(v)
    dims=[n,cur.rank]
    while cur.rank:
        nxt=Span(3,n)
        for v in cur.basis():
            for perm in perms:
                out=[0]*n
                for i,c in enumerate(v):
                    if c:
                        out[perm[i]]=(out[perm[i]]+c)%3
                        out[i]=(out[i]-c)%3
                nxt.add(out)
        cur=nxt;dims.append(cur.rank)
    return dims,[dims[i]-dims[i+1] for i in range(len(dims)-1)]

def payload():
    kd,kl=augmentation_dims()
    assert kd==[81,80,77,70,60,47,34,21,11,4,1,0]
    assert kl==[1,3,7,10,13,13,13,10,7,3,1]
    old=json.loads((ROOT/"data/PART_W33_PASS5108_U81_JENNINGS_MEMORY.json").read_text())
    ud=old["U81"]["augmentation_power_dimensions"];ul=old["U81"]["successive_layers"]
    assert ud==[81,80,78,74,69,62,54,45,36,27,19,12,7,3,1,0]
    assert ul==[1,2,4,5,7,8,9,9,9,8,7,5,4,2,1]

    return {
      "schema":"w33.pass11066.jennings-clock-promotion.v1",
      "status":"PASS_SPLIT_K81_TO_CHAMBER_U81_PROMOTES_ONE_DEGREE1_TERNARY_COORDINATE_TO_ROOT_HEIGHT3",
      "headline":"The split scheduler and the class-three chamber group have the same order 81 but different modular memory depth. Over F3, K81=H27 x C3 has Jennings Hilbert series (1+t+t^2)^3(1+t^2+t^4), with degrees 1,1,1,2. U81 has (1+t+t^2)^2(1+t^2+t^4)(1+t^3+t^6), with C2 root heights 1,1,2,3. Passing from K81 to U81 exactly replaces one degree-1 ternary factor by the degree-3 highest-root factor.",
      "K81":{
        "augmentation_power_dimensions":kd,"successive_layers":kl,
        "Hilbert_series":"(1+t+t^2)^3(1+t^2+t^4)",
        "Jennings_degrees":[1,1,1,2],"last_nonzero_augmentation_power":10,"augmentation_nilpotency_index":11
      },
      "U81":{
        "augmentation_power_dimensions":ud,"successive_layers":ul,
        "Hilbert_series":"(1+t+t^2)^2(1+t^2+t^4)(1+t^3+t^6)",
        "Jennings_degrees":[1,1,2,3],"last_nonzero_augmentation_power":14,"augmentation_nilpotency_index":15
      },
      "exact_promotion":{
        "removed_factor":"(1+t+t^2), degree 1",
        "inserted_factor":"(1+t^3+t^6), degree 3",
        "dimension_preserved":81,
        "interpretation":"the independent external C3 scheduler coordinate is replaced algebraically by the class-three highest-root coordinate"
      },
      "boundary":"Jennings degree is algebraic nilpotence depth, not laboratory time. The 'clock promotion' phrase names the exact filtration replacement and does not by itself prove a physical arrow of time.",
      "parents":["data/w33_pass11065_split_nonsplit_order81_bridge.json","data/PART_W33_PASS5108_U81_JENNINGS_MEMORY.json"],
      "checks":{"K_profile_exact":True,"U_profile_reused_exact":True,"degree1_to_degree3_replacement":True,"dimension81_preserved":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args()
    p=payload();text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"K_last":10,"U_last":14},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
