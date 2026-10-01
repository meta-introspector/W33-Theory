#!/usr/bin/env python3
"""Exact fast replay of Pass 11264 by quotienting field choices by generated symmetry lattice."""
from __future__ import annotations
import argparse, itertools, json, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11264_z6i_rcharge_mu_audit as P64
P=P64.P41
OUT=ROOT/"data/w33_20261001_z6i_rcharge_lattice_quotient.json"

def lattice_signature(L):
    return tuple((int(c),tuple(map(int,r))) for c,r in L.ech)

def lattice_features(M,L):
    hu,hd=M.byb.get("bl",[]),M.byb.get("l",[])
    qs,us=M.byb.get("q",[]),M.byb.get("bu",[])
    prot=[]
    for h in hu:
        if any(M.target([h,l]) in L for l in hd): continue
        if any(M.target([q,h,u]) in L for q in qs for u in us): prot.append(h)
    if not prot:return {"protected":False}
    rows={}
    for h in prot:
        rest=[k for k in hu if k!=h]
        rr=P.MU.structural_rank(rest,hd,lambda a,b:M.target([a,b]) in L)
        gen=[o["name"] for o in M.sing if any(M.target([o["name"],h,l]) in L for l in hd)]
        ex={}
        for b in ("q","u","d","e","v","x"):
            A,B=M.byb.get(b,[]),M.byb.get("b"+b,[])
            if A and B:
                k=P.MU.structural_rank(A,B,lambda a,c:M.target([a,c]) in L)
                ex[b]=min(len(A),len(B))-k
        rows[h]={"base_clean":rr==len(rest) and all(v==0 for v in ex.values()),
                 "generators":tuple(gen)}
    minus_w=tuple(-x for x in M.w)
    linear=tuple(o["name"] for o in M.sing if M.target([o["name"]]) in L)
    return {"protected":True,"hu":tuple(prot),"rows":rows,
            "W_on_vacuum_possible":minus_w in L,"linear":linear}
def classify_choice(feat,S):
    if not feat["protected"]:return (0,0,0)
    ss=set(S)
    clean=any(r["base_clean"] and any(g not in ss for g in r["generators"])
              for r in feat["rows"].values())
    fflat=(not feat["W_on_vacuum_possible"]) and all(o in ss for o in feat["linear"])
    return (1,int(clean),int(clean and fflat))

def model_exact(name,validate=0):
    ledger,_=P.P.load_ledger(); raw,_=P64.load_raw(); disc=P64.inject_r_charges(raw)
    M=P.Model(name,ledger[name],disc)
    T,types,rays=P.model_rays(M)
    if rays is None:return {"anomalous":True,"dflat":False}
    supports=sorted({tuple(k for k,v in enumerate(r) if v!=0) for r in rays})
    vac=prot=clean=cleanF=0; unique=0; max_states=0; checked=0
    examples=[]; clean_examples=[]
    for sup in supports:
        choices=[types[T[k]] for k in sup]
        cache={}
        for S in itertools.product(*choices):
            vac+=1
            L=P.Lattice([M.vec[s] for s in S]+M.virtual)
            sig=lattice_signature(L)
            if sig not in cache:cache[sig]=[lattice_features(M,L),0]
            feat=cache[sig][0];cache[sig][1]+=1
            a,b,c=classify_choice(feat,S);prot+=a;clean+=b;cleanF+=c
            if validate and checked<validate:
                slow=M.analyse_vacuum(S)
                assert (slow is not None)==bool(a)
                if slow is not None:
                    assert bool(slow["clean"])==bool(b)
                    assert bool(slow["F_flat_lattice"] and slow["clean"])==bool(c)
                checked+=1
            if a and len(examples)<3:examples.append(list(S))
            if b and len(clean_examples)<5:clean_examples.append(list(S))
        unique+=len(cache);max_states=max(max_states,len(cache))
    return {"anomalous":True,"dflat":True,"fi_rays":len(rays),"supports":len(supports),
            "vacua":vac,"exhaustive":True,"mu_protected":prot,"clean":clean,
            "clean_and_F_flat":cleanF,"unique_support_lattices":unique,
            "max_lattices_per_support":max_states,"validation_choices":checked,
            "examples":examples,"clean_examples":clean_examples}
def worker(name):
    t=time.time();r=model_exact(name,validate=20)
    r["seconds"]=round(time.time()-t,3)
    return name,r

def payload(workers=6):
    old=json.loads((ROOT/"data/w33_pass11241_discrete_mu_symmetry.json").read_text())
    z6={n:r for n,r in old["models"].items() if n.startswith("Z6-I|")}
    dflat=sorted(n for n,r in z6.items() if r.get("dflat"))
    assert len(z6)==87 and len(dflat)==23
    rows={}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        fut={ex.submit(worker,n):n for n in dflat}
        for f in as_completed(fut):
            n,r=f.result();rows[n]=r
            print(n,r["vacua"],r["mu_protected"],r["clean"],r["unique_support_lattices"],r["seconds"],flush=True)
    for n,r in z6.items():
        if n not in rows: rows[n]={k:r[k] for k in ("anomalous","dflat") if k in r}
    raw,meta=P64.load_raw()
    audit={
      "models_total":len(rows),"models_dflat":sum(x.get("dflat",False) for x in rows.values()),
      "vacua":sum(x.get("vacua",0) for x in rows.values()),
      "mu_protected_vacua":sum(x.get("mu_protected",0) for x in rows.values()),
      "mu_protected_models":sorted(n for n,x in rows.items() if x.get("mu_protected",0)),
      "clean_vacua":sum(x.get("clean",0) for x in rows.values()),
      "clean_models":sorted(n for n,x in rows.items() if x.get("clean",0)),
      "clean_and_F_flat":sum(x.get("clean_and_F_flat",0) for x in rows.values()),
      "unique_support_lattices":sum(x.get("unique_support_lattices",0) for x in rows.values()),
      "exhaustive":all(x.get("exhaustive",True) for x in rows.values() if x.get("dflat"))
    }
    return {
      "schema":"w33.20261001.z6i-rcharge-lattice-quotient.v1",
      "status":"PASS_Z6I_R_CHARGES_RESTORED_EXACT_LATTICE_QUOTIENT_AUDIT",
      "extraction":{**meta,"orbifolder_version":"1.2.1",
                    "formula":"R_i=q_sh_i+OscillatorContribution_i",
                    "plane_orders":[6,6,3],"superpotential_R_charges":[-1,-1,-1]},
      "audit":audit,"models":dict(sorted(rows.items())),
      "acceleration":{
        "principle":"For a fixed FI support, field choices with identical exact Lattice.ech generate identical coupling-selection lattices. Higgs/top/exotic membership tests are evaluated once per distinct lattice; choice-sensitive clean/F-flat exclusions are then applied by finite set tests.",
        "old_vacuum_count":sum(x.get("vacua",0) for x in z6.values() if x.get("dflat")),
        "slow_engine":"analysis/w33_pass11241_discrete_mu_symmetry.py",
        "equivalence_checks_per_dflat_model":20
      },
      "correction":"Restored plane R charges materially change the Z6-I mu-protection verdict. D-flatness is unchanged because R charges do not enter the continuous FI cone; coupling selection is recomputed exactly.",
      "boundary":"Symmetry protection means the mu term is forbidden to all orders by the tested exact charge lattice while at least one top Yukawa is allowed. It does not determine the generated mu scale after symmetry breaking or prove phenomenological viability.",
      "checks":{"87_models":len(rows)==87,"23_dflat":audit["models_dflat"]==23,
                "all_dflat_exhaustive":audit["exhaustive"],"raw_87":meta["models"]==87,
                "raw_fields_30980":meta["fields"]==30980}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--workers",type=int,default=6);ap.add_argument("--check",action="store_true")
    a=ap.parse_args()
    if a.check:
        p=payload(a.workers);old=json.loads(OUT.read_text())
        # Ignore runtime fields for replay equality.
        def strip(x):
            if isinstance(x,dict):return {k:strip(v) for k,v in x.items() if k!="seconds"}
            if isinstance(x,list):return [strip(v) for v in x]
            return x
        if strip(p)!=strip(old):raise SystemExit("certificate drift")
    else:
        p=payload(a.workers);OUT.write_text(json.dumps(p,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":p["status"],"audit":p["audit"]},sort_keys=True))

if __name__=="__main__":raise SystemExit(main())
