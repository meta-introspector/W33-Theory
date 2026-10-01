#!/usr/bin/env python3
"""Sensitivity probe, not a corrected physical selection-rule census.
Three stored examples in each of 23 D-flat models, under relaxed charge rules.
Requires the exact polyhedral dependency used by Pass 11241 (Windows py -3).
"""
import sys,json,copy,time
from pathlib import Path
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'analysis'))
import w33_20261001_z6i_rcharge_lattice_quotient as Q
P=Q.P;ledger,_=P.P.load_ledger();raw,_=Q.P64.load_raw();base=Q.P64.inject_r_charges(raw)
cert=json.loads((root/'data/w33_20261001_z6i_rcharge_lattice_quotient.json').read_text());rows=[]
for policy in ['third_plane_only','no_R']:
 d=copy.deepcopy(base)
 for n,v in d.items():
  if n.startswith('Z6I_'):
   v['R']=[v['R'][2]] if policy=='third_plane_only' else []
   for f in v['fields'].values():f['R']=[f['R'][2]] if policy=='third_plane_only' else []
 for name,r in cert['models'].items():
  if not r.get('dflat'):continue
  m=P.Model(name,ledger[name],d)
  for s in r['examples']:
   out=m.analyse_vacuum(s);row=dict(policy=policy,model=name,support=s,protected=out is not None,clean=bool(out and out['clean']))
   rows.append(row);print(json.dumps(row),flush=True)
(root/'data/w33_20261001_z6i_prime_plane_sensitivity.json').write_text(json.dumps(rows,indent=2))
