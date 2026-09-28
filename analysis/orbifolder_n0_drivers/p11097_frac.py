"""Pass 11097 driver: can the fractionally charged fermions of the 104 SO(16)xSO(16) A8 models get masses from
SM-neutral scalar VEVs?  usage: p11097_frac.py our_v2.dump their_sm.dump out.json"""
import json
import sys
from collections import Counter

import networkx as nx

sys.path.insert(0, sys.path[0])
from so16_couplings import (F, Lattice, charge_vector, conj_token, dimof, parse_ours, parse_theirs, periods,  # noqa: E402
                            sm_data, sm_info)

SELFCONJ = {'A1'}


def conj(tok, alg):
    if tok.endswith('adj') or dimof(tok) == 1:
        return tok
    if alg in ('A1', 'E7', 'E8') or (alg.startswith('D') and int(alg[1:]) % 2 == 0):
        return tok
    return tok[1:] if tok.startswith('-') else '-' + tok


def analyse(us, th):
    c, oc, ow, qcol = sm_data(th, us)
    algs = [a for a, _ in us['factors']]
    hid = [j for j in range(len(algs)) if j not in (oc, ow)]
    nu1 = len(us['u1'])
    ferm = [f for f in us['fields'] if f['m'] == 2]
    scal = [f for f in us['fields'] if f['m'] == 6]
    for f in ferm + scal:
        f['Y'], f['col'], f['w'], f['t'], f['frac'] = sm_info(f, c, oc, ow, qcol)
        f['hidden_singlet'] = all(dimof(f['dim'][j]) == 1 for j in hid)
        f['mult'] = 1
        for x in f['dim']:
            f['mult'] *= dimof(x)
    smneutral = [s for s in scal if dimof(s['col']) == 1 and s['w'] == 1 and s['Y'] == 0]
    vev = {'A_hidden_unbroken': [s for s in smneutral if s['hidden_singlet']], 'B_hidden_broken': smneutral}
    frac = [f for f in ferm if f['frac']]
    out = dict(fractional_fields=len(frac), fractional_states=sum(f['mult'] for f in frac),
               fractional_hidden_singlet_states=sum(f['mult'] for f in frac if f['hidden_singlet']),
               sm_neutral_scalars=len(smneutral), results={})
    for variant, S in vev.items():
        for rules in ('gauge+Z2W+PG+SG', 'gauge+Z2W+PG+SG+R'):
            useR = rules.endswith('+R')
            vecs = [charge_vector(s, useR) for s in S] + periods(nu1, useR)
            L = Lattice(vecs)
            G = nx.Graph()
            for i, a in enumerate(frac):
                G.add_node(i)
            for i, a in enumerate(frac):
                for j in range(i, len(frac)):
                    b = frac[j]
                    if b['col'] != conj(a['col'], 'A2') or b['w'] != a['w'] or a['Y'] + b['Y'] != 0:
                        continue
                    if variant.startswith('A'):
                        if any(b['dim'][h] != conj(a['dim'][h], algs[h]) for h in hid):
                            continue
                    va, vb = charge_vector(a, useR), charge_vector(b, useR)
                    if not L.contains([x + y for x, y in zip(va, vb)]):
                        continue
                    if i == j:
                        G.add_edge(i, ('self', i), weight=a['mult'])
                    else:
                        G.add_edge(i, j, weight=2 * a['mult'])
            M = nx.max_weight_matching(G, maxcardinality=False, weight='weight')
            massive = set()
            for u, v in M:
                massive.add(u)
                massive.add(v)
            light = [frac[i] for i in range(len(frac)) if i not in massive]
            out['results'][f'{variant}|{rules}'] = dict(
                vev_scalars=len(S), allowed_edges=G.number_of_edges(),
                light_fractional_states=sum(f['mult'] for f in light),
                light_fractional_hidden_singlet_states=sum(f['mult'] for f in light if f['hidden_singlet']),
                light_reps=dict(Counter(f"({f['col']},{f['w']})_{f['Y']}" + ('' if f['hidden_singlet'] else '*hidden') for f in light)))
    return out


def main(ours, theirs, outp):
    US = parse_ours(ours)
    TH = parse_theirs(theirs)
    res = {}
    for i in sorted(US):
        r = analyse(US[i], TH[i])
        res[i] = dict(label=US[i]['label'], **r)
        short = {k.split('|')[0][0] + ('R' if k.endswith('+R') else ''): (v['light_fractional_states'], v['light_fractional_hidden_singlet_states'])
                 for k, v in r['results'].items()}
        print(i, US[i]['label'], 'frac states', r['fractional_states'], '(hidden-singlet', r['fractional_hidden_singlet_states'], ')',
              'light (all, hidden-singlet):', short, flush=True)
    json.dump(res, open(outp, 'w'), indent=1, default=str)
    for key in next(iter(res.values()))['results']:
        vals = [v['results'][key]['light_fractional_hidden_singlet_states'] for v in res.values()]
        allv = [v['results'][key]['light_fractional_states'] for v in res.values()]
        print(key, 'models with 0 light hidden-singlet fractional states:', sum(1 for x in vals if x == 0),
              '| with 0 light fractional states at all:', sum(1 for x in allv if x == 0), '| min', min(vals), 'max', max(vals))


if __name__ == '__main__':
    main(*sys.argv[1:4])
