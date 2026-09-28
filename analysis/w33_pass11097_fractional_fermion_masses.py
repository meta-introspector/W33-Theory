#!/usr/bin/env python3
"""Pass 11097: the fractionally charged fermions of the SO(16)xSO(16) W(3,3) A8 models -- which can get masses?

Every one of the 104 tachyon-free three-generation models of Pass 11095 carries 88-230 vector-like fractionally
charged fermion states (Q not in Z for colour singlets, 3Q + t not 0 mod 3 for triplets).  In Z3xZ3 the analogous states
proved unremovable (Passes 11024, 11087-11090).  Here, with the selection rules of analysis/w33_so16_selection_rules.py
(gauge, Witten Z2, point group, space group, H-momentum):

  A. symmetry test (all orders, existence): with VEVs for the SM-neutral scalars -- A: hidden singlets only (hidden
     gauge group unbroken), B: also hidden-charged ones (hidden group broken, hidden quantum numbers ignored -- the
     most permissive) -- is a mass term psi_i psi_j x (monomial) allowed?  Maximum matching over the fractional
     fermions (Majorana-type self pairing for self-conjugate reps) gives the number that must stay massless;
  B. orders: for the models where everything can pair (A, all rules), the exact minimal order of each mass term
     (order 1 = cubic psi psi phi with EXACT H-momentum; higher orders with the mod-3 rule; integer programs
     validated 240/240 against w33_exact_monomial_orders) and the bottleneck order T* at which ALL fractional
     fermions can be massive simultaneously.

Frozen full-run results (104 models; 33 for B) in data/w33_pass11097_*.json; three sample models
(data/w33_pass11097_so16_a8_sample_models.json) are re-run from their field dumps here.
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

import networkx as nx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_so16_selection_rules as S  # noqa: E402

SAMPLE = ROOT / "data" / "w33_pass11097_so16_a8_sample_models.json"
FRAC = ROOT / "data" / "w33_pass11097_symmetry_test_104.json"
ORDERS = ROOT / "data" / "w33_pass11097_orders_33.json"
OUT = ROOT / "data" / "w33_pass11097_fractional_fermion_masses.json"
KEY = "A_hidden_unbroken|gauge+Z2W+PG+SG+R"


def conj(tok, alg):
    if tok.endswith('adj') or S.dimof(tok) == 1:
        return tok
    if alg in ('A1', 'E7', 'E8') or (alg.startswith('D') and int(alg[1:]) % 2 == 0):
        return tok
    return tok[1:] if tok.startswith('-') else '-' + tok


def prepare(us, th):
    c, oc, ow, qcol = S.sm_data(th, us)
    algs = [a for a, _ in us['factors']]
    hid = [j for j in range(len(algs)) if j not in (oc, ow)]
    for f in us['fields']:
        f['Y'], f['col'], f['w'], f['t'], f['frac'] = S.sm_info(f, c, oc, ow, qcol)
        f['hidden_singlet'] = all(S.dimof(f['dim'][j]) == 1 for j in hid)
        f['mult'] = 1
        for x in f['dim']:
            f['mult'] *= S.dimof(x)
    return algs, hid, len(us['u1'])


def pairs_for(frac, algs, hid, hidden_unbroken):
    out = []
    for i, a in enumerate(frac):
        for j in range(i, len(frac)):
            b = frac[j]
            if b['col'] != conj(a['col'], 'A2') or b['w'] != a['w'] or a['Y'] + b['Y'] != 0:
                continue
            if hidden_unbroken and any(b['dim'][h] != conj(a['dim'][h], algs[h]) for h in hid):
                continue
            out.append((i, j))
    return out


def light(frac, edges):
    G = nx.Graph()
    G.add_nodes_from(range(len(frac)))
    for i, j in edges:
        if i == j:
            G.add_edge(i, ('self', i), weight=frac[i]['mult'])
        else:
            G.add_edge(i, j, weight=2 * frac[i]['mult'])
    M = nx.max_weight_matching(G, weight='weight')
    massive = {u for e in M for u in e}
    rest = [frac[k] for k in range(len(frac)) if k not in massive]
    return sum(f['mult'] for f in rest), sum(f['mult'] for f in rest if f['hidden_singlet'])


def symmetry_test(us, th):
    algs, hid, nu1 = prepare(us, th)
    ferm = [f for f in us['fields'] if f['m'] == 2]
    scal = [f for f in us['fields'] if f['m'] == 6]
    frac = [f for f in ferm if f['frac']]
    smn = [s for s in scal if S.dimof(s['col']) == 1 and s['w'] == 1 and s['Y'] == 0]
    out = {}
    for variant, V in (('A_hidden_unbroken', [s for s in smn if s['hidden_singlet']]), ('B_hidden_broken', smn)):
        for rules in ('gauge+Z2W+PG+SG', 'gauge+Z2W+PG+SG+R'):
            useR = rules.endswith('+R')
            L = S.Lattice([S.charge_vector(s, useR) for s in V] + S.periods(nu1, useR))
            ok = [(i, j) for i, j in pairs_for(frac, algs, hid, variant.startswith('A'))
                  if L.contains([x + y for x, y in zip(S.charge_vector(frac[i], useR), S.charge_vector(frac[j], useR))])]
            lt, lh = light(frac, ok)
            out[f'{variant}|{rules}'] = dict(light_fractional_states=lt, light_fractional_hidden_singlet_states=lh)
    return dict(fractional_states=sum(f['mult'] for f in frac), results=out)


def cubic_ok(vs, nu1):
    tot = [sum(x) for x in zip(*vs)]
    return all(x == 0 for x in tot[:nu1]) and tot[nu1] % 2 == 0 and all(x % 3 == 0 for x in tot[nu1 + 1:nu1 + 5]) \
        and all(x == 0 for x in tot[nu1 + 5:nu1 + 8])


def cubic_T1(us, th):
    """is every fractional fermion massive already through CUBIC couplings (exact rules, hidden unbroken)?"""
    algs, hid, nu1 = prepare(us, th)
    frac = [f for f in us['fields'] if f['m'] == 2 and f['frac']]
    V = [S.charge_vector(s, True) for s in us['fields'] if s['m'] == 6 and S.dimof(s['col']) == 1 and s['w'] == 1
         and s['Y'] == 0 and s['hidden_singlet']]
    edges = [(i, j) for i, j in pairs_for(frac, algs, hid, True)
             if any(cubic_ok([S.charge_vector(frac[i], True), S.charge_vector(frac[j], True), v], nu1) for v in V)]
    return light(frac, edges)


def sample_models():
    d = json.loads(SAMPLE.read_text())
    with tempfile.TemporaryDirectory() as tmp:
        po, pt = os.path.join(tmp, 'o.dump'), os.path.join(tmp, 't.dump')
        open(po, 'w').write(''.join(d['ours'].values()))
        open(pt, 'w').write(''.join(d['theirs'].values()))
        US, TH = S.parse_ours(po), S.parse_theirs(pt)
    return US, TH


def main():
    US, TH = sample_models()
    frac = json.loads(FRAC.read_text())
    orders = json.loads(ORDERS.read_text())
    recheck = {}
    for i in US:
        st = symmetry_test(US[i], TH[i])
        recheck[str(i)] = dict(symmetry_test=st, frozen=frac[str(i)]['results'],
                               agrees=all(st['results'][k]['light_fractional_states'] == frac[str(i)]['results'][k]['light_fractional_states']
                                          for k in st['results']),
                               cubic_light=cubic_T1(US[i], TH[i]))
        print(i, recheck[str(i)]['agrees'], st['results'][KEY], 'cubic light', recheck[str(i)]['cubic_light'], flush=True)
    per = {k: v['results'] for k, v in frac.items()}
    summ = {}
    for key in next(iter(per.values())):
        summ[key] = dict(models_all_fractional_massive=sum(1 for v in per.values() if v[key]['light_fractional_states'] == 0),
                         models_no_light_hidden_singlet=sum(1 for v in per.values() if v[key]['light_fractional_hidden_singlet_states'] == 0))
    Ts = Counter(str(v['T_star']) for v in orders.values())
    res = dict(pass_id=11097, models=len(frac), fractional_states_range=[min(v['fractional_states'] for v in frac.values()),
                                                                         max(v['fractional_states'] for v in frac.values())],
               symmetry_test=summ, bottleneck_order_T_star=dict(Ts), orders_models=len(orders),
               sample_recheck=recheck)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True, default=str))
    print(json.dumps({k: v for k, v in res.items() if k != 'sample_recheck'}, indent=1, default=str))


if __name__ == "__main__":
    main()
