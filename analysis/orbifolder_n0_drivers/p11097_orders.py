"""Pass 11097: exact minimal order of every allowed fractional-fermion mass term, and the smallest order T* at which
ALL fractional fermions can be paired (bottleneck matching), for the models where the symmetry test allows it.
usage: p11097_orders.py our_v2.dump their_sm.dump frac.json out.json"""
import json
import sys
from fractions import Fraction as F
from multiprocessing import Pool

import networkx as nx

sys.path.insert(0, sys.path[0])
sys.path.insert(0, r'C:\Repos\Theory of Everything\analysis')
import w33_exact_monomial_orders as EXO  # noqa: E402
from p11097_frac import conj  # noqa: E402
from so16_couplings import charge_vector, dimof, parse_ours, parse_theirs, sm_data, sm_info  # noqa: E402

KEY = 'A_hidden_unbroken|gauge+Z2W+PG+SG+R'
US = TH = None


def mod1(v, nu1):
    """continuous U(1) part, then the discrete coordinates as fractions mod 1: k/2, l/3, classes/3, R/3"""
    per = [2, 3, 3, 3, 3, 3, 3, 3]
    return list(v[:nu1]) + [F(x) / p for x, p in zip(v[nu1:], per)]


def job(i):
    us, th = US[i], TH[i]
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
    S = [s for s in scal if dimof(s['col']) == 1 and s['w'] == 1 and s['Y'] == 0 and s['hidden_singlet']]
    from fast_order import FastOrders
    E = FastOrders([mod1(charge_vector(s, True), nu1) for s in S], nu1, [1] * 8)
    frac = [f for f in ferm if f['frac']]
    w0 = [F(0)] * (nu1 + 8)
    edges, timeouts = [], []
    pairs = []
    for a_i, a in enumerate(frac):
        for b_i in range(a_i, len(frac)):
            b = frac[b_i]
            if b['col'] != conj(a['col'], 'A2') or b['w'] != a['w'] or a['Y'] + b['Y'] != 0:
                continue
            if any(b['dim'][h] != conj(a['dim'][h], algs[h]) for h in hid):
                continue
            pairs.append((a_i, b_i))
    SV = [charge_vector(s, True) for s in S]
    cubic = set()
    for a_i, b_i in pairs:
        va, vb = charge_vector(frac[a_i], True), charge_vector(frac[b_i], True)
        for vs in SV:
            tot = [x + y + z for x, y, z in zip(va, vb, vs)]
            if all(x == 0 for x in tot[:nu1]) and tot[nu1] % 2 == 0 and all(x % 3 == 0 for x in tot[nu1 + 1:nu1 + 5])                     and all(x == 0 for x in tot[nu1 + 5:nu1 + 8]):
                cubic.add((a_i, b_i))
                edges.append((a_i, b_i, 1))
                break
    total0 = sum(f['mult'] for f in frac)

    def light_edges(E_, T):
        G = nx.Graph()
        G.add_nodes_from(range(len(frac)))
        for a_i, b_i, o in E_:
            if o <= T:
                if a_i == b_i:
                    G.add_edge(a_i, ('self', a_i), weight=frac[a_i]['mult'])
                else:
                    G.add_edge(a_i, b_i, weight=2 * frac[a_i]['mult'])
        M = nx.max_weight_matching(G, weight='weight')
        massive = {u for e in M for u in e}
        return sum(frac[k]['mult'] for k in range(len(frac)) if k not in massive)
    Tstar, profile = None, {}
    lt = light_edges(edges, 1)
    profile[1] = lt
    if lt == 0:
        Tstar = 1
    else:
        rest = [(a_i, b_i) for a_i, b_i in pairs if (a_i, b_i) not in cubic]
        known = {}
        for T in range(2, 9):
            cache = {}
            for a_i, b_i in rest:
                if (a_i, b_i) in known:
                    continue
                va, vb = mod1(charge_vector(frac[a_i], True), nu1), mod1(charge_vector(frac[b_i], True), nu1)
                key = tuple(-(x + y) for x, y in zip(va, vb))
                if key not in cache:
                    cache[key] = E.feasible(list(key), T)
                if cache[key] is True:
                    known[(a_i, b_i)] = T
                    edges.append((a_i, b_i, T))
                elif cache[key] == 'timeout':
                    timeouts.append((a_i, b_i, T))
            lt = light_edges(edges, T)
            profile[T] = lt
            if lt == 0:
                Tstar = T
                break
    return i, dict(label=us['label'], vev_scalars=len(S), fractional_states=sum(f['mult'] for f in frac),
                   cubic_edges=len(cubic), T_star=Tstar, timeouts=len(timeouts), light_profile={str(k): v for k, v in profile.items()})


def init(ours, theirs):
    global US, TH
    US = parse_ours(ours)
    TH = parse_theirs(theirs)


def main(ours, theirs, fracp, outp):
    fr = json.load(open(fracp))
    todo = sorted(int(i) for i, v in fr.items() if v['results'][KEY]['light_fractional_states'] == 0)
    print('models:', len(todo), flush=True)
    with Pool(8, initializer=init, initargs=(ours, theirs)) as pool:
        res = {}
        for i, r in pool.imap_unordered(job, todo):
            res[i] = r
            print(i, r['label'], 'T* =', r['T_star'], 'timeouts', r['timeouts'], 'profile', r['light_profile'], flush=True)
    json.dump(res, open(outp, 'w'), indent=1, default=str)


if __name__ == '__main__':
    main(*sys.argv[1:5])
