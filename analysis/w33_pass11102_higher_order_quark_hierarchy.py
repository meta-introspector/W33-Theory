#!/usr/bin/env python3
"""Pass 11102: does the scalar vacuum split charm from up in the 12 surviving SO(16)xSO(16) A8 models?

Pass 11101: in the 12 survivors a localized light Higgs H gives one O(1) top; charm and up sit at the other two
theta-fixed points of the same SU(3) torus, equidistant from H (the three fixed points of a Z3 torus are an affine line
of AG(3,3); the point reflection about the top's fixed point swaps the other two), so m_c = m_u at tree level.

Here every up-type (and down-type) Yukawa entry q_i ubar_j H is extended to all orders with the SM-neutral hidden-singlet
scalars phi that give the fractional fermions their masses (Pass 11097): the minimal order n_ij of q_i ubar_j H phi^n
(n = 0: cubic, EXACT rules; n >= 1: H-momentum mod 3; integer programs of w33_so16_selection_rules).  Cubic entries with
different fixed points carry the geometric factor eps_geo^(tori differing).  Singular-value exponents are read off with
eps_VEV = eps and eps_geo = eps^g (g = 1/2, 1, 2) at eps = 1e-3, 1e-6 (random O(1) coefficients; exponents only).
The question: are the exponents of the second and third up-type singular values DIFFERENT (a charm/up hierarchy)?
Needs the full field dumps (WSL); results frozen in data/w33_pass11102_higher_order_quark_hierarchy.json.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11097_fractional_fermion_masses as P97  # noqa: E402
import w33_pass11098_yukawa_textures as Y  # noqa: E402
import w33_so16_selection_rules as S  # noqa: E402

OUT = ROOT / "data" / "w33_pass11102_higher_order_quark_hierarchy.json"
BEST = [2, 10, 13, 14, 15, 35, 53, 57, 69, 77, 78, 102]


def mod1(v, nu1):
    per = [2, 3, 3, 3, 3, 3, 3, 3]
    return list(v[:nu1]) + [F(x) / p for x, p in zip(v[nu1:], per)]


def order_matrix(A, B, h, Sv, E, nu1, hid, algs):
    """n_ij (None = forbidden at every order), and the geometric tori-differing vector for cubic entries"""
    n, geo = {}, {}
    w0 = [F(0)] * (nu1 + 8)
    for i, a in enumerate(A):
        for j, b in enumerate(B):
            if not Y.hidden_ok([a, b, h], hid, algs):
                continue
            vs = [S.charge_vector(x, True) for x in (a, b, h)]
            if P97.cubic_ok(vs, nu1):
                if all(x['k'] == 0 and x['l'] == 0 for x in (a, b, h)):
                    if Y.eps(Y.plane(a), Y.plane(b), Y.plane(h)) == 0:
                        continue
                    geo[(i, j)] = 0
                else:
                    ca, cb, ch = S.classes(a['n']), S.classes(b['n']), S.classes(h['n'])
                    geo[(i, j)] = sum(1 for t in range(3) if not (ca[t] == cb[t] == ch[t]))
                n[(i, j)] = 0
                continue
            o = E.coupling([mod1(v, nu1) for v in vs], w0, min_total=1)
            if isinstance(o, int):
                n[(i, j)] = o
    return n, geo


def exponents(n, geo, n1, n2, g, rng):
    coef = {k: rng.standard_normal() + 1j * rng.standard_normal() for k in n}
    out = []
    for eps in (1e-3, 1e-6):
        M = np.zeros((n1, n2), dtype=complex)
        for (i, j), o in n.items():
            M[i, j] = coef[(i, j)] * eps ** o * (eps ** (g * geo[(i, j)]) if o == 0 else 1.0)
        out.append(np.sort(np.linalg.svd(M, compute_uv=False))[::-1])
    return [None if x < 1e-250 or y < 1e-250 else round(float(np.log(y / x) / np.log(1e-3)), 2) for x, y in zip(*out)]


def analyse(us, th, rng):
    algs, hid, nu1 = P97.prepare(us, th)
    c, oc, ow, qcol = S.sm_data(th, us)
    qb = P97.conj(qcol, 'A2')
    ferm = [f for f in us['fields'] if f['m'] == 2]
    scal = [f for f in us['fields'] if f['m'] == 6]
    Sv = [s for s in scal if S.dimof(s['col']) == 1 and s['w'] == 1 and s['Y'] == 0 and s['hidden_singlet']]
    E = S.FastOrders([mod1(S.charge_vector(s, True), nu1) for s in Sv], nu1, [1] * 8)
    pick = lambda L, col, w, Yv: [f for f in L if f['col'] == col and f['w'] == w and f['Y'] == Yv]
    Q, U, D = pick(ferm, qcol, 2, F(1, 6)), pick(ferm, qb, 1, F(-2, 3)), pick(ferm, qb, 1, F(1, 3))
    Hu, Hd = pick(scal, '1', 2, F(1, 2)), pick(scal, '1', 2, F(-1, 2))
    res = {}
    for name, A, B, H in (('up', Q, U, Hu), ('down', Q, D, Hd)):
        per = []
        for k, h in enumerate(H):
            n, geo = order_matrix(A, B, h, Sv, E, nu1, hid, algs)
            if not any(o == 0 for o in n.values()):
                continue                      # only Higgs doublets with a tree-level (cubic) coupling can be the light one
            per.append(dict(k=k, orders={f"{i},{j}": o for (i, j), o in sorted(n.items())},
                            geo={f"{i},{j}": v for (i, j), v in sorted(geo.items())},
                            exponents={str(g): exponents(n, geo, len(A), len(B), g, rng) for g in (0.5, 1.0, 2.0)}))
        res[name] = per
    return res


def reflection_symmetric(h):
    """is the 2x2 block of the two non-top generations symmetric under exchanging them (orders, and geometric factors of
    the cubic entries) -- the point reflection of the generation line about the top's fixed point?"""
    n = {tuple(map(int, k.split(','))): o for k, o in h['orders'].items()}
    g = {tuple(map(int, k.split(','))): o for k, o in h['geo'].items()}
    top = [k for k, o in n.items() if o == 0 and g.get(k) == 0]
    if len(top) != 1:
        return None
    i0, j0 = top[0]
    I = [i for i in range(3) if i != i0]
    J = [j for j in range(3) if j != j0]
    cost = lambda k: (n.get(k), g.get(k) if n.get(k) == 0 else None)
    return any(cost((I[0], J[a])) == cost((I[1], J[1 - a])) and cost((I[0], J[1 - a])) == cost((I[1], J[a])) for a in (0, 1))


def summarize():
    d = json.loads(OUT.read_text())
    from collections import Counter
    sym = Counter(str(reflection_symmetric(h)) for v in d['models'].values() for s in ('up', 'down') for h in v[s])
    pats = Counter(str(h['exponents']['1.0']) for v in d['models'].values() for s in ('up', 'down') for h in v[s])
    d['summary'] = dict(higgs_cases=sum(sym.values()), reflection_symmetric=sym.get('True', 0), exponent_patterns=dict(pats))
    OUT.write_text(json.dumps(d, indent=1, sort_keys=True, default=str))
    print(d['summary'])


def main(ours, theirs):
    rng = np.random.default_rng(11102)
    US, TH = S.parse_ours(ours), S.parse_theirs(theirs)
    out = {}
    for m in BEST:
        r = analyse(US[m], TH[m], rng)
        out[str(m)] = dict(label=US[m]['label'], **r)
        for h in r['up'][:2]:
            print(m, 'up H', h['k'], 'orders', h['orders'], 'exp', h['exponents'], flush=True)
    OUT.write_text(json.dumps(dict(pass_id=11102, models=out), indent=1, sort_keys=True, default=str))


if __name__ == "__main__":
    if len(sys.argv) > 2:
        main(*sys.argv[1:3])
    summarize()
