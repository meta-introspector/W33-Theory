#!/usr/bin/env python3
"""Pass 11103: no alignment of the Standard-Model-neutral scalar VEVs splits charm from up in the 12 surviving
SO(16)xSO(16) A8 models -- the Delta(54) protection is not just the reflection.

Pass 11102 counted ORDERS with every scalar VEV of the same size; the point reflection that swaps charm and up (the -I of
the flavour group Delta(54) = H27 : <-I>, Pass 11105) then maps the whole VEV set to itself.  A vacuum can break it
spontaneously: the twisted scalars sit at definite fixed points of the family torus, and their VEVs need not be equal.
Here the VEV of every scalar is eps^w, with w depending on the GEOMETRIC fixed point g = l c of the scalar in the family
torus, relative to the point p0 of the light Higgs (theta^2 class labels are inverted; untwisted scalars w = 1):

    sym (1,1,1)   p2_small (1,1,2)   split12 (1,1.5,3)   only_p0 (1,-,-)   only_p1 (-,1,-)   p0p1 (1,1,-)
    [ (w_p0, w_p1, w_p2);  '-' = zero VEV ]

For every Yukawa entry q_i ubar_j H (and q_i dbar_j H_d) the minimal WEIGHTED order sum_s w_s e_s is an integer program
(HiGHS) under all selection rules (gauge, Witten Z2, point group, space group, H-momentum mod 3; cubic entries exact);
cubic entries carry the geometric factor eps_geo^(tori differing) with eps_geo = eps^g, g = 1 (small torus area) and
g = 4 (large area: VEV terms can dominate).  Singular-value exponents from eps = 1e-3, 1e-6 with random O(1)
coefficients.  Two VEV sets: A (hidden singlets only, hidden group unbroken, 120-126 scalars) and B (also the hidden-
charged SM-neutral scalars, hidden group broken, 126-138).  Results frozen in data/w33_pass11103_*.json.

Result (12 models x 3 light-Higgs cases, up and down, 6 alignments x 2 regimes x 2 VEV sets):
  * UP: charm and up are at the SAME order in 36/36 cases for EVERY alignment, both regimes, both VEV sets.
    The leading diagonal entries always need one scalar at the Higgs point p0 and one at a non-Higgs point, and are
    reached at the same weighted cost for charm and up (the theta and theta^2 scalars at a point carry opposite class
    charges); VEVs only off p0 kill the diagonal entries and leave the reflection-symmetric off-diagonal block.
  * DOWN: a parametric strange/down split (exponent gap 1 or 3) appears at large torus area for reflection-breaking
    alignments in exactly the 8 models with three extra dbar-like states (24/36 cases), never in the 4 minimal models:
    it runs through the mixing with the vector-like dbar states (the 3x6 Q.dbar block; the vector-like mass matrix is
    not included).
  * VEVs aligned at the Higgs point (only_p0) never split anything, as Schur's lemma demands (Pass 11105 D).
Reading: m_c/m_u ~ 600 cannot come from any VEV alignment of the SM-neutral scalars in these models; it needs an O(1)
coefficient tuning or physics outside the scalar vacuum (moduli-dependent Yukawas / non-perturbative effects).
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from scipy.optimize import LinearConstraint, milp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
OUT = ROOT / "data" / "w33_pass11103_reflection_breaking_vacua.json"
FROZEN = {"A": ROOT / "data" / "w33_pass11103_alignments_hidden_unbroken.json",
          "B": ROOT / "data" / "w33_pass11103_alignments_hidden_broken.json"}
GRID = [('sym', (1, 1, 1)), ('p2_small', (1, 1, 2)), ('split12', (1, 1.5, 3)), ('only_p0', (1, None, None)),
        ('only_p1', (None, 1, None)), ('p0p1', (1, 1, None))]
GS = (1.0, 4.0)


def geo_pt(f, t, classes):
    """geometric fixed point in torus t (theta^2 class labels are inverted)"""
    return (f['l'] * classes(f['n'])[t]) % 3


def weighted_orders_class():
    import w33_so16_selection_rules as S

    class Weighted(S.FastOrders):
        """minimal weighted order: sum_s w_s e_s (w_s = None: zero VEV)"""

        def __init__(self, vecs, ws, nq, periods, emax=60):
            best = {}
            for v, w in zip(vecs, ws):
                if w is not None:
                    k = tuple(F(x) for x in v)
                    best[k] = min(best.get(k, 1e9), w)
            S.FastOrders.__init__(self, [list(k) for k in best], nq, periods, emax)
            self.cobj = np.concatenate([np.array([best[tuple(v)] for v in self.vecs], float), np.zeros(self.nd)])

        def order(self, target, min_total=0):
            t = [F(x) * self.den for x in target]
            if any(x.denominator != 1 for x in t):
                return None
            b = np.array([float(int(x)) for x in t])
            cons = [LinearConstraint(self.A, b, b)]
            if min_total:
                cons.append(LinearConstraint(np.concatenate([np.ones(self.n), np.zeros(self.nd)]).reshape(1, -1),
                                             min_total, np.inf))
            r = milp(self.cobj, constraints=cons, integrality=np.ones(self.n + self.nd), bounds=self.bounds,
                     options={"time_limit": 30})
            if r.status == 2:
                return None
            return 'timeout' if r.status != 0 else round(float(r.fun), 6)
    return Weighted


def exps(n, geo, n1, n2, g, coef):
    out = []
    for eps in (1e-3, 1e-6):
        M = np.zeros((n1, n2), dtype=complex)
        for (i, j), o in n.items():
            M[i, j] = coef[(i, j)] * eps ** o * (eps ** (g * geo.get((i, j), 0)) if o == 0 else 1.0)
        out.append(np.sort(np.linalg.svd(M, compute_uv=False))[::-1])
    return [None if x < 1e-250 or y < 1e-250 else round(float(np.log(y / x) / np.log(1e-3)), 2) for x, y in zip(*out)]


def run(m, ours, theirs, hidden_broken=False):
    import w33_pass11097_fractional_fermion_masses as P97
    import w33_pass11098_yukawa_textures as Y
    import w33_pass11102_higher_order_quark_hierarchy as H
    import w33_so16_selection_rules as S
    Weighted = weighted_orders_class()
    US, TH = S.parse_ours(ours), S.parse_theirs(theirs)
    us, th = US[m], TH[m]
    rng = np.random.default_rng(11103 + m)
    algs, hid, nu1 = P97.prepare(us, th)
    c, oc, ow, qcol = S.sm_data(th, us)
    qb = P97.conj(qcol, 'A2')
    ferm = [f for f in us['fields'] if f['m'] == 2]
    scal = [f for f in us['fields'] if f['m'] == 6]
    pick = lambda L, col, w, Yv: [f for f in L if f['col'] == col and f['w'] == w and f['Y'] == Yv]
    Q, U, D = pick(ferm, qcol, 2, F(1, 6)), pick(ferm, qb, 1, F(-2, 3)), pick(ferm, qb, 1, F(1, 3))
    Hu, Hd = pick(scal, '1', 2, F(1, 2)), pick(scal, '1', 2, F(-1, 2))
    Sv = [s for s in scal if S.dimof(s['col']) == 1 and s['w'] == 1 and s['Y'] == 0 and (hidden_broken or s['hidden_singlet'])]
    tstar = [t for t in range(3) if len({S.classes(q['n'])[t] for q in Q}) == 3]
    rec = dict(label=us['label'], tstar=tstar, n_scalars=len(Sv), nQ=len(Q), nU=len(U), nD=len(D))
    if len(tstar) != 1:
        return m, rec
    t = tstar[0]
    vecs = [H.mod1(S.charge_vector(s, True), nu1) for s in Sv]
    w0 = [F(0)] * (nu1 + 8)

    def omat(A, B, h, E):
        n, geo = {}, {}
        for i, a in enumerate(A):
            for j, b in enumerate(B):
                if not hidden_broken and not Y.hidden_ok([a, b, h], hid, algs):
                    continue
                vs = [S.charge_vector(x, True) for x in (a, b, h)]
                if P97.cubic_ok(vs, nu1):
                    if all(x['k'] == 0 and x['l'] == 0 for x in (a, b, h)):
                        if Y.eps(Y.plane(a), Y.plane(b), Y.plane(h)) == 0:
                            continue
                        geo[(i, j)] = 0
                    else:
                        ca, cb, ch = S.classes(a['n']), S.classes(b['n']), S.classes(h['n'])
                        geo[(i, j)] = sum(1 for u in range(3) if not (ca[u] == cb[u] == ch[u]))
                    n[(i, j)] = 0
                    continue
                o = E.coupling([H.mod1(v, nu1) for v in vs], w0, min_total=1)
                if isinstance(o, (int, float)):
                    n[(i, j)] = o
        return n, geo
    res = {}
    E0 = S.FastOrders(vecs, nu1, [1] * 8)
    for name, A, B, HH in (('up', Q, U, Hu), ('down', Q, D, Hd)):
        per = []
        for k, h in enumerate(HH):
            n, geo = omat(A, B, h, E0)
            if not any(o == 0 for o in n.values()) or not h['l']:
                continue
            p0 = geo_pt(h, t, S.classes)
            coef = {kk: rng.standard_normal() + 1j * rng.standard_normal() for kk in n}
            row = dict(k=k, p0=p0, generation_points=[geo_pt(q, t, S.classes) for q in Q],
                       geo={f"{i},{j}": v for (i, j), v in sorted(geo.items())})
            for wname, (a, b, cc) in GRID:
                rel = {p0: a, (p0 + 1) % 3: b, (p0 + 2) % 3: cc}
                E = Weighted(vecs, [rel[geo_pt(s, t, S.classes)] if s['l'] else 1.0 for s in Sv], nu1, [1] * 8)
                nw, _ = omat(A, B, h, E)
                for kk in nw:
                    coef.setdefault(kk, rng.standard_normal() + 1j * rng.standard_normal())
                row[wname] = dict(orders={f"{i},{j}": o for (i, j), o in sorted(nw.items())},
                                  exps={str(g): exps(nw, geo, len(A), len(B), g, coef) for g in GS})
            per.append(row)
        res[name] = per
    rec['res'] = res
    return m, rec


def summarize():
    out = {}
    for var, path in FROZEN.items():
        d = json.loads(path.read_text())
        tab = {}
        for g in map(str, GS):
            for w, _ in GRID:
                c = Counter()
                for r in d.values():
                    for s in ('up', 'down'):
                        for row in r['res'][s]:
                            e = row[w]['exps'][g]
                            c[f"{s}|gap={None if None in e else round(e[2] - e[1], 1)}"] += 1
                tab[f"g={g}|{w}"] = dict(sorted(c.items()))
        nd = {m: max(int(k.split(',')[1]) for row in r['res']['down'] for k in row['sym']['orders']) + 1
              for m, r in d.items()}
        split_models = sorted((m for m, r in d.items()
                               if any(row['only_p1']['exps']['4.0'][2] != row['only_p1']['exps']['4.0'][1]
                                      for row in r['res']['down'])), key=int)
        out[var] = dict(models=len(d), n_scalars={m: r['n_scalars'] for m, r in d.items()}, table=tab,
                        down_split_models=split_models, models_with_6_dbar=sorted((m for m, v in nd.items() if v == 6), key=int),
                        up_split_cases=sum(v for k, x in tab.items() for kk, v in x.items()
                                           if kk.startswith('up|') and kk != 'up|gap=0.0'))
    res = dict(pass_id=11103, alignments={w: list(x) for w, x in GRID}, geo_regimes=list(GS), results=out)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True, default=str))
    for var, r in out.items():
        print(var, 'up split cases', r['up_split_cases'], 'down split models', r['down_split_models'],
              '6-dbar models', r['models_with_6_dbar'])
    return res


def main(ours, theirs):
    from multiprocessing import Pool
    import w33_pass11102_higher_order_quark_hierarchy as H
    for var, broken in (("A", False), ("B", True)):
        with Pool(6) as p:
            out = dict(p.starmap(run, [(m, ours, theirs, broken) for m in H.BEST]))
        FROZEN[var].write_text(json.dumps({str(k): v for k, v in out.items()}, indent=1, default=str))


if __name__ == "__main__":
    if len(sys.argv) > 2:
        main(*sys.argv[1:3])
    summarize()
