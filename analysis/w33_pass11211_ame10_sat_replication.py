#!/usr/bin/env python3
"""Pass 11211 (independent replication of Pass 11190): the sign structure of AME(10,3) stabilizer states, and a proof of the five-qutrit pattern conjecture.

Pass 11186 conjectured that a perfect five-qutrit Clifford gate has orientation pattern 'C10' or 'cross+perm' only
('C4+C6', 'double cross', 'all' never occur).  Its Choi state is an AME(10,3) stabilizer state: a Lagrangian L in
W_1 + ... + W_10 (W_p the symplectic plane of party p) with L cap W_T = 0 for every 5-set T.

SIGNS.  For a 4-set K, V_K = L cap W_{K^c} is 2-dimensional and every projection pi_p : V_K -> W_p (p outside K) is an
isomorphism.  With a basis (a, b) of V_K, s_p = omega_p(a_p, b_p) = det pi_p = +-1; the partition of the 6-set U = K^c
into {s = +1}, {s = -1} is basis-free.  Isotropy of L gives sum_p s_p = 0 mod 3, so U is UNIFORM (6 + 0) or split 3 + 3.
DICTIONARY (from the Choi Lagrangian {(u, S T u)}, T = (x, z) -> (x, -z) on inputs):
        det S[o, i] = -1   <=>   i and o lie in the same sign class of the 6-set {i} + O          (checked 6300/6300)
so the pattern of the gate with inputs I is read off the six-sets {i} + O.
THE KLEIN QUADRIC.  For a 3-set A, V_A = L cap W_{A^c} is 4-dimensional; beta_r = pi_r^* omega_r (r outside A) are 7
rank-2 alternating forms on V_A -- points of the Klein quadric Q+(5,3) in Lambda^2 V_A^* -- with radicals V_{A+r}
pairwise transverse (AME), so beta_r ^ beta_r' != 0, and sum_r beta_r = 0 (isotropy).  With G_rr' = beta_r ^ beta_r'
(a volume form fixes G up to a global sign):
        sigma_{A^c - r'}(r, r'') = G_rr' G_r''r'          and every row of G sums to 0 mod 3,
so the SIGN GRAPH Gamma_A (edges where G = -1, up to complement) has all degrees in {0, 3, 6}: 1052 labelled graphs in
four classes up to complement -- 7K1, K4+3K1, prism+K1, K33+K1.  Realisability as a Gram matrix in the 6-dimensional
hyperbolic Klein space kills 7K1 and K33+K1 (both give a nondegenerate form of elliptic type); prism+K1 (rank 4) and
K4+3K1 survive locally.
DUALITY (isotropy between V_K and V_K' for disjoint 4-sets K, K' leaving {i, j}):  sigma_{K^c}(i, j) = sigma_{K'^c}(i, j).
THE SAT PROOF.  Variables: the sign partition of each of the 210 six-sets (11 values).  Constraints: every 7-set carries
a Klein-realisable sign graph, and duality (1575 pairs).  CaDiCaL (confirmed by Glucose4):
  * Q1  a K4+3K1 sign graph on {0..6}                                   UNSAT
    => every Gamma_A is prism+K1, whose unique isolated vertex is the only uniform 4-set through A:
       THE UNIFORM 4-SETS OF EVERY AME(10,3) STABILIZER STATE FORM A STEINER SYSTEM S(3,4,10).
  * Q2  a pattern other than C10 / cross+perm at inputs {0..4}           UNSAT  (controls: C10 SAT, cross+perm SAT)
    => THEOREM: every perfect five-qutrit Clifford gate has pattern C10 or cross+perm (the Pass 11186 conjecture).
Since cross+perm <=> some I - {i} is uniform <=> I contains a block, and a 5-set holds at most one block, every AME(10,3)
state gives cross+perm for exactly the 180 input sets containing a block and C10 for the other 72: the 2 : 5 ratio of
the F9 census (6 291 456 : 15 728 640) and of the 71 tabu graphs (5112 : 12780) is forced.
Positive control: the 71 AME(10,3) graphs of Pass 11186 satisfy every constraint (zero-sum, prism+K1 on all 120
seven-sets, duality on all 1575 pairs, S(3,4,10)).
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11175_scan_five_qutrit as SC  # noqa: E402

TABU = ROOT / "data" / "w33_pass11186_ame10_tabu.json"
OUT = ROOT / "data" / "w33_pass11211_ame10_sat_replication.json"
NQ = 10
SIX = list(itertools.combinations(range(NQ), 6))
SEVEN = list(itertools.combinations(range(NQ), 7))


# ---------------- linear algebra mod 3 ----------------
def nullspace_mod3(A):
    A = A.copy() % 3
    m, n = A.shape
    piv, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, m) if A[i, c]), None)
        if p is None:
            continue
        A[[r, p]] = A[[p, r]]
        A[r] = (A[r] * pow(int(A[r, c]), -1, 3)) % 3
        for i in range(m):
            if i != r and A[i, c]:
                A[i] = (A[i] - A[i, c] * A[r]) % 3
        piv.append(c)
        r += 1
    basis = []
    for f in (c for c in range(n) if c not in piv):
        v = np.zeros(n, np.int64)
        v[f] = 1
        for i, c in enumerate(piv):
            v[c] = (-A[i, f]) % 3
        basis.append(v)
    return np.array(basis)


def form_rank_disc(A):
    """rank and discriminant (of the nondegenerate part) of a symmetric form mod 3"""
    A = np.array(A, np.int64) % 3
    d, r, disc = len(A), 0, 1
    for k in range(d):
        idx = list(range(k, d))
        if not any(A[i, i] for i in idx):
            pq = next(((p, q) for p in idx for q in idx if p < q and A[p, q]), None)
            if pq is None:
                break
            p, q = pq
            A[p] = (A[p] + A[q]) % 3
            A[:, p] = (A[:, p] + A[:, q]) % 3
        j = next(i for i in idx if A[i, i])
        A[[k, j]] = A[[j, k]]
        A[:, [k, j]] = A[:, [j, k]]
        piv = int(A[k, k])
        for i in range(k + 1, d):
            f = (A[i, k] * pow(piv, -1, 3)) % 3
            A[i] = (A[i] - f * A[k]) % 3
            A[:, i] = (A[:, i] - f * A[:, k]) % 3
        r += 1
        disc = (disc * piv) % 3
    return r, disc


# ---------------- signs of a concrete AME(10,3) graph state ----------------
def stab_matrix(G):
    L = np.zeros((NQ, 2 * NQ), np.int64)
    for v in range(NQ):
        L[v, 2 * v] = 1
        for u in range(NQ):
            L[v, 2 * u + 1] = G[v, u] % 3
    return L


def signs(G):
    """U -> tuple of +-1 (relative to U[0]) for the 210 six-sets"""
    L = stab_matrix(G)
    out = {}
    for U in SIX:
        K = [p for p in range(NQ) if p not in U]
        N = nullspace_mod3(L[:, [2 * p + t for p in K for t in range(2)]].T)
        assert len(N) == 2
        a, b = (N[0] @ L) % 3, (N[1] @ L) % 3
        s = [int((a[2 * p] * b[2 * p + 1] - a[2 * p + 1] * b[2 * p]) % 3) for p in U]
        assert all(s)
        out[U] = tuple(1 if x == s[0] else -1 for x in s)
    return out


def to_value(U, s):
    """the sign structure as one of 11 values: None (uniform) or the 3-class containing U[0]"""
    if len(set(s)) == 1:
        return None
    return frozenset(u for u, x in zip(U, s) if x == s[0])


# ---------------- local classification of sign graphs ----------------
def sign_graphs():
    """labelled graphs on 7 vertices with all degrees in {0,3,6}, with their class and Klein realisability"""
    out = []
    for bits in itertools.product([0, 1], repeat=21):
        E = [e for e, b in zip(itertools.combinations(range(7), 2), bits) if b]
        deg = [0] * 7
        for a, b in E:
            deg[a] += 1
            deg[b] += 1
        if not all(x in (0, 3, 6) for x in deg):
            continue
        G = np.ones((7, 7), np.int64) - np.eye(7, dtype=np.int64)
        for a, b in E:
            G[a, b] = G[b, a] = 2
        r, disc = form_rank_disc(G[:6, :6])          # the Gram form on F_3^7 / <1>
        realisable = (6 - r) >= 1 or disc == 2        # t = 0 needs hyperbolic type: -disc a square
        g = nx.Graph()
        g.add_nodes_from(range(7))
        g.add_edges_from(E)
        h = g if len(E) <= 10 else nx.complement(g)
        comps = sorted(len(c) for c in nx.connected_components(h))
        name = {0: '7K1', 6: 'K4+3K1'}.get(h.number_of_edges())
        if name is None:
            six = h.subgraph(max(nx.connected_components(h), key=len))
            name = 'prism+K1' if nx.is_isomorphic(six, nx.circular_ladder_graph(3)) else 'K33+K1'
        out.append(dict(E=E, graph=g, name=name, rank=r, disc=disc, realisable=realisable, comps=comps))
    return out


def induced_tuple(g, Y):
    out = []
    for rpos, r in enumerate(Y):
        U = tuple(y for y in Y if y != r)
        nb = {Y[q] for q in g.neighbors(rpos)}
        cls = {u for u in U if u in nb}
        if len(cls) in (0, 6):
            out.append((U, None))
        else:
            out.append((U, frozenset(cls if U[0] in cls else set(U) - cls)))
    return tuple(out)


def part_values(U):
    return [None] + [frozenset((U[0],) + c) for c in itertools.combinations(U[1:], 2)]


def same(val, i, j):
    return True if val is None else (i in val) == (j in val)


def duality_pairs():
    for i, j in itertools.combinations(range(NQ), 2):
        rest = [p for p in range(NQ) if p not in (i, j)]
        for K in itertools.combinations(rest, 4):
            if rest[0] in K:
                Kp = tuple(p for p in rest if p not in K)
                yield tuple(sorted(Kp + (i, j))), tuple(sorted(K + (i, j))), i, j


# ---------------- gate patterns ----------------
def shape(P):
    r, c = P.sum(1), P.sum(0)
    if (r == 5).all() and (c == 5).all():
        return 'all'
    fr, fc = int((r == 5).sum()), int((c == 5).sum())
    if fr == 0 and fc == 0 and (r == 2).all() and (c == 2).all():
        g = nx.Graph([(('o', a), ('i', b)) for a in range(5) for b in range(5) if P[a, b]])
        cyc = sorted(len(x) for x in nx.connected_components(g))
        return 'C10' if cyc == [10] else ('C4+C6' if cyc == [4, 6] else 'other')
    if fr == 1 and fc == 1:
        Q = np.delete(np.delete(P, int(np.argmax(r == 5)), 0), int(np.argmax(c == 5)), 1)
        if (Q.sum(0) == 1).all() and (Q.sum(1) == 1).all():
            return 'cross+perm'
    if fr == 2 and fc == 2 and (r[r < 5] == 2).all() and (c[c < 5] == 2).all():
        return 'double cross'
    return 'other'


def pattern_from_values(I, vals):
    O = [p for p in range(NQ) if p not in I]
    P = np.zeros((5, 5), np.int64)
    for ci, i in enumerate(I):
        for ro, o in enumerate(O):
            P[ro, ci] = same(vals[ci], i, o)
    return P


# ---------------- SAT model ----------------
def build_cnf(graphs, kinds):
    var = {}

    def v(key):
        if key not in var:
            var[key] = len(var) + 1
        return var[key]
    cls = []
    for U in SIX:
        lits = [v(('x', U, val)) for val in part_values(U)]
        cls.append(lits)
        cls += [[-a, -b] for a, b in itertools.combinations(lits, 2)]
    ykind = {}
    for Y in SEVEN:
        ys, seen = [], set()
        for gd in graphs:
            if gd['name'] not in kinds:
                continue
            t = induced_tuple(gd['graph'], Y)
            if t in seen:
                continue
            seen.add(t)
            y = v(('y', Y, len(ys)))
            ys.append(y)
            ykind[y] = gd['name']
            cls += [[-y, v(('x', U, val))] for U, val in t]
        cls.append(ys)
    for U, Up, i, j in duality_pairs():
        for a in part_values(U):
            for b in part_values(Up):
                if same(a, i, j) != same(b, i, j):
                    cls.append([-v(('x', U, a)), -v(('x', Up, b))])
    return var, cls, ykind


def solve(cls, extra=(), solver='cadical'):
    from pysat.solvers import Cadical153, Glucose4
    S = {'cadical': Cadical153, 'glucose': Glucose4}[solver]
    with S(bootstrap_with=list(cls) + list(extra)) as s:
        return bool(s.solve())


def exists_clauses(var, Us, combos):
    """clauses asserting that the six-sets Us take one of the value tuples in combos"""
    extra, zs, nxt = [], [], len(var) + 1
    for vals in combos:
        zs.append(nxt)
        extra += [[-nxt, var[('x', U, val)]] for U, val in zip(Us, vals)]
        nxt += 1
    extra.append(zs)
    return extra


def positive_control(graphs_local, max_graphs=None, dictionary_graphs=2):
    data = json.loads(TABU.read_text())['graphs']
    realisable = {induced_tuple(gd['graph'], tuple(range(7))): gd['name'] for gd in graphs_local}
    res = dict(graphs=0, zero_sum=True, all_prism=True, duality=True, steiner=True, dictionary_checked=0,
               dictionary_ok=True, pattern_counts=Counter())
    for gi, w in enumerate(data[:max_graphs]):
        G = SC.to_mat(np.array(w))
        sg = signs(G)
        vals = {U: to_value(U, s) for U, s in sg.items()}
        res['graphs'] += 1
        res['zero_sum'] &= all(sum(s) % 3 == 0 for s in sg.values())
        for Y in SEVEN:
            relabel = {y: k for k, y in enumerate(Y)}
            key = tuple((tuple(relabel[u] for u in U),
                         None if vals[U] is None else frozenset(relabel[u] for u in vals[U]))
                        for U in (tuple(y for y in Y if y != r) for r in Y))
            res['all_prism'] &= realisable.get(key) == 'prism+K1'
        res['duality'] &= all(same(vals[U], i, j) == same(vals[Up], i, j) for U, Up, i, j in duality_pairs())
        blocks = [tuple(p for p in range(NQ) if p not in U) for U in SIX if vals[U] is None]
        cover = Counter(t for B in blocks for t in itertools.combinations(B, 3))
        res['steiner'] &= len(blocks) == 30 and len(cover) == 120 and set(cover.values()) == {1}
        for I in itertools.combinations(range(NQ), 5):
            O = [p for p in range(NQ) if p not in I]
            P = pattern_from_values(I, [vals[tuple(sorted([i] + O))] for i in I])
            res['pattern_counts'][shape(P)] += 1
            if gi < dictionary_graphs:
                S = SC.gate_from_graph(G, list(I))
                B = S.reshape(5, 2, 5, 2).transpose(0, 2, 1, 3)
                det = (B[..., 0, 0] * B[..., 1, 1] - B[..., 0, 1] * B[..., 1, 0]) % 3
                res['dictionary_ok'] &= bool(((det == 2) == (P == 1)).all())
                res['dictionary_checked'] += 25
    res['pattern_counts'] = dict(res['pattern_counts'])
    return res


def summarize(max_graphs=None):
    t0 = time.time()
    graphs = sign_graphs()
    classes = Counter((g['name'], g['realisable'], g['rank'], g['disc']) for g in graphs)
    local = dict(labelled_graphs=len(graphs),
                 classes={f"{n}|realisable={r}|rank={k}|disc={d}": c for (n, r, k, d), c in sorted(classes.items())},
                 realisable_classes=sorted({g['name'] for g in graphs if g['realisable']}))
    control = positive_control(graphs, max_graphs)
    var, cls, ykind = build_cnf(graphs, ('prism+K1', 'K4+3K1'))
    Y0 = tuple(range(7))
    k4 = [i for k, i in var.items() if k[0] == 'y' and k[1] == Y0 and ykind[i] == 'K4+3K1']
    I0 = (0, 1, 2, 3, 4)
    O0 = tuple(p for p in range(NQ) if p not in I0)
    Us = [tuple(sorted((i,) + O0)) for i in I0]
    combos = {}
    for vals in itertools.product(*[part_values(U) for U in Us]):
        combos.setdefault(shape(pattern_from_values(I0, vals)), []).append(vals)
    bad = [v for k, vs in combos.items() if k not in ('C10', 'cross+perm') for v in vs]
    sat = dict(variables=len(var), clauses=len(cls),
               base=solve(cls),
               Q1_K4_3K1_somewhere=solve(cls, [k4]),
               Q2_forbidden_pattern=solve(cls, exists_clauses(var, Us, bad)),
               control_C10=solve(cls, exists_clauses(var, Us, combos['C10'])),
               control_cross_perm=solve(cls, exists_clauses(var, Us, combos['cross+perm'])),
               glucose_Q1=solve(cls, [k4], 'glucose'),
               glucose_Q2=solve(cls, exists_clauses(var, Us, bad), 'glucose'),
               value_combinations={k: len(v) for k, v in combos.items()})
    res = dict(pass_id=11211, local=local, positive_control=control, sat=sat,
               theorem_steiner=not sat['Q1_K4_3K1_somewhere'] and not sat['glucose_Q1'],
               theorem_patterns=(not sat['Q2_forbidden_pattern'] and not sat['glucose_Q2'] and sat['control_C10']
                                 and sat['control_cross_perm']),
               seconds=round(time.time() - t0, 1))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
