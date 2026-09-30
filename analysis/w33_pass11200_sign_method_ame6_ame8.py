"""Pass 11200: the sign method of Pass 11190 on AME(6,3) and AME(8,3).

For an AME(2m,3) stabilizer state and an (m-1)-set K, the stabilisers trivial on K form a 2-dim space restricting
isomorphically to each of the m+1 other parties; the pulled-back symplectic signs sum to 0 mod 3 (isotropy).
  m = 2 (AME(4,3)):  3-sets are uniform            -> two-qutrit perfect gates have all block dets -1;
  m = 3 (AME(6,3)):  4-sets split 2|2              -> three-qutrit perfect gates have permutation patterns;
  m = 4 (AME(8,3)):  5-sets split 4|1 (never 5|0)  ;
  m = 5 (AME(10,3)): 6-sets split 3|3 or 6|0       (Pass 11190).
For an (m-2)-set A the signs of the (m+1)-sets inside A^c assemble into a symmetric +-1 matrix G on the m+2 parties
outside A (zero diagonal, rows summing to 0 mod 3; Pluecker points of pairwise skew lines on the Klein quadric).

  AME(6,3): the local graphs on 5 parties are 2-regular: pentagons.  Sign structures of random AME(6,3) graph states
            are computed and compared (automorphism group, isomorphism).
  AME(8,3): local graphs on 6 parties have all degrees in {1,4}.  CP-SAT over the 56 five-set states with the
            28 local-graph tables and the isotropy duality decides whether ANY sign structure exists.  AME(8,3) is known
            not to exist (shadow inequalities, Huber-Eltschka-Siewert-Guehne, arXiv:1708.06298); an infeasible model
            is an independent, elementary proof for stabilizer states.
"""

from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11175_scan_five_qutrit as F  # noqa: E402

OUT = ROOT / "data" / "w33_pass11200_sign_method_ame6_ame8.json"


def stabiliser(G):
    n = len(G)
    L = np.zeros((n, 2 * n), np.int64)
    for v in range(n):
        L[v, 2 * v] = 1
        for u in range(n):
            L[v, 2 * u + 1] = G[v, u] % 3
    return L


def cols(parties):
    return [2 * q + t for q in parties for t in range(2)]


def is_ame(G):
    n = len(G)
    for S in itertools.combinations(range(n), n // 2):
        if 0 not in S:
            continue
        T = [v for v in range(n) if v not in S]
        if F.inv_mod3(G[np.ix_(list(S), T)] % 3) is None:
            return False
    return True


def signs(L, K):
    n = L.shape[0]
    U = [q for q in range(n) if q not in K]
    LR = L[:, cols(list(K) + U[: n // 2 - len(K)])]
    Ri = F.inv_mod3(LR)
    W = Ri[-2:, :] @ L % 3
    assert not (W[:, cols(K)] % 3).any()
    return {q: int((W[0, 2 * q] * W[1, 2 * q + 1] - W[0, 2 * q + 1] * W[1, 2 * q]) % 3) for q in U}


def random_ame(n, rng, tries=2_000_000):
    iu = np.triu_indices(n, 1)
    for _ in range(tries):
        G = np.zeros((n, n), np.int64)
        G[iu] = rng.integers(0, 3, len(iu[0]))
        G = (G + G.T) % 3
        if is_ame(G):
            return G
    return None


def structure(G, m):
    n = 2 * m
    L = stabiliser(G)
    out = {}
    for K in itertools.combinations(range(n), m - 1):
        c = signs(L, K)
        assert all(v for v in c.values()) and sum(c.values()) % 3 == 0
        U = tuple(sorted(c))
        out[U] = frozenset(q for q in U if c[q] == c[U[0]])      # class of min(U)
    return out


def local_degrees(struct, m):
    """(-1)-graph degrees of G for each (m-2)-set A, from the splits (G symmetric up to row rescaling)"""
    n = 2 * m
    res = Counter()
    for A in itertools.combinations(range(n), m - 2):
        P_ = [q for q in range(n) if q not in A]
        row = {}
        for d in P_:
            U = tuple(q for q in P_ if q != d)
            cls = struct[U]
            row[d] = {p: (1 if (p in cls) == (U[0] in cls) else -1) for p in U}
        mu = {P_[0]: 1}
        for d in P_[1:]:
            mu[d] = row[P_[0]][d] * row[d][P_[0]]
        Gm = np.array([[0 if p == d else mu[d] * row[d][p] for d in P_] for p in P_])
        assert np.array_equal(Gm, Gm.T)
        E = (Gm == -1).astype(int)
        if E.sum() // 2 > len(P_) * (len(P_) - 1) // 4:
            E = (Gm == 1).astype(int)
        res[tuple(sorted(E.sum(0)))] += 1
    return dict((str(k), v) for k, v in res.items())


def aut_order(struct, n):
    """automorphism group order of a sign structure (splits of (m+1)-sets), brute force over S_n (n <= 8)"""
    cnt = 0
    items = list(struct.items())
    for p in itertools.permutations(range(n)):
        ok = True
        for U, cls in items:
            V = tuple(sorted(p[q] for q in U))
            img = frozenset(p[q] for q in cls)
            tgt = struct[V]
            if not (img == tgt or img == frozenset(V) - tgt):
                ok = False
                break
        cnt += ok
    return cnt


def csp_ame8():
    """CP-SAT: does any sign structure on 8 parties satisfy the local {1,4}-degree graphs and the isotropy duality?"""
    import networkx as nx
    from ortools.sat.python import cp_model
    n, m = 8, 4
    FIVE = list(itertools.combinations(range(n), 5))
    FI = {U: i for i, U in enumerate(FIVE)}
    # graphs on 6 labelled vertices with all degrees in {1,4}
    E6 = list(itertools.combinations(range(6), 2))
    graphs = []
    for mask in range(1 << 15):
        deg = [0] * 6
        for k, (a, b) in enumerate(E6):
            if mask >> k & 1:
                deg[a] += 1
                deg[b] += 1
        if all(d in (1, 4) for d in deg):
            graphs.append([E6[k] for k in range(15) if mask >> k & 1])
    classes = []
    for g in graphs:
        H = nx.Graph(g)
        H.add_nodes_from(range(6))
        if not any(nx.is_isomorphic(H, R) or nx.is_isomorphic(nx.complement(H), R) for R in classes):
            classes.append(H)
    model = cp_model.CpModel()
    s = [model.NewIntVar(0, 4, "") for _ in FIVE]          # state = which of the 5 is the singleton
    for A in itertools.combinations(range(n), 2):
        P_ = [q for q in range(n) if q not in A]
        rows = set()
        for g in graphs:
            adj = {v: set() for v in range(6)}
            for a, b in g:
                adj[a].add(b)
                adj[b].add(a)
            tup = []
            for d in range(6):
                U = tuple(P_[v] for v in range(6) if v != d)
                # row d of G: -1 on neighbours; the singleton class is the unique party whose sign differs
                minus = {P_[v] for v in adj[d]}
                single = minus if len(minus) == 1 else set(U) - minus
                assert len(single) == 1
                tup.append(U.index(next(iter(single))))
            rows.add(tuple(tup))
        model.AddAllowedAssignments([s[FI[tuple(q for q in P_ if q != d)]] for d in P_], sorted(rows))
    same = {}
    for U in FIVE:
        for i, j in itertools.combinations(U, 2):
            v = model.NewBoolVar("")
            model.AddAllowedAssignments([s[FI[U]], v],
                                        [(x, int(U[x] not in (i, j))) for x in range(5)])
            same[U, i, j] = v
    for i, j in itertools.combinations(range(n), 2):
        rest = [q for q in range(n) if q not in (i, j)]
        for K1 in itertools.combinations(rest, 3):
            if rest[0] in K1:
                K2 = tuple(q for q in rest if q not in K1)
                model.Add(same[tuple(sorted(K2 + (i, j))), i, j] == same[tuple(sorted(K1 + (i, j))), i, j])
    sol = cp_model.CpSolver()
    sol.parameters.num_workers = 8
    st_dual = sol.StatusName(sol.Solve(model))
    # enumerate every combinatorial solution (ONE worker) and measure its symmetry
    sols = []

    class Cb(cp_model.CpSolverSolutionCallback):
        def __init__(self):
            super().__init__()

        def on_solution_callback(self):
            sols.append([self.Value(v) for v in s])
    en = cp_model.CpSolver()
    en.parameters.num_workers = 1
    en.parameters.enumerate_all_solutions = True
    en.Solve(model, Cb())
    auts = Counter()
    for sv in sols[:40]:
        single = {U: U[sv[i]] for i, U in enumerate(FIVE)}
        a = sum(all(p[single[U]] == single[tuple(sorted(p[q] for q in U))] for U in FIVE)
                for p in itertools.permutations(range(n)))
        auts[a] += 1
    steiner = True
    for sv in sols:
        blocks = {tuple(q for q in U if q != U[sv[i]]) for i, U in enumerate(FIVE)}
        steiner &= len(blocks) == 14 and all(sum(set(T) <= set(B) for B in blocks) == 1
                                             for T in itertools.combinations(range(n), 3))
    return dict(local_graph_classes=[(R.number_of_edges(), sorted(d for _, d in R.degree())) for R in classes],
                labelled_graphs=len(graphs), with_duality=st_dual, n_solutions=len(sols),
                aut_orders_first40=dict(auts), four_sets_minus_singleton_form_S3_4_8=bool(steiner))


def run():
    rng = np.random.default_rng(11200)
    res = dict(pass_id=11200)
    # AME(6,3)
    structs, degs = [], Counter()
    for _ in range(12):
        G = random_ame(6, rng)
        st = structure(G, 3)
        structs.append(st)
        degs.update(local_degrees(st, 3))
    auts = [aut_order(st, 6) for st in structs[:3]]
    iso = []
    for st in structs[1:]:
        found = False
        for p in itertools.permutations(range(6)):
            if all((frozenset(p[q] for q in cls) in (structs[0][tuple(sorted(p[q] for q in U))],
                                                     frozenset(p[q] for q in U) - structs[0][tuple(sorted(p[q] for q in U))]))
                   for U, cls in st.items()):
                found = True
                break
        iso.append(found)
    res["ame6"] = dict(samples=len(structs), local_graph_degrees=dict(degs), aut_orders=auts,
                       all_isomorphic=all(iso))
    # AME(8,3)
    res["ame8"] = csp_ame8()
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
