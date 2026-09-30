"""Pass 11190: the sign structure of AME(10,3) stabilizer states, and a proof of the five-qutrit pattern conjecture.

Pass 11186 found that the 17 892 perfect five-qutrit Clifford gates coming from 71 AME(10,3) graph states realise only
two of the five determinant patterns admitted by the determinant law (10-cycle and cross+permutation) and conjectured
that C4+C6, double cross and all-(-1) never occur.  This pass proves it, by exposing a rigid combinatorial structure
that every AME(10,3) stabilizer state carries.

1. Signs.  For a 4-set K of parties, the stabilisers trivial on K form a 2-dim space W_K, and restriction to each of
   the six remaining parties is an isomorphism (anything also trivial there has weight <= 5).  Pulling back each
   party's symplectic form gives signs c_p = +-1 with sum_p c_p = 0 (mod 3) by isotropy, so every 6-set of parties is
   uniform (6|0) or split 3|3.  The det of a gate block is d_ij = -c_i c_j computed in U = outputs + {j}.
2. Local graphs.  For a 3-set A, W_A is 4-dim; the kernels of the seven restrictions are pairwise skew lines whose
   Pluecker points l_p satisfy sum lambda_p l_p = 0 on the Klein quadric, and G_pq = lambda_p lambda_q B(l_p, l_q)
   is a symmetric +-1 matrix whose row d is the split of the 6-set A^c - {d}.  Its (-1)-graph has all degrees in
   {0,3,6}: up to complement it is empty, K4+3K1, K33+K1 or prism+K1 (exhaustive over 2^21 graphs).  Empty and
   K33+K1 have rank 6 with the Gram matrix of an ELLIPTIC O-(6,3) space, impossible inside the hyperbolic Klein
   quadric O+(6,3); K4+3K1 is O+, prism+K1 has rank 4.
3. Duality.  For disjoint 4-sets K1, K2 with K1 u K2 = complement of {i,j}, isotropy between W_K1 and W_K2 forces
   det phi^K1_{j<-i} = det phi^K2_{j<-i}.
4. CP-SAT over all 210 six-set states (11 states each) with the 120 local-graph tables and the 1575 duality
   equalities: 'some local graph is K4+3K1' is INFEASIBLE, 'prism+K1 everywhere' is feasible.  So every local graph
   is prism+K1: each 3-set has exactly one uniform extension, i.e. the 30 uniform 4-sets form a Steiner system
   S(3,4,10) (the inversive plane of order 3).
5. Consequences.  A 5-set contains at most one block of an S(3,4,10), so a gate has at most one all-(-1) column:
   double cross and all-(-1) are impossible.  In a prism two columns with two (-1)s never share both, so no 4-cycle:
   C4+C6 is impossible.  Hence: cross+permutation iff the input set contains a block (180 of 252 splits), 10-cycle
   otherwise (72) -- exactly the 12 780 / 5 112 census of Pass 11186.
"""

from __future__ import annotations

import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11175_scan_five_qutrit as F  # noqa: E402

OUT = ROOT / "data" / "w33_pass11190_ame10_sign_structure.json"
TABU = ROOT / "data" / "w33_pass11186_ame10_tabu.json"
NQ = 10
SIX = list(itertools.combinations(range(NQ), 6))
SIXI = {U: i for i, U in enumerate(SIX)}
TRI = list(itertools.combinations(range(NQ), 3))


# ------------------------------------------------------------------ signs of real codes
def stabiliser(G):
    L = np.zeros((NQ, 2 * NQ), np.int64)
    for v in range(NQ):
        L[v, 2 * v] = 1
        for u in range(NQ):
            L[v, 2 * u + 1] = G[v, u] % 3
    return L


def _cols(parties):
    return [2 * q + t for q in parties for t in range(2)]


def signs(L, K):
    """c_q (q outside the 4-set K): symplectic determinant of restricting W_K to party q (basis fixed by one party)"""
    U = [q for q in range(NQ) if q not in K]
    LR = L[:, _cols(list(K) + [U[0]])]
    Ri = F.inv_mod3(LR)
    assert Ri is not None and np.array_equal((LR @ Ri) % 3, np.eye(10, dtype=np.int64))
    W = Ri[-2:, :] @ L % 3
    assert not (W[:, _cols(K)] % 3).any()
    c = {}
    for q in U:
        x1, z1, x2, z2 = W[0, 2 * q], W[0, 2 * q + 1], W[1, 2 * q], W[1, 2 * q + 1]
        c[q] = int((x1 * z2 - z1 * x2) % 3)
    return c


def split_state(U, c):
    """0 = uniform, else 1 + index of the pair that joins min(U) in its sign class"""
    if len({c[q] for q in U}) == 1:
        return 0
    m = U[0]
    rest = tuple(q for q in U[1:] if c[q] == c[m])
    return 1 + list(itertools.combinations(U[1:], 2)).index(rest)


def same_class(U, state, i, j):
    if state == 0:
        return True
    cls = {U[0], *list(itertools.combinations(U[1:], 2))[state - 1]}
    return (i in cls) == (j in cls)


def local_graph(states, A):
    P7 = [q for q in range(NQ) if q not in A]
    row = {}
    for d in P7:
        U = tuple(q for q in P7 if q != d)
        st = states[U]
        row[d] = {p: (1 if same_class(U, st, p, U[0]) else -1) for p in U}
    mu = {P7[0]: 1}
    for d in P7[1:]:
        mu[d] = row[P7[0]][d] * row[d][P7[0]]
    G = np.zeros((7, 7), np.int64)
    for a, p in enumerate(P7):
        for b, d in enumerate(P7):
            if p != d:
                G[a, b] = mu[d] * row[d][p]
    return G


def graph_type(G):
    E = (G == -1).astype(int)
    if E.sum() // 2 > 10:
        E = (G == 1).astype(int)
    ne, tri = int(E.sum() // 2), int(np.trace(E @ E @ E) // 6)
    return {(0, 0): "empty", (6, 4): "K4+3K1", (9, 0): "K33+K1", (9, 2): "prism+K1"}[(ne, tri)]


def analyse_code(G):
    L = stabiliser(G)
    C = {K: signs(L, K) for K in itertools.combinations(range(NQ), 4)}
    ok_sum = all(all(v != 0 for v in c.values()) and sum(c.values()) % 3 == 0 for c in C.values())
    states = {tuple(q for q in range(NQ) if q not in K): split_state(tuple(q for q in range(NQ) if q not in K), c)
              for K, c in C.items()}
    types, sym = Counter(), True
    for A in TRI:
        Gm = local_graph(states, A)
        sym &= bool(np.array_equal(Gm, Gm.T))
        types[graph_type(Gm)] += 1
    blocks = [K for K in itertools.combinations(range(NQ), 4) if states[tuple(q for q in range(NQ) if q not in K)] == 0]
    steiner = all(sum(set(T) <= set(B) for B in blocks) == 1 for T in TRI)
    dual_bad = 0
    for i, j in itertools.combinations(range(NQ), 2):
        rest = [q for q in range(NQ) if q not in (i, j)]
        for K1 in itertools.combinations(rest, 4):
            if rest[0] in K1:
                K2 = tuple(q for q in rest if q not in K1)
                dual_bad += (C[K1][i] * C[K1][j] - C[K2][i] * C[K2][j]) % 3 != 0
    # gate patterns from the signs vs the direct gate computation
    agree, pats = True, Counter()
    for J in itertools.combinations(range(NQ), 5):
        O = [q for q in range(NQ) if q not in J]
        S = F.gate_from_graph(G, list(J))
        Bk = S.reshape(5, 2, 5, 2).transpose(0, 2, 1, 3)
        det = (Bk[..., 0, 0] * Bk[..., 1, 1] - Bk[..., 0, 1] * Bk[..., 1, 0]) % 3
        pred = np.zeros((5, 5), np.int64)
        for b, j in enumerate(J):
            K = tuple(q for q in J if q != j)
            for a, i in enumerate(O):
                pred[a, b] = (-C[K][i] * C[K][j]) % 3
        agree &= bool(np.array_equal(pred, det))
        has_block = any(set(B) <= set(J) for B in blocks)
        pats[("cross+perm" if has_block else "no-block", F.canon((det == 2).astype(int)))] += 1
    return dict(sign_sums_ok=ok_sum, local_graphs_symmetric=sym, local_graph_types=dict(types),
                uniform_six_sets=len(blocks), uniform_four_sets_are_S3_4_10=steiner, duality_violations=int(dual_bad),
                signs_predict_gate_dets=agree, patterns={f"{k[0]}|{k[1]}": v for k, v in pats.items()})


# ------------------------------------------------------------------ the local lemma
def degree_lemma():
    import networkx as nx
    E7 = list(itertools.combinations(range(7), 2))
    bits = np.arange(1 << 21, dtype=np.int64)
    deg = np.zeros((1 << 21, 7), np.int64)
    for k, (a, b) in enumerate(E7):
        e = (bits >> k) & 1
        deg[:, a] += e
        deg[:, b] += e
    cands = bits[np.isin(deg, [0, 3, 6]).all(1)]
    reps = []
    for m in cands:
        G = nx.Graph()
        G.add_nodes_from(range(7))
        G.add_edges_from(e for k, e in enumerate(E7) if (m >> k) & 1)
        if not any(nx.is_isomorphic(R, G) or nx.is_isomorphic(R, nx.complement(G)) for R in reps):
            reps.append(G)
    return len(cands), sorted((R.number_of_edges(), int(sum(nx.triangles(R).values()) // 3)) for R in reps)


BASE = {
    "empty": [],
    "K4+3K1": list(itertools.combinations(range(4), 2)),
    "K33+K1": [(a, b) for a in range(3) for b in range(3, 6)],
    "prism+K1": [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (0, 3), (1, 4), (2, 5)],
}


def quadric_types():
    """rank of G over F3 and, when 6, the O+- type of the 6-dim span: O+ iff (-1)^3 det is a square (det = 2)"""
    from sympy import GF
    from sympy.polys.matrices import DomainMatrix
    out = {}
    for name, edges in BASE.items():
        G = np.ones((7, 7), np.int64) - np.eye(7, dtype=np.int64)
        for a, b in edges:
            G[a, b] = G[b, a] = -1
        r = DomainMatrix.from_list((G % 3).tolist(), GF(3)).rank()
        dets = {int(DomainMatrix.from_list((G[np.ix_(idx, idx)] % 3).tolist(), GF(3)).det()) % 3
                for idx in ([i for i in range(7) if i != d] for d in range(7))}
        kind = "O+" if dets == {2} else "O-" if dets == {1} else "degenerate"
        out[name] = dict(rank=int(r), det6=sorted(dets), span_type=kind,
                         allowed=(kind == "O+" or r < 6))
    assert not out["empty"]["allowed"] and not out["K33+K1"]["allowed"]
    return out


# ------------------------------------------------------------------ the CP-SAT relaxation
def _tables(kinds):
    import networkx as nx
    lab = {}
    for kind in kinds:
        base = nx.Graph(BASE[kind])
        base.add_nodes_from(range(7))
        seen = set()
        for perm in itertools.permutations(range(7)):
            E = frozenset(frozenset((perm[a], perm[b])) for a, b in base.edges())
            if E not in seen:
                seen.add(E)
        lab[kind] = seen
    return lab


def build_model(kinds, need=None, fix_pattern=None, fix_blocks=None, force_at=None):
    from ortools.sat.python import cp_model
    lab = _tables(kinds)
    m = cp_model.CpModel()
    s = [m.NewIntVar(0, 10, f"s{i}") for i in range(len(SIX))]
    ind = {}
    for A in TRI:
        P7 = [q for q in range(NQ) if q not in A]
        vs = [s[SIXI[tuple(q for q in P7 if q != d)]] for d in P7]
        rows = set()
        for k, kind in enumerate(kinds):
            for E in lab[kind]:
                adj = {v: {w for e in E if v in e for w in e if w != v} for v in range(7)}
                tup = []
                for d in range(7):
                    U = tuple(P7[v] for v in range(7) if v != d)
                    if len(adj[d]) in (0, 6):
                        tup.append(0)
                    else:
                        c = {q: (1 if P7.index(q) in adj[d] else 2) for q in U}
                        tup.append(split_state(U, c))
                rows.add(tuple(tup) + (k,))
        b = m.NewIntVar(0, len(kinds) - 1, "")
        ind[A] = b
        m.AddAllowedAssignments(vs + [b], sorted(rows))
    same = {}
    for U in SIX:
        for i, j in itertools.combinations(U, 2):
            v = m.NewBoolVar("")
            m.AddAllowedAssignments([s[SIXI[U]], v], [(x, int(same_class(U, x, i, j))) for x in range(11)])
            same[U, i, j] = v
    for i, j in itertools.combinations(range(NQ), 2):
        rest = [q for q in range(NQ) if q not in (i, j)]
        for K1 in itertools.combinations(rest, 4):
            if rest[0] in K1:
                K2 = tuple(q for q in rest if q not in K1)
                m.Add(same[tuple(sorted(K2 + (i, j))), i, j] == same[tuple(sorted(K1 + (i, j))), i, j])
    if need is not None:
        m.Add(sum(ind.values()) >= 1)
    if force_at is not None:
        m.Add(ind[force_at] == need)
    if fix_pattern is not None:              # inputs 0..4, outputs 5..9, pattern[a][b] = 1 where d_{5+a, b} = -1
        J, O = list(range(5)), list(range(5, 10))
        for b, j in enumerate(J):
            U = tuple(sorted(O + [j]))
            for a, i in enumerate(O):
                m.Add(same[U, min(i, j), max(i, j)] == int(fix_pattern[a][b]))
    if fix_blocks is not None:
        for U in SIX:
            K = tuple(q for q in range(NQ) if q not in U)
            if K in fix_blocks:
                m.Add(s[SIXI[U]] == 0)
            else:
                m.Add(s[SIXI[U]] != 0)
    return m, s


def solve(m, workers=8, enumerate_all=False, s=None):
    from ortools.sat.python import cp_model
    sol = cp_model.CpSolver()
    sol.parameters.max_time_in_seconds = 7200
    if enumerate_all:
        sol.parameters.num_workers = 1               # enumeration needs ONE worker
        sol.parameters.enumerate_all_solutions = True

        class Cb(cp_model.CpSolverSolutionCallback):
            def __init__(self):
                super().__init__()
                self.n = 0

            def on_solution_callback(self):
                self.n += 1
        cb = Cb()
        st = sol.Solve(m, cb)
        return sol.StatusName(st), cb.n
    sol.parameters.num_workers = workers
    st = sol.Solve(m)
    return sol.StatusName(st), None


PATTERNS = {   # labelled representatives; rows = outputs, cols = inputs, 1 = det -1
    "C10": [[1, 1, 0, 0, 0], [0, 1, 1, 0, 0], [0, 0, 1, 1, 0], [0, 0, 0, 1, 1], [1, 0, 0, 0, 1]],
    "cross+perm": [[1, 1, 1, 1, 1], [1, 1, 0, 0, 0], [1, 0, 1, 0, 0], [1, 0, 0, 1, 0], [1, 0, 0, 0, 1]],
    "C4+C6": [[1, 1, 0, 0, 0], [1, 1, 0, 0, 0], [0, 0, 1, 1, 0], [0, 0, 0, 1, 1], [0, 0, 1, 0, 1]],
    "double-cross": [[1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 0, 0, 0], [1, 1, 0, 0, 0], [1, 1, 0, 0, 0]],
    "all": [[1] * 5 for _ in range(5)],
}


def run():
    t0 = time.time()
    D = json.load(open(TABU))
    graphs = [F.to_mat(np.array(w)) for w in D["graphs"]]
    codes = [analyse_code(G) for G in graphs]
    agg = Counter()
    for c in codes:
        agg.update(c["patterns"])
    res = dict(pass_id=11190, n_graphs=len(graphs),
               all_sign_sums_zero=all(c["sign_sums_ok"] for c in codes),
               all_local_graphs_symmetric=all(c["local_graphs_symmetric"] for c in codes),
               local_graph_types=dict(sum((Counter(c["local_graph_types"]) for c in codes), Counter())),
               uniform_six_sets={int(k): v for k, v in Counter(c["uniform_six_sets"] for c in codes).items()},
               all_S3_4_10=all(c["uniform_four_sets_are_S3_4_10"] for c in codes),
               duality_violations=sum(c["duality_violations"] for c in codes),
               signs_predict_all_gate_dets=all(c["signs_predict_gate_dets"] for c in codes),
               pattern_census=dict(agg))
    n, classes = degree_lemma()
    res["degree_lemma"] = dict(graphs_with_degrees_0_3_6=int(n), classes_up_to_complement=classes)
    res["quadric_types"] = quadric_types()
    kinds = ["prism+K1", "K4+3K1"]
    res["csp_some_K4"] = solve(build_model(kinds, need=1)[0])[0]
    res["csp_K4_at_012_one_worker"] = solve(build_model(kinds, need=1, force_at=(0, 1, 2))[0], workers=1)[0]
    res["csp_prism_only_control"] = solve(build_model(["prism+K1"])[0])[0]
    res["csp_patterns"] = {name: solve(build_model(["prism+K1"], fix_pattern=P)[0])[0] for name, P in PATTERNS.items()}
    blocks0 = set()
    L0 = stabiliser(graphs[0])
    for K in itertools.combinations(range(NQ), 4):
        c = signs(L0, K)
        if len(set(c.values())) == 1:
            blocks0.add(K)
    m, s = build_model(["prism+K1"], fix_blocks=blocks0)
    res["structures_with_fixed_steiner_system"] = solve(m, enumerate_all=True, s=s)
    res["seconds"] = round(time.time() - t0)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    for k, v in res.items():
        print(k, ":", v)


if __name__ == "__main__":
    main()
