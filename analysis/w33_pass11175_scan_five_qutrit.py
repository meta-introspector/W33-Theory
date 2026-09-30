"""Pass 11175 scan: perfect five-qutrit Clifford gates from two independent sources.

(a) F9-linear: every unitary 5x5 matrix over F9 with all entries and all 2x2 minors nonzero (exhaustive clique search,
    rows = unit vectors with all entries nonzero, edges = orthogonal with all ten 2x2 minors nonzero).
(b) Graph states: random-restart hill climbing for weighted qutrit graphs on 10 vertices with all 126 balanced blocks
    Gamma[S, S^c] invertible (AME(10,3)); each such graph gives, for every choice of 5 input legs, a perfect gate
    S = M T (the Choi Lagrangian is {(u, M u)} with M anti-symplectic; T = (x, z) -> (x, -z) on the inputs).
Records the determinant patterns (-1 positions, canonical under row and column permutations).
Output: data/w33_pass11175_five_qutrit_scan.json.  Usage: py -3 ... [graphs_wanted] [workers]
"""
import itertools
import json
import sys
import time
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "data" / "w33_pass11175_five_qutrit_scan.json"
F9 = np.array([(a, b) for a in range(3) for b in range(3)])


def mul(a, b):
    return np.stack([(a[..., 0] * b[..., 0] - a[..., 1] * b[..., 1]) % 3, (a[..., 0] * b[..., 1] + a[..., 1] * b[..., 0]) % 3], -1)


def canon(P):
    n = P.shape[0]
    best = None
    for r in itertools.permutations(range(n)):
        Q = P[list(r)]
        cols = sorted(tuple(Q[:, j]) for j in range(n))
        key = tuple(itertools.chain(*zip(*cols)))
        if best is None or key < best:
            best = key
    return ''.join(map(str, best))


def f9_perfect(m=5, chunk=800):
    nz = [tuple(x) for x in F9 if tuple(x) != (0, 0)]
    norm = {x: (x[0] ** 2 + x[1] ** 2) % 3 for x in nz}
    R = np.array([r for r in itertools.product(nz, repeat=m) if sum(norm[x] for x in r) % 3 == 1], dtype=np.int16)
    n = len(R)
    a, b = R[..., 0], R[..., 1]
    Adj = np.zeros((n, n), bool)
    for s in range(0, n, chunk):
        A_, B_ = a[s:s + chunk], b[s:s + chunk]
        re = (A_[:, None, :] * a[None, :, :] + B_[:, None, :] * b[None, :, :]).sum(-1) % 3
        im = (A_[:, None, :] * b[None, :, :] - B_[:, None, :] * a[None, :, :]).sum(-1) % 3
        ok = (re == 0) & (im == 0)
        for i, j in itertools.combinations(range(m), 2):
            p_re = A_[:, None, i] * a[None, :, j] - B_[:, None, i] * b[None, :, j]
            p_im = A_[:, None, i] * b[None, :, j] + B_[:, None, i] * a[None, :, j]
            q_re = A_[:, None, j] * a[None, :, i] - B_[:, None, j] * b[None, :, i]
            q_im = A_[:, None, j] * b[None, :, i] + B_[:, None, j] * a[None, :, i]
            ok &= ~(((p_re - q_re) % 3 == 0) & ((p_im - q_im) % 3 == 0))
        Adj[s:s + chunk] = ok
    np.fill_diagonal(Adj, False)
    cliques = []

    def grow(cl, cand):
        if len(cl) == m:
            cliques.append(tuple(cl))
            return
        for v in cand:
            nxt = cand[cand > v]
            grow(cl + [int(v)], nxt[Adj[v, nxt]])
    import math
    from fractions import Fraction
    # orbit-stabiliser: column monomial unitaries (permutations x unit phases) act on the rows with orbits labelled by
    # the positions of norm-1 entries up to permutation, i.e. by the NUMBER of norm-1 entries.  Count the m-cliques
    # through one representative row of each orbit and weight by the orbit size; each clique has m rows.
    nnorm1 = np.array([sum(1 for x in R[i] if norm[tuple(x)] == 1) for i in range(n)])
    total, pats = Fraction(0), Counter()
    for k in sorted(set(nnorm1.tolist())):
        orbit = np.flatnonzero(nnorm1 == k)
        r = int(orbit[0])
        cliques.clear()
        nb = np.flatnonzero(Adj[r])

        def grow_any(cl, cand):
            if len(cl) == m:
                cliques.append(tuple(cl))
                return
            for v in cand:
                nxt = cand[cand > v]
                grow_any(cl + [int(v)], nxt[Adj[v, nxt]])
        grow_any([r], nb)
        weight = Fraction(len(orbit), m)
        total += weight * len(cliques)
        for cl in cliques:
            normm = np.array([[norm[tuple(R[q][c])] for c in range(m)] for q in cl])
            pats[canon((normm == 2).astype(int))] += weight
    assert total.denominator == 1
    return dict(rows=n, row_sets=int(total), unitary_matrices=int(total) * math.factorial(m),
                degree_mean=float(Adj.sum(1).mean()),
                patterns={k: int(v) for k, v in pats.items()}, patterns_integral=all(v.denominator == 1 for v in pats.values()))


# ---------- graph states ----------
NQ = 10
PAIRS = list(itertools.combinations(range(NQ), 2))
SPLITS = [S for S in itertools.combinations(range(NQ), NQ // 2) if 0 in S]
ROWS = np.array(SPLITS)
COLS = np.array([[v for v in range(NQ) if v not in S] for S in SPLITS])


def to_mat(w):
    G = np.zeros((NQ, NQ), np.int64)
    for (i, j), x in zip(PAIRS, w):
        G[i, j] = G[j, i] = x
    return G


def bad_count(G):
    B = G[:, ROWS[:, :, None], COLS[:, None, :]]
    d = np.rint(np.linalg.det(B.astype(float))).astype(np.int64) % 3
    return (d == 0).sum(-1)


def hill(rng, max_steps=5000):
    w = rng.integers(0, 3, len(PAIRS))
    cur = bad_count(to_mat(w)[None])[0]
    for step in range(max_steps):
        if cur == 0:
            return w
        nbrs = []
        for e in range(len(PAIRS)):
            for dv in (1, 2):
                v = w.copy()
                v[e] = (v[e] + dv) % 3
                nbrs.append(v)
        nbrs = np.array(nbrs)
        costs = bad_count(np.stack([to_mat(v) for v in nbrs]))
        best = costs.min()
        if best > cur or (best == cur and rng.random() < 0.3):
            v = nbrs[rng.integers(len(nbrs))]
            w, cur = v, bad_count(to_mat(v)[None])[0]
        else:
            w, cur = nbrs[rng.choice(np.flatnonzero(costs == best))], best
    return None


def inv_mod3(A):
    n = A.shape[0]
    M = np.concatenate([A % 3, np.eye(n, dtype=np.int64)], 1)
    r = 0
    for c in range(n):
        p = next((i for i in range(r, n) if M[i, c] % 3), None)
        if p is None:
            return None
        M[[r, p]] = M[[p, r]]
        M[r] = (M[r] * pow(int(M[r, c]), -1, 3)) % 3
        for i in range(n):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % 3
        r += 1
    return M[:, n:]


def gate_from_graph(G, inputs):
    """symplectic gate S (per-qutrit (x,z) ordering) whose Choi state is the graph state, inputs -> outputs"""
    outputs = [v for v in range(NQ) if v not in inputs]
    L = np.zeros((NQ, 2 * NQ), np.int64)          # rows = generators X_v Z^{G_v}
    for v in range(NQ):
        L[v, 2 * v] = 1
        for u in range(NQ):
            L[v, 2 * u + 1] = G[v, u] % 3
    cols_in = [2 * q + t for q in inputs for t in range(2)]
    cols_out = [2 * q + t for q in outputs for t in range(2)]
    A = L[:, cols_in].T        # 10 x 10: input coordinates of the generators (as columns)
    B = L[:, cols_out].T
    Ai = inv_mod3(A)
    if Ai is None:
        return None
    Mm = (B @ Ai) % 3
    T = np.diag([1, -1] * (NQ // 2)) % 3
    return (Mm @ T) % 3


def gate_patterns(G):
    n = NQ // 2
    J = np.zeros((NQ, NQ), np.int64)
    for k in range(n):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, -1
    pats, ok_all = Counter(), True
    for inputs in itertools.combinations(range(NQ), n):
        S = gate_from_graph(G, list(inputs))
        ok_all &= S is not None and np.array_equal((S.T @ J @ S) % 3, J % 3)
        B = S.reshape(n, 2, n, 2).transpose(0, 2, 1, 3)
        det = (B[..., 0, 0] * B[..., 1, 1] - B[..., 0, 1] * B[..., 1, 0]) % 3
        ok_all &= bool((det != 0).all())
        pats[canon((det == 2).astype(int))] += 1
    return pats, ok_all


def circulant_stabiliser():
    """order of {R in SL(2,3)^5 : S R S^-1 local} for the explicit F9 circulant perfect gate (Pass 11170), and the
    fraction of Sp(10,3) filled by its local orbit (24^10 / stabiliser)"""
    c = np.array([F9[i] for i in (1, 1, 5, 6, 5)])
    K = np.array([[c[(s - r) % 5] for s in range(5)] for r in range(5)])
    S = np.zeros((NQ, NQ), np.int64)
    for i in range(5):
        for j in range(5):
            x, y = K[i, j]
            S[2 * i:2 * i + 2, 2 * j:2 * j + 2] = [[x, -y], [y, x]]
    S %= 3
    J = np.zeros((NQ, NQ), np.int64)
    for k in range(5):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, -1
    Si = (-J @ S.T @ J) % 3
    G = np.array([np.array([[p, q], [r, s]]) for p, q, r, s in itertools.product(range(3), repeat=4) if (p * s - q * r) % 3 == 1])
    mask = np.ones((NQ, NQ), bool)
    for k in range(5):
        mask[2 * k:2 * k + 2, 2 * k:2 * k + 2] = False
    idx = np.array(list(itertools.product(range(24), repeat=5)))
    count = 0
    for s0 in range(0, len(idx), 200000):
        ch = idx[s0:s0 + 200000]
        R = np.zeros((len(ch), NQ, NQ), np.int64)
        for k in range(5):
            R[:, 2 * k:2 * k + 2, 2 * k:2 * k + 2] = G[ch[:, k]]
        T = np.einsum('ij,mjk,kl->mil', S, R, Si) % 3
        count += int((T[:, mask] == 0).all(axis=1).sum())
    sp10 = 3 ** 25 * 8 * 80 * 728 * 6560 * 59048
    return dict(stabiliser=count, local_orbit=24 ** 10 // count, sp10=str(sp10), orbit_fraction=(24 ** 10 // count) / sp10)


def worker(args):
    seed, want = args
    rng = np.random.default_rng(seed)
    found = []
    t0 = time.time()
    while len(found) < want and time.time() - t0 < 1500:
        w = hill(rng)
        if w is not None:
            found.append([int(x) for x in w])
    return found


if __name__ == '__main__':
    want = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    t0 = time.time()
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    if 'f9' not in res:
        res['f9'] = f9_perfect()
        OUT.write_text(json.dumps(res, indent=1))
    print('f9', {k: v for k, v in res['f9'].items() if k != 'patterns'}, res['f9']['patterns'], round(time.time() - t0), flush=True)
    with Pool(W) as p:
        graphs = [g for part in p.map(worker, [(11175 + s, want) for s in range(W)]) for g in part]
    allp, ok = Counter(), True
    for w in graphs:
        pats, o = gate_patterns(to_mat(np.array(w)))
        allp.update(pats)
        ok &= o
    res.update(graphs=graphs, graph_gate_patterns=dict(allp), graph_gates_all_perfect_symplectic=bool(ok))
    OUT.write_text(json.dumps(res, indent=1))
    print('graphs', len(graphs), dict(allp), ok, round(time.time() - t0), flush=True)
