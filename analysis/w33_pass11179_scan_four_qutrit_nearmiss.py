"""Pass 11179 scan: how close can four qutrits come to a perfect tick?  AME(8,3) does not exist (shadow inequalities), so
every 8-qutrit stabilizer state (equivalently, up to local Cliffords, every weighted qutrit graph on 8 vertices) has at
least one balanced cut S|S^c (|S| = 4) with Gamma[S, S^c] singular.  Tabu search minimises the number of such cuts
(35 in all) and records the best graphs with their deficient cuts and ranks.
Output: data/w33_pass11179_nearmiss_scan.json.  Usage: py -3 ... [restarts_per_worker] [steps] [workers]
"""
import itertools
import json
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "data" / "w33_pass11179_nearmiss_scan.json"
NQ = 8
PAIRS = list(itertools.combinations(range(NQ), 2))
SPLITS = [S for S in itertools.combinations(range(NQ), 4) if 0 in S]
ROWS = np.array(SPLITS)
COLS = np.array([[v for v in range(NQ) if v not in S] for S in SPLITS])
# neighbour index tables: move e -> value
MOVES = [(e, dv) for e in range(len(PAIRS)) for dv in (1, 2)]


def to_mats(W):
    G = np.zeros((W.shape[0], NQ, NQ), np.int64)
    for k, (i, j) in enumerate(PAIRS):
        G[:, i, j] = G[:, j, i] = W[:, k]
    return G


def bad(W):
    G = to_mats(W)
    B = G[:, ROWS[:, :, None], COLS[:, None, :]]
    d = np.rint(np.linalg.det(B.astype(float))).astype(np.int64) % 3
    return (d == 0).sum(-1)


R3 = [(list(r), list(c)) for r in itertools.combinations(range(4), 3) for c in itertools.combinations(range(4), 3)]
R2 = [(list(r), list(c)) for r in itertools.combinations(range(4), 2) for c in itertools.combinations(range(4), 2)]


def deficit(W):
    """total entanglement deficit sum over the 35 balanced cuts of (4 - rank_F3 Gamma[S, S^c])"""
    G = to_mats(W)
    B = G[:, ROWS[:, :, None], COLS[:, None, :]]                       # m, 35, 4, 4
    d4 = np.rint(np.linalg.det(B.astype(float))).astype(np.int64) % 3 != 0
    d3 = np.zeros(d4.shape, bool)
    for r, c in R3:
        d3 |= np.rint(np.linalg.det(B[:, :, r][:, :, :, c].astype(float))).astype(np.int64) % 3 != 0
    d2 = np.zeros(d4.shape, bool)
    for r, c in R2:
        d2 |= np.rint(np.linalg.det(B[:, :, r][:, :, :, c].astype(float))).astype(np.int64) % 3 != 0
    d1 = (B % 3 != 0).any(axis=(-1, -2))
    rank = np.where(d4, 4, np.where(d3, 3, np.where(d2, 2, np.where(d1, 1, 0))))
    return (4 - rank).sum(-1)


OBJECTIVE = {'count': None, 'deficit': None}


def tabu(rng, steps, tenure=7, objective='count'):
    f = bad if objective == 'count' else deficit
    w = rng.integers(0, 3, len(PAIRS))
    cur = f(w[None])[0]
    best, best_w = cur, w.copy()
    tabu_until = np.zeros(len(PAIRS), int)
    for it in range(steps):
        N = np.repeat(w[None], len(MOVES), 0)
        for k, (e, dv) in enumerate(MOVES):
            N[k, e] = (N[k, e] + dv) % 3
        c = f(N)
        allowed = np.array([tabu_until[e] <= it or c[k] < best for k, (e, dv) in enumerate(MOVES)])
        c_eff = np.where(allowed, c, 999)
        m = c_eff.min()
        k = rng.choice(np.flatnonzero(c_eff == m))
        e = MOVES[k][0]
        w, cur = N[k], c[k]
        tabu_until[e] = it + tenure + rng.integers(0, 4)
        if cur < best:
            best, best_w = cur, w.copy()
            if best == 0:
                break
    return int(best), best_w


def rank3(A):
    A = A.copy() % 3
    r = 0
    rows, cols = A.shape
    for c in range(cols):
        p = next((i for i in range(r, rows) if A[i, c]), None)
        if p is None:
            continue
        A[[r, p]] = A[[p, r]]
        A[r] = (A[r] * pow(int(A[r, c]), -1, 3)) % 3
        for i in range(rows):
            if i != r and A[i, c]:
                A[i] = (A[i] - A[i, c] * A[r]) % 3
        r += 1
    return r


def profile(w):
    G = to_mats(np.array(w)[None])[0]
    cuts = []
    for S in SPLITS:
        Sc = [v for v in range(NQ) if v not in S]
        r = rank3(G[np.ix_(list(S), Sc)])
        if r < 4:
            cuts.append((list(S), r))
    three_uniform = all(rank3(G[np.ix_(list(T), [v for v in range(NQ) if v not in T])]) == 3
                        for T in itertools.combinations(range(NQ), 3))
    return cuts, three_uniform


def worker(args):
    seed, restarts, steps, objective = args
    rng = np.random.default_rng(seed)
    dist, bests = Counter(), []
    for _ in range(restarts):
        b, w = tabu(rng, steps, objective=objective)
        dist[b] += 1
        if b <= (2 if objective == 'count' else 4) and len(bests) < 30:
            bests.append((b, [int(x) for x in w]))
    return dist, bests


if __name__ == '__main__':
    R = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
    W = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    objective = sys.argv[4] if len(sys.argv) > 4 else 'count'
    if objective == 'deficit':
        OUT = OUT.with_name('w33_pass11179_nearmiss_deficit_scan.json')
    with Pool(W) as p:
        res = p.map(worker, [(11179 + 100 * (objective == 'deficit') + s, R, steps, objective) for s in range(W)])
    dist, bests = Counter(), []
    for d, b in res:
        dist.update(d)
        bests += b
    mn = min(dist)
    prof = []
    for b, w in sorted(bests)[:20]:
        cuts, u3 = profile(w)
        prof.append(dict(bad=b, weights=w, deficient_cuts=cuts, three_uniform=u3))
    out = dict(objective=objective, restarts=sum(dist.values()), steps=steps, distribution={str(k): v for k, v in sorted(dist.items())},
               minimum=int(mn), best=prof)
    OUT.write_text(json.dumps(out, indent=1))
    print(dict(sorted(dist.items())), 'min', mn, flush=True)
    for p_ in prof[:5]:
        print(p_['bad'], p_['deficient_cuts'], p_['three_uniform'], flush=True)
