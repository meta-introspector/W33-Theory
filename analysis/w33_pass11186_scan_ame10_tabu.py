"""Pass 11186 scan: AME(10,3) weighted qutrit graph states by tabu search (the plain hill climb of Pass 11175 found none),
and the orientation patterns of the perfect five-qutrit gates they give (every choice of five input legs).
Output: data/w33_pass11186_ame10_tabu.json.  Usage: py -3 ... [restarts_per_worker] [steps] [workers] [minutes]
"""
import itertools
import json
import sys
import time
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w33_pass11175_scan_five_qutrit as F  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "data" / "w33_pass11186_ame10_tabu.json"
NQ = 10
PAIRS = list(itertools.combinations(range(NQ), 2))
MOVES = [(e, dv) for e in range(len(PAIRS)) for dv in (1, 2)]


def to_mats(W):
    G = np.zeros((W.shape[0], NQ, NQ), np.int64)
    for k, (i, j) in enumerate(PAIRS):
        G[:, i, j] = G[:, j, i] = W[:, k]
    return G


def bad(W):
    G = to_mats(W)
    B = G[:, F.ROWS[:, :, None], F.COLS[:, None, :]]
    d = np.rint(np.linalg.det(B.astype(float))).astype(np.int64) % 3
    return (d == 0).sum(-1)


def tabu(rng, steps, deadline, tenure=9):
    w = rng.integers(0, 3, len(PAIRS))
    cur = bad(w[None])[0]
    best, best_w = cur, w.copy()
    until = np.zeros(len(PAIRS), int)
    for it in range(steps):
        if time.time() > deadline:
            break
        N = np.repeat(w[None], len(MOVES), 0)
        for k, (e, dv) in enumerate(MOVES):
            N[k, e] = (N[k, e] + dv) % 3
        c = bad(N)
        allowed = np.array([until[e] <= it or c[k] < best for k, (e, dv) in enumerate(MOVES)])
        ce = np.where(allowed, c, 999)
        k = rng.choice(np.flatnonzero(ce == ce.min()))
        w, cur = N[k], c[k]
        until[MOVES[k][0]] = it + tenure + rng.integers(0, 5)
        if cur < best:
            best, best_w = cur, w.copy()
            if best == 0:
                break
    return int(best), best_w


def worker(args):
    seed, restarts, steps, minutes = args
    rng = np.random.default_rng(seed)
    deadline = time.time() + 60 * minutes
    dist, found = Counter(), []
    for _ in range(restarts):
        if time.time() > deadline:
            break
        b, w = tabu(rng, steps, deadline)
        dist[b] += 1
        if b == 0:
            found.append([int(x) for x in w])
    return dist, found


if __name__ == '__main__':
    R = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
    W = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    minutes = float(sys.argv[4]) if len(sys.argv) > 4 else 40
    with Pool(W) as p:
        res = p.map(worker, [(11186 + s, R, steps, minutes) for s in range(W)])
    dist, found = Counter(), []
    for d, f in res:
        dist.update(d)
        found += f
    pats, ok = Counter(), True
    for w in found:
        pp, o = F.gate_patterns(F.to_mat(np.array(w)))
        pats.update(pp)
        ok &= o
    out = dict(restarts=sum(dist.values()), steps=steps, distribution={str(k): v for k, v in sorted(dist.items())},
               graphs=found, gate_patterns=dict(pats), gates_all_perfect_symplectic=bool(ok))
    OUT.write_text(json.dumps(out, indent=1))
    print(dict(sorted(dist.items())), 'found', len(found), dict(pats), ok, flush=True)
