"""Pass 11175 scan: uniform-ish samples of Sp(10,3) (products of 200 random symplectic transvections, vectorised); a
five-qutrit Clifford gate is perfect (Choi state AME(10,3)) iff all 25 party blocks and all 100 two-party 4x4 minors are
invertible (the 3- and 4-party minors follow: det S[I,J] = +- det S[I^c,J^c] for symplectic S, by Jacobi's identity and
S^{-1} = -J S^T J).  Records the fraction and the determinant patterns (-1 positions) of the perfect samples.
Output: data/w33_pass11175_sp10_3_sample.json.  Usage: py -3 ... [n_per_worker] [workers]
"""
import itertools
import json
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "data" / "w33_pass11175_sp10_3_sample.json"
N = 5
D = 2 * N
J = np.zeros((D, D), np.int64)
for k in range(N):
    J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, -1
PAIRS = list(itertools.combinations(range(N), 2))


def sample(rng, m, steps=200):
    M = np.broadcast_to(np.eye(D, dtype=np.int64), (m, D, D)).copy()
    for _ in range(steps):
        v = rng.integers(0, 3, (m, D))
        c = rng.integers(1, 3, (m, 1))
        w = np.einsum('mi,ij,mjk->mk', v, J, M) % 3              # (v^T J M)
        M = (M + c[:, :, None] * v[:, :, None] * w[:, None, :]) % 3
    return M


def block_dets(M):
    B = M.reshape(-1, N, 2, N, 2).transpose(0, 1, 3, 2, 4)        # m, i, j, 2, 2
    return (B[..., 0, 0] * B[..., 1, 1] - B[..., 0, 1] * B[..., 1, 0]) % 3


def perfect_mask(M):
    Db = block_dets(M)
    ok = (Db != 0).all(axis=(1, 2))
    for (i1, i2) in PAIRS:
        rows = [2 * i1, 2 * i1 + 1, 2 * i2, 2 * i2 + 1]
        for (j1, j2) in PAIRS:
            cols = [2 * j1, 2 * j1 + 1, 2 * j2, 2 * j2 + 1]
            sub = M[:, rows][:, :, cols].astype(float)
            d = np.rint(np.linalg.det(sub)).astype(np.int64) % 3
            ok &= d != 0
    return ok, Db


def canon(P):
    """canonical form of a 0/1 matrix under row and column permutations"""
    best = None
    for r in itertools.permutations(range(N)):
        Q = P[list(r)]
        cols = sorted(tuple(Q[:, j]) for j in range(N))
        key = tuple(itertools.chain(*zip(*cols)))
        if best is None or key < best:
            best = key
    return best


def worker(args):
    seed, n = args
    rng = np.random.default_rng(seed)
    tot, perf, pats, reps = 0, 0, Counter(), []
    while tot < n:
        m = min(20000, n - tot)
        M = sample(rng, m)
        ok, Db = perfect_mask(M)
        tot += m
        perf += int(ok.sum())
        for idx in np.flatnonzero(ok):
            P = (Db[idx] == 2).astype(int)
            pats[canon(P)] += 1
            if len(reps) < 20:
                reps.append(M[idx].tolist())
    return tot, perf, pats, reps


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    with Pool(W) as p:
        res = p.map(worker, [(11175 + s, n) for s in range(W)])
    tot = sum(r[0] for r in res)
    perf = sum(r[1] for r in res)
    pats = Counter()
    for r in res:
        pats.update(r[2])
    reps = [x for r in res for x in r[3]][:40]
    out = dict(samples=tot, perfect=perf, fraction=perf / tot,
               patterns={''.join(map(str, k)): v for k, v in pats.most_common()}, representatives=reps)
    OUT.write_text(json.dumps(out))
    print(tot, perf, perf / tot, len(pats), flush=True)
    for k, v in pats.most_common(10):
        print(k, v, flush=True)
