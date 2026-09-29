"""Pass 11171 scan: temporal CGLMP (one qudit, rank-1 Lueders measurements A_x then B_y) to large d.

For fixed bases the optimal state is the top eigenvector of W = sum_x A_x diag(w_x) A_x^dagger with
w_x[a] = sum_{y,b} C_xy[a,b] |<b_y|a_x>|^2, so I_d = max over the four bases of lambda_max(W).  Riemannian gradient ascent on
U(d)^4 with analytic (Wirtinger) gradients and Armijo backtracking; the envelope theorem gives the gradient at the top
eigenvector.  Validated against the committed optima for d <= 10 (Passes 11147, 11155, 11160).
Output: data/w33_pass11171_temporal_cglmp_large_d.json.  Usage: py -3 ... [dlist] [restarts]
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.linalg import expm

OUT = Path(__file__).resolve().parents[1] / "data" / "w33_pass11171_temporal_cglmp_large_d.json"


def coeff(d):
    C = np.zeros((2, 2, d, d))
    for k in range(d // 2):
        c = 1 - 2 * k / (d - 1)
        for a in range(d):
            for b in range(d):
                if a == (b + k) % d: C[0, 0, a, b] += c
                if b == (a + k + 1) % d: C[1, 0, a, b] += c
                if a == (b + k) % d: C[1, 1, a, b] += c
                if b == (a + k) % d: C[0, 1, a, b] += c
                if a == (b - k - 1) % d: C[0, 0, a, b] -= c
                if b == (a - k) % d: C[1, 0, a, b] -= c
                if a == (b - k - 1) % d: C[1, 1, a, b] -= c
                if b == (a - k - 1) % d: C[0, 1, a, b] -= c
    return C


def value_grad(U, C):
    A, B = U[:2], U[2:]
    T = [[np.abs(B[y].conj().T @ A[x]).T ** 2 for y in range(2)] for x in range(2)]        # T[x][y][a,b]
    w = [sum((C[x, y] * T[x][y]).sum(1) for y in range(2)) for x in range(2)]
    W = sum(A[x] @ np.diag(w[x]) @ A[x].conj().T for x in range(2))
    ev, V = np.linalg.eigh(W)
    f, psi = ev[-1], V[:, -1]
    p = [np.abs(A[x].conj().T @ psi) ** 2 for x in range(2)]
    G = [None] * 4
    for x in range(2):
        g = np.outer(psi, psi.conj()) @ A[x] @ np.diag(w[x])
        for y in range(2):
            Wxy = C[x, y] * p[x][:, None]                                                 # [a,b]
            g = g + B[y] @ (Wxy.T * (B[y].conj().T @ A[x]))
        G[x] = g
    for y in range(2):
        g = 0
        for x in range(2):
            Wxy = C[x, y] * p[x][:, None]
            g = g + A[x] @ (Wxy * (A[x].conj().T @ B[y]))
        G[2 + y] = g
    Om = [(U[i].conj().T @ G[i] - G[i].conj().T @ U[i]) / 2 for i in range(4)]
    return f, Om


def run(args):
    d, seed, iters = args
    rng = np.random.default_rng(seed)
    C = coeff(d)
    U = [np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))[0] for _ in range(4)]
    f, Om = value_grad(U, C)
    eta = 0.5
    for it in range(iters):
        gn = sum(np.linalg.norm(o) ** 2 for o in Om)
        if gn < 1e-20:
            break
        while True:
            U2 = [U[i] @ expm(eta * Om[i]) for i in range(4)]
            f2, Om2 = value_grad(U2, C)
            if f2 >= f + 1e-4 * eta * gn or eta < 1e-10:
                break
            eta *= 0.5
        if f2 < f:
            break
        U, f, Om = U2, f2, Om2
        eta = min(eta * 1.5, 4.0)
    return d, seed, float(f)


if __name__ == '__main__':
    ds = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [3, 4, 5, 6, 8, 10, 12, 16, 20, 24]
    R = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    iters = int(sys.argv[3]) if len(sys.argv) > 3 else 6000
    jobs = [(d, 1000 * d + s, iters) for d in ds for s in range(R)]
    with Pool(8) as pool:
        res = pool.map(run, jobs, chunksize=1)
    best = {}
    for d, s, v in res:
        best.setdefault(d, []).append(v)
    out = {str(d): dict(best=max(v), values=sorted(v)) for d, v in sorted(best.items())}
    OUT.write_text(json.dumps(out, indent=1))
    for d, v in sorted(best.items()):
        print(d, max(v), flush=True)
