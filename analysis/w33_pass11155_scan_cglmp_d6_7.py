"""Pass 11155 scan: temporal vs spatial CGLMP optima for d = 6, 7 (8 L-BFGS-B restarts each); frozen as data/w33_pass11155_cglmp_optima_d6_7.json"""
import numpy as np, json, sys
from scipy.linalg import expm
from scipy.optimize import minimize
from multiprocessing import Pool
def I_d(P, d):
    def pr(x, y, rel):
        return sum(P[x][y][a, b] for a in range(d) for b in range(d) if rel(a, b))
    tot = 0.0
    for k in range(d // 2):
        c = 1 - 2 * k / (d - 1)
        plus = (pr(0, 0, lambda a, b: a == (b + k) % d) + pr(1, 0, lambda a, b: b == (a + k + 1) % d)
                + pr(1, 1, lambda a, b: a == (b + k) % d) + pr(0, 1, lambda a, b: b == (a + k) % d))
        minus = (pr(0, 0, lambda a, b: a == (b - k - 1) % d) + pr(1, 0, lambda a, b: b == (a - k) % d)
                 + pr(1, 1, lambda a, b: a == (b - k - 1) % d) + pr(0, 1, lambda a, b: b == (a - k - 1) % d))
        tot += c * (plus - minus)
    return tot
def U(p, d):
    H = np.zeros((d, d), complex); iu = np.triu_indices(d, 1); m = len(iu[0])
    H[iu] = p[:m] + 1j * p[m:2 * m]; H = H + H.conj().T + np.diag(p[2 * m:2 * m + d]); return expm(1j * H)
def unpack(p, d, kind):
    n = d * d
    Us = [U(p[i * n:(i + 1) * n], d) for i in range(4)]
    if kind == 'temporal':
        v = p[4 * n:4 * n + d] + 1j * p[4 * n + d:4 * n + 2 * d]
    else:
        v = p[4 * n:4 * n + d * d] + 1j * p[4 * n + d * d:4 * n + 2 * d * d]
    return Us, v / np.linalg.norm(v)
def value(p, d, kind):
    (A0, A1, B0, B1), v = unpack(p, d, kind); A = [A0, A1]; B = [B0, B1]
    P = [[np.zeros((d, d)) for _ in range(2)] for _ in range(2)]
    for x in range(2):
        for y in range(2):
            if kind == 'temporal':
                pa = np.abs(A[x].conj().T @ v) ** 2
                T = np.abs(B[y].conj().T @ A[x]) ** 2          # T[b, a]
                P[x][y] = (pa[None, :] * T).T
            else:
                psi = v.reshape(d, d)
                P[x][y] = np.abs(A[x].conj().T @ psi @ B[y].conj()) ** 2
    return I_d(P, d)
def run(args):
    d, kind, seed = args
    n = 4 * d * d + (2 * d if kind == 'temporal' else 2 * d * d)
    r = minimize(lambda p: -value(p, d, kind), np.random.default_rng(seed).normal(size=n), method='L-BFGS-B',
                 options=dict(maxiter=6000, ftol=1e-15, gtol=1e-11))
    return d, kind, seed, -r.fun, r.x.tolist()
if __name__ == '__main__':
    jobs = [(d, kind, s) for d in (6, 7) for kind in ('temporal', 'spatial') for s in range(8)]
    with Pool(5) as p:
        res = p.map(run, jobs, chunksize=1)
    best = {}
    for d, kind, s, v, x in res:
        k = f"{d}|{kind}"
        if k not in best or v > best[k][0]: best[k] = (v, x, s)
    json.dump({k: dict(value=v[0], seed=v[2], params=v[1]) for k, v in best.items()}, open('tcgd67.json', 'w'))
    for k, v in sorted(best.items()): print(k, v[0], flush=True)
