"""Pass 11160 scan: temporal vs spatial CGLMP for d = 8..10 with a batched evaluator and central-difference gradients (5 restarts each); frozen as data/w33_pass11160_cglmp_optima_d8_10.json"""
import numpy as np, json, sys
from scipy.optimize import minimize
from multiprocessing import Pool
def herm(p, d):
    """batch (m, d*d) -> Hermitian (m, d, d)"""
    m = p.shape[0]; H = np.zeros((m, d, d), complex); iu = np.triu_indices(d, 1); k = len(iu[0])
    H[:, iu[0], iu[1]] = p[:, :k] + 1j * p[:, k:2 * k]
    H = H + np.conj(np.transpose(H, (0, 2, 1)))
    H[:, np.arange(d), np.arange(d)] = p[:, 2 * k:2 * k + d]
    return H
def unit(p, d):
    w, V = np.linalg.eigh(herm(p, d))
    return V @ (np.exp(1j * w)[:, :, None] * np.conj(np.transpose(V, (0, 2, 1))))
def coeff(d):
    """C[x][y] (d x d): I = sum_{x,y,a,b} C[x][y][a,b] P(a,b|x,y)"""
    C = [[np.zeros((d, d)) for _ in range(2)] for _ in range(2)]
    for k in range(d // 2):
        c = 1 - 2 * k / (d - 1)
        for a in range(d):
            for b in range(d):
                if a == (b + k) % d: C[0][0][a, b] += c
                if b == (a + k + 1) % d: C[1][0][a, b] += c
                if a == (b + k) % d: C[1][1][a, b] += c
                if b == (a + k) % d: C[0][1][a, b] += c
                if a == (b - k - 1) % d: C[0][0][a, b] -= c
                if b == (a - k) % d: C[1][0][a, b] -= c
                if a == (b - k - 1) % d: C[1][1][a, b] -= c
                if b == (a - k - 1) % d: C[0][1][a, b] -= c
    return C
def batch_value(P, d, kind, C):
    n = d * d
    Us = [unit(P[:, i * n:(i + 1) * n], d) for i in range(4)]
    A, B = Us[:2], Us[2:]
    tot = np.zeros(P.shape[0])
    if kind == 'temporal':
        v = P[:, 4 * n:4 * n + d] + 1j * P[:, 4 * n + d:4 * n + 2 * d]; v /= np.linalg.norm(v, axis=1, keepdims=True)
        for x in range(2):
            pa = np.abs(np.einsum('mia,mi->ma', np.conj(A[x]), v)) ** 2          # (m, a)
            for y in range(2):
                T = np.abs(np.einsum('mib,mia->mab', np.conj(B[y]), A[x])) ** 2   # (m, a, b)
                tot += np.einsum('ma,mab,ab->m', pa, T, C[x][y])
    else:
        v = P[:, 4 * n:4 * n + n] + 1j * P[:, 4 * n + n:4 * n + 2 * n]; v /= np.linalg.norm(v, axis=1, keepdims=True)
        psi = v.reshape(-1, d, d)
        for x in range(2):
            for y in range(2):
                amp = np.einsum('mia,mij,mjb->mab', np.conj(A[x]), psi, np.conj(B[y]))
                tot += np.einsum('mab,ab->m', np.abs(amp) ** 2, C[x][y])
    return tot
def run(args):
    d, kind, seed = args
    C = coeff(d); n = 4 * d * d + (2 * d if kind == 'temporal' else 2 * d * d)
    h = 1e-6
    def fg(p):
        E = np.eye(n) * h
        P = np.vstack([p[None, :], p + E, p - E])
        vals = batch_value(P, d, kind, C)
        return -vals[0], -(vals[1:n + 1] - vals[n + 1:]) / (2 * h)
    r = minimize(fg, np.random.default_rng(seed).normal(size=n), jac=True, method='L-BFGS-B', options=dict(maxiter=4000, ftol=1e-14, gtol=1e-9))
    return d, kind, seed, -r.fun
if __name__ == '__main__':
    # method check against the committed d = 3 values
    jobs = [(3, 'temporal', 0), (3, 'spatial', 0)] + [(d, kind, s) for d in (8, 9, 10) for kind in ('temporal', 'spatial') for s in range(5)]
    with Pool(6) as p:
        res = p.map(run, jobs, chunksize=1)
    best = {}
    for d, kind, s, v in res:
        k = f"{d}|{kind}"; best[k] = max(best.get(k, -9), v)
    json.dump(best, open('tcgd810.json', 'w'), indent=1)
    for k, v in sorted(best.items()): print(k, v, flush=True)
