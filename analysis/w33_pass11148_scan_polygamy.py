"""Pass 11148: spatial monogamy vs temporal polygamy of qutrit negativity"""
import numpy as np, json
from scipy.optimize import minimize
from multiprocessing import Pool
def neg(rho, d=3):
    pt = rho.reshape(d, d, d, d).transpose(0, 3, 2, 1).reshape(d * d, d * d)
    ev = np.linalg.eigvalsh((pt + pt.conj().T) / 2)
    return float((np.abs(ev).sum() - 1) / 2)
def pair(psi, keep):
    T = psi.reshape(3, 3, 3)
    if keep == (0, 1): M = T.reshape(9, 3); return M @ M.conj().T
    if keep == (0, 2): M = T.transpose(0, 2, 1).reshape(9, 3); return M @ M.conj().T
    M = T.transpose(1, 2, 0).reshape(9, 3); return M @ M.conj().T
def f(p, which):
    v = p[:27] + 1j * p[27:]; v /= np.linalg.norm(v)
    tot = sum(neg(pair(v, k)) for k in which)
    return -tot
def run(args):
    which, seed = args
    r = minimize(f, np.random.default_rng(seed).normal(size=54), args=(which,), method='Nelder-Mead',
                 options=dict(maxiter=40000, maxfev=40000, xatol=1e-10, fatol=1e-12))
    r = minimize(f, r.x, args=(which,), method='Powell', options=dict(maxiter=20000, xtol=1e-10, ftol=1e-13))
    return str(which), -r.fun
if __name__ == '__main__':
    AB_AC = ((0, 1), (0, 2)); ALL = ((0, 1), (0, 2), (1, 2))
    jobs = [(AB_AC, s) for s in range(16)] + [(ALL, s) for s in range(16)]
    with Pool(2) as p:
        res = p.map(run, jobs)
    best = {}
    for k, v in res: best[k] = max(best.get(k, -1), v)
    # antisymmetric state and GHZ as references
    A = np.zeros(27)
    import itertools
    for perm in itertools.permutations(range(3)):
        sgn = np.linalg.det(np.eye(3)[list(perm)])
        A[perm[0] * 9 + perm[1] * 3 + perm[2]] = sgn
    A /= np.linalg.norm(A)
    G = np.zeros(27); G[[0, 13, 26]] = 1; G /= np.linalg.norm(G)
    refs = dict(antisymmetric=[neg(pair(A, k)) for k in ALL], ghz=[neg(pair(G, k)) for k in ALL])
    # temporal: identity channel, maximally mixed start: every pair of times has PDM = SWAP/3
    SW = np.zeros((9, 9))
    for i in range(3):
        for j in range(3): SW[i * 3 + j, j * 3 + i] = 1
    tneg = float((np.abs(np.linalg.eigvalsh(SW / 3)).sum() - 1) / 2)
    out = dict(spatial_max=best, refs=refs, temporal_pair_negativity=tneg, temporal_AB_AC=2 * tneg, temporal_all=3 * tneg)
    json.dump(out, open('poly.json', 'w'), indent=1)
    print(json.dumps(out, indent=1))
