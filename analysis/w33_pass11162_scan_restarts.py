"""Pass 11162 scan: 96 random restarts (Nelder-Mead then Powell, random scales) of N_AB + N_AC over three-qutrit pure
states; frozen output data/w33_pass11162_restarts.json (about 3 minutes on 4 cores)."""
import json

import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize


def neg(rho, d=3):
    pt = rho.reshape(d, d, d, d).transpose(0, 3, 2, 1).reshape(d * d, d * d)
    ev = np.linalg.eigvalsh((pt + pt.conj().T) / 2)
    return float((np.abs(ev).sum() - 1) / 2)


def pair(psi, keep):
    T = psi.reshape(3, 3, 3)
    M = T.reshape(9, 3) if keep == (0, 1) else T.transpose(0, 2, 1).reshape(9, 3)
    return M @ M.conj().T


def f(p, which):
    v = p[:27] + 1j * p[27:]
    v /= np.linalg.norm(v)
    return -sum(neg(pair(v, k)) for k in which)


def run(seed):
    x0 = np.random.default_rng(5000 + seed).normal(size=54) * np.random.default_rng(seed).uniform(0.2, 3)
    r = minimize(f, x0, args=(((0, 1), (0, 2)),), method='Nelder-Mead', options=dict(maxiter=50000, maxfev=50000, xatol=1e-11, fatol=1e-13))
    r = minimize(f, r.x, args=(((0, 1), (0, 2)),), method='Powell', options=dict(maxiter=30000, xtol=1e-11, ftol=1e-14))
    return -r.fun


if __name__ == '__main__':
    with Pool(4) as p:
        vals = p.map(run, range(96))
    json.dump(dict(values=vals, best=max(vals), target=4 / np.sqrt(15)), open('restarts.json', 'w'), indent=1)
    print(max(vals), 4 / np.sqrt(15), sum(v > 4 / np.sqrt(15) + 1e-9 for v in vals), flush=True)
