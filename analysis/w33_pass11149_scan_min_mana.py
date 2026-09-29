"""Pass 11149 scan: minimum mana of a two-qutrit state reaching CGLMP value 2 + eps with the optimal Pauli-measurement Bell operator (40 Nelder-Mead + Powell restarts per eps); frozen as data/w33_pass11149_min_mana_curve.json"""
import sys; sys.path.insert(0, r'C:\Repos\Theory of Everything\analysis')
import numpy as np, itertools, json
from scipy.optimize import minimize
from multiprocessing import Pool
import w33_pass11149_magic_for_spookiness as P
AOPS = np.array([np.kron(P.A1[u], P.A1[v]) for u in P.PTS for v in P.PTS])   # (81, 9, 9)
def mana_fast(v):
    w = np.real(np.einsum('i,uij,j->u', v.conj(), AOPS, v)) / 9
    return float(np.log(np.abs(w).sum()))
def bestB():
    PB = P.pauli_bases(); zb = PB[0]; best = None
    for a1, b0, b1 in itertools.product(range(len(PB)), repeat=3):
        B = P.bell_operator(zb, PB[a1], PB[b0], PB[b1]); B = (B + B.conj().T) / 2
        ev = np.linalg.eigvalsh(B)
        if best is None or ev[-1] > best[0] + 1e-12: best = (ev[-1], B)
    return best[1]
def run(eps):
    B = bestB()
    def st(p):
        v = p[:9] + 1j * p[9:]; return v / np.linalg.norm(v)
    def f(p):
        v = st(p); val = float(np.real(v.conj() @ B @ v))
        return mana_fast(v) + 500 * max(0.0, 2 + eps - val) ** 2
    rng = np.random.default_rng(int(eps * 1000) + 7); bestm = None
    for s in range(40):
        r = minimize(f, rng.normal(size=18), method='Nelder-Mead', options=dict(maxiter=40000, xatol=1e-10, fatol=1e-12))
        r = minimize(f, r.x, method='Powell', options=dict(maxiter=20000, xtol=1e-10, ftol=1e-13))
        v = st(r.x); val = float(np.real(v.conj() @ B @ v)); m = mana_fast(v)
        if val >= 2 + eps - 2e-3 and (bestm is None or m < bestm[0]): bestm = (m, val)
    return eps, bestm
if __name__ == '__main__':
    # check mana_fast against the reference implementation
    v = np.random.default_rng(0).normal(size=9) + 1j * np.random.default_rng(1).normal(size=9); v /= np.linalg.norm(v)
    assert abs(mana_fast(v) - P.mana(np.outer(v, v.conj()), 2)) < 1e-10
    with Pool(4) as p:
        res = p.map(run, [0.02, 0.1, 0.3, 0.6])
    json.dump({str(e): b for e, b in res}, open('minmana2.json', 'w'), indent=1)
    for e, b in res: print(e, b, flush=True)
