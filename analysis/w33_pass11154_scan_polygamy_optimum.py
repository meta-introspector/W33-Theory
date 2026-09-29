"""Pass 11154 scan (poly.py = analysis/w33_pass11148_scan_polygamy.py in the working directory): 30-restart optimum of N_AB + N_AC over three-qutrit pure states and its invariants"""
import numpy as np, json
from scipy.optimize import minimize
exec(open('poly.py').read().split("if __name__ == '__main__':")[0])
AB_AC = ((0, 1), (0, 2))
best = None
for s in range(30):
    r = minimize(f, np.random.default_rng(1000 + s).normal(size=54), args=(AB_AC,), method='Nelder-Mead', options=dict(maxiter=60000, maxfev=60000, xatol=1e-12, fatol=1e-14))
    r = minimize(f, r.x, args=(AB_AC,), method='Powell', options=dict(maxiter=40000, xtol=1e-12, ftol=1e-15))
    if best is None or r.fun < best.fun: best = r
v = best.x[:27] + 1j * best.x[27:]; v /= np.linalg.norm(v)
print('value', -best.fun, 4 / np.sqrt(15))
T = v.reshape(3, 3, 3)
sA = np.linalg.svd(T.reshape(3, 9), compute_uv=False) ** 2
sB = np.linalg.svd(T.transpose(1, 0, 2).reshape(3, 9), compute_uv=False) ** 2
sC = np.linalg.svd(T.transpose(2, 0, 1).reshape(3, 9), compute_uv=False) ** 2
print('A spectrum', np.round(sA, 6), 'B', np.round(sB, 6), 'C', np.round(sC, 6))
for k in ((0, 1), (0, 2), (1, 2)):
    r_ = pair(v, k); pt = r_.reshape(3, 3, 3, 3).transpose(0, 3, 2, 1).reshape(9, 9)
    print(k, 'neg', round(neg(r_), 8), 'rho eig', np.round(np.linalg.eigvalsh(r_), 6), 'pt eig', np.round(np.linalg.eigvalsh((pt + pt.conj().T) / 2), 6))
np.save('poly_best.npy', v)
