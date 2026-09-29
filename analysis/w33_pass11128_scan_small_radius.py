"""Pass 11128 scan (needs the Pass 11108 orbifolder dumps): the neutral torus shrunk toward zero (other WL torus at 3i, family at rho) -- tachyon depth, count, charges; frozen as data/w33_pass11128_small_radius_six.json"""
import sys, json, re, math
import numpy as np
from functools import lru_cache
from multiprocessing import Pool
sys.path.insert(0, r'C:\Repos\Theory of Everything\analysis')
import w33_so16_selection_rules as S
import w33_pass11107_winding_tachyons as T7
import w33_so16_one_loop as E
R = r'C:\Repos\Theory of Everything\data'
D2 = (r'C:/Users/wiljd/AppData/Local/Temp/rescan2/all387_our.dump', r'C:/Users/wiljd/AppData/Local/Temp/rescan2/all387_theirs.dump')
SIX = {'A8SM_20260942_36621_TF': 0, 'A8SM_20260942_46043_TF': 1, 'A8SM_20260951_24165_TF': 1, 'A8SM_20260971_40521_TF': 0,
       'A8SM_20260971_5904_TF': 0, 'A8SM_20260972_17224_TF': 0}
YS = [1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2]
@lru_cache(maxsize=None)
def table(T, o, N, nmax, mmax):
    return T7.torus_table(T, np.array([o, o]), N, nmax, mmax)
def tachyons(rows, Tw, Nr):
    V0, V, W = E.parse_model(rows)
    wlt = [i for i in (0, 2, 4) if np.any(W[i] != 0)]
    Wt = [W[i] for i in wlt]
    free = table(T7.RHO, 0.0, None, 3, 6)
    found = []
    for N1 in range(-Nr[0], Nr[0] + 1):
        for N2 in range(-Nr[1], Nr[1] + 1):
            w = V0 + N1 * Wt[0] + N2 * Wt[1]
            A = [T7.e8_near(w[:8], 2.0), T7.e8_near(w[8:], 2.0)]
            for x1 in A[0]:
                for x2 in A[1]:
                    pi = np.concatenate([x1, x2]); l = pi + w; l2 = l @ l
                    if l2 >= 2 - 1e-9: continue
                    if abs(np.exp(2j * np.pi * (pi @ V0 + V0 @ V0)) + 1) > 1e-6: continue
                    pp = pi + V0
                    offs = [round(float((pp @ Wv + 0.5 * Wv @ (N1 * Wt[0] + N2 * Wt[1])) % 1.0), 9) for Wv in Wt]
                    t1 = table(complex(Tw[0]), offs[0], N1, max(3, Nr[0] + 2), 6)
                    t2 = table(complex(Tw[1]), offs[1], N2, max(3, Nr[1] + 2), 6)
                    target = round(1 - l2, 6); best = None
                    for k1, (r1, s1) in t1.items():
                        for k2, (r2, s2) in t2.items():
                            k3 = round(target - k1 - k2, 6)
                            if k3 in free:
                                tot = r1 + r2 + free[k3][0]
                                if tot < 1 - 1e-9 and (best is None or tot < best): best = tot
                    if best is not None:
                        found.append(dict(Delta=-0.5 + best / 2, l=l.tolist()))
    return found
def job(a):
    lab, idx, rows, t, y = a
    Tw = [y * 1j, 3j] if t == 0 else [3j, y * 1j]
    n = max(3, math.ceil(2.5 / y))
    Nr = (n, 3) if t == 0 else (3, n)
    return lab, y, tachyons(rows, Tw, Nr)
if __name__ == '__main__':
    resc = json.load(open(R + r'\w33_pass11108_rescan_levels.json'))
    US, TH = S.parse_ours(D2[0]), S.parse_theirs(D2[1])
    idx = {lab: int(i) for i, lab in (re.match(r'MODEL (\d+) (\S+)', l).groups() for l in open(D2[0]) if l.startswith('MODEL'))}
    jobs = [(lab, idx[lab], resc[lab]['rows'], t, y) for lab, t in SIX.items() for y in YS]
    with Pool(4) as p:
        res = p.map(job, jobs, chunksize=1)
    out = {}
    for lab, y, found in res:
        us, th = US[idx[lab]], TH[idx[lab]]
        c, oc, ow, qcol = S.sm_data(th, us)
        U = np.array([[float(x) for x in u] for u in us['u1']])
        fac = [np.array([[float(x) for x in r] for r in roots]) for _, roots in us['factors']]
        cl = {}
        for r in found:
            l = np.array(r['l']); Y = float(sum(float(ci) * qi for ci, qi in zip(c, U @ l)))
            neu = abs(Y) < 1e-9 and np.all(np.abs(fac[ow] @ l) < 1e-9) and np.all(np.abs(fac[oc] @ l) < 1e-9)
            k = 'neutral' if neu else 'charged'
            cl.setdefault(k, [0, 0.0]); cl[k][0] += 1; cl[k][1] = min(cl[k][1], r['Delta'])
        out.setdefault(lab, {})[str(y)] = cl
        print(lab[14:], y, cl, flush=True)
    json.dump(out, open('smallr.json', 'w'), indent=1)
