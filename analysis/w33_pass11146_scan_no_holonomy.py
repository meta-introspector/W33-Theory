"""Pass 11146 scan (needs smallr.py = analysis/w33_pass11128_scan_small_radius.py in the working directory): Witten-sector tachyons on the diagonal of both Wilson-line tori with and without the Wilson lines, 21 survivors; frozen as data/w33_pass11146_no_holonomy_21.json"""
import sys, json, math
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, r'C:\Repos\Theory of Everything\analysis')
exec(open('smallr.py').read().split("def job(a):")[0])
labs = json.load(open(R + r'\w33_pass11119_neutral_tori_across_491.json'))['passing_gauntlet']
resc = json.load(open(R + r'\w33_pass11108_rescan_levels.json'))
sh = json.load(open(R + r'\w33_pass11106_survivor_shifts.json'))['models']
YS = [1.5, 1.0, 0.6, 0.3, 0.2]
def rows_of(lab):
    return resc[lab]['rows'] if lab in resc else next(v['rows'] for v in sh.values() if v['label'] == lab)
def tach_W(V0, Wt, Tw, Nr):
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
                    target = round(1 - l2, 6)
                    for k1, (r1, s1) in t1.items():
                        for k2, (r2, s2) in t2.items():
                            k3 = round(target - k1 - k2, 6)
                            if k3 in free and r1 + r2 + free[k3][0] < 1 - 1e-9:
                                found.append(-0.5 + (r1 + r2 + free[k3][0]) / 2)
    return found
def job(lab):
    V0, V, W = E.parse_model(rows_of(lab))
    wlt = [i for i in (0, 2, 4) if np.any(W[i] != 0)]
    out = {}
    for y in YS:
        n = max(3, math.ceil(2.5 / y))
        on = tach_W(V0, [W[i] for i in wlt], [y * 1j, y * 1j], (n, n))
        off = tach_W(V0, [0 * W[i] for i in wlt], [y * 1j, y * 1j], (n, n))
        out[str(y)] = dict(with_wl=len(on), without_wl=len(off), min_with=min(on) if on else None)
    return lab, out
if __name__ == '__main__':
    with Pool(2) as p:
        res = dict(p.map(job, labs, chunksize=1))
    json.dump(res, open('nohol.json', 'w'), indent=1)
    print('any tachyon without WL:', any(v[y]['without_wl'] for v in res.values() for y in v), flush=True)
    for k, v in res.items(): print(k, v, flush=True)
