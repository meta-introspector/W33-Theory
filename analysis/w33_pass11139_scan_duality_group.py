"""Pass 11139: which modular maps of the Wilson-line modulus T leave the beta-sector partition function invariant"""
import sys, json
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, r'C:\Repos\Theory of Everything\analysis')
import w33_so16_one_loop as E
R = r'C:\Repos\Theory of Everything\data'
resc = json.load(open(R + r'\w33_pass11108_rescan_levels.json'))
RHO = 0.5 + 0.8660254037844386j
MODELS = ['A8SM_20260942_36621_TF', 'A8SM_20260971_40521_TF']
MAPS = {'T+1': (1, 1, 0, 1), 'T+3': (1, 3, 0, 1), 'S': (0, -1, 1, 0), 'Fricke3': 'F',
        'G0(3):T/(3T+1)': (1, 0, 3, 1), 'G^0(3):T/(T+1)': (1, 0, 1, 1), 'ell3:(-T+1)/(-3T+2)': (-1, 1, -3, 2)}
POINTS = {'T+1': [0.21 + 1.3j], 'T+3': [0.21 + 1.3j], 'S': [0.2 + 0.98j, 0.35 + 0.95j], 'Fricke3': [0.1 + 0.57j, 0.05 + 0.58j],
          'G0(3):T/(3T+1)': [-0.33 + 0.33j, -0.3 + 0.36j], 'G^0(3):T/(T+1)': [-0.5 + 0.9j, -0.45 + 0.95j],
          'ell3:(-T+1)/(-3T+2)': [0.55 + 0.3j, 0.45 + 0.31j]}
TAUS = [0.13 + 1.1j, -0.27 + 1.45j, 0.4 + 0.95j]
def act(m, T):
    if m == 'F': return -1 / (3 * T)
    a, b, c, d = m; return (a * T + b) / (c * T + d)
def model(lab):
    V0, V, W = E.parse_model(resc[lab]['rows'])
    wlt = [i for i in (0, 2, 4) if np.any(W[i] != 0)]
    return V0, [W[i] for i in wlt]
def z(args):
    lab, T, tau, K = args
    V0, Wt = model(lab)
    wl = E.WLPart(V0, Wt, [T, T], K=K)
    return E.Z_beta(tau, wl, RHO)
if __name__ == '__main__':
    jobs, keys = [], []
    for lab in MODELS:
        for name, m in MAPS.items():
            for T in POINTS[name]:
                for tau in TAUS:
                    for K in (3, 4):
                        jobs.append((lab, T, tau, K)); keys.append((lab, name, T, tau, K, 'T'))
                        jobs.append((lab, act(m, T), tau, K)); keys.append((lab, name, T, tau, K, 'gT'))
    with Pool(3) as p:
        vals = p.map(z, jobs, chunksize=2)
    out = {}
    for k, v in zip(keys, vals):
        lab, name, T, tau, K, which = k
        out.setdefault(lab, {}).setdefault(name, {}).setdefault(f"{T}|{tau}|K{K}", {})[which] = [v.real, v.imag]
    json.dump(out, open('dualgrp.json', 'w'), indent=1)
    for lab, d in out.items():
        for name, e in d.items():
            devs = []; kconv = []
            for key, v in e.items():
                zT = complex(*v['T']); zg = complex(*v['gT'])
                devs.append((key.split('|')[2], abs(zg / zT - 1)))
            print(lab[14:], name, 'K3 max rel dev', max(d_ for k_, d_ in devs if k_ == 'K3'), 'K4', max(d_ for k_, d_ in devs if k_ == 'K4'), flush=True)
