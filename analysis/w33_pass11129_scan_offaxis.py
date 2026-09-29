"""Pass 11129 scan: beta-sector integrands off the B = 0 axis (x = 0.25, 0.5) for 36621 and 40521 at y_c + 0.07 and 2.0; frozen as data/w33_pass11129_offaxis_beta_integrands.json"""
import numpy as np, json, sys
from multiprocessing import Pool
sys.path.insert(0, r'C:\Repos\Theory of Everything\analysis')
import w33_so16_one_loop as E
R = json.load(open(r'C:\Repos\Theory of Everything\data\w33_pass11108_rescan_levels.json'))
CRIT = json.load(open(r'C:\Repos\Theory of Everything\data\w33_pass11119_diagonal_onset_21.json'))
SIX = ['A8SM_20260942_36621_TF', 'A8SM_20260942_46043_TF', 'A8SM_20260951_24165_TF', 'A8SM_20260971_40521_TF',
       'A8SM_20260971_5904_TF', 'A8SM_20260972_17224_TF']
RHO = 0.5 + 0.8660254037844386j
RUNS = [(lab, x, y) for lab in ('A8SM_20260942_36621_TF', 'A8SM_20260971_40521_TF') for y in (round(CRIT[lab]['crit'] + 0.07, 3), 2.0) for x in (0.25, 0.5)]
MODELS = {}
def model(lab):
    if lab not in MODELS:
        V0, V, W = E.parse_model(R[lab]['rows'])
        wlt = [i for i in (0, 2, 4) if np.any(W[i] != 0)]
        MODELS[lab] = (V0, [W[i] for i in wlt])
    return MODELS[lab]
def f(args):
    k, tau = args
    lab, x, y = RUNS[k]
    V0, Wt = model(lab)
    wl = E.WLPart(V0, Wt, [x + y * 1j, x + y * 1j], K=3)
    return E.Z_beta(tau, wl, RHO)
if __name__ == '__main__':
    pts = E.grid()
    with Pool(8) as p:
        vals = p.map(f, [(k, q[0]) for k in range(len(RUNS)) for q in pts], chunksize=16)
    n = len(pts)
    out = [dict(model=RUNS[k][0], x=RUNS[k][1], y=RUNS[k][2], vals=[[x.real, x.imag] for x in vals[k * n:(k + 1) * n]]) for k in range(len(RUNS))]
    json.dump(dict(pts=[[q[0].real, q[0].imag, q[1], q[2], q[3]] for q in pts], runs=out), open('offaxis_runs.json', 'w'))
    print('done')
