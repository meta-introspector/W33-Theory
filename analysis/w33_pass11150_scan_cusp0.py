"""Pass 11150 scan: beta-sector integrands for 40521 at Im T = 0.2, 0.17 (K = 5) and at their Fricke images (K = 3); frozen as data/w33_pass11150_cusp0_beta_integrands.json"""
import numpy as np, json, sys
from multiprocessing import Pool
sys.path.insert(0, r'C:\Repos\Theory of Everything\analysis')
import w33_so16_one_loop as E
R = json.load(open(r'C:\Repos\Theory of Everything\data\w33_pass11108_rescan_levels.json'))
RHO = 0.5 + 0.8660254037844386j
LAB = 'A8SM_20260971_40521_TF'
RUNS = [(0.2j, 5), (0.17j, 5), (-1 / (3 * 0.2j), 3), (-1 / (3 * 0.17j), 3)]
MODEL = {}
def model():
    if 'm' not in MODEL:
        V0, V, W = E.parse_model(R[LAB]['rows'])
        wlt = [i for i in (0, 2, 4) if np.any(W[i] != 0)]
        MODEL['m'] = (V0, [W[i] for i in wlt])
    return MODEL['m']
def f(args):
    k, tau = args
    T, K = RUNS[k]
    V0, Wt = model()
    return E.Z_beta(tau, E.WLPart(V0, Wt, [T, T], K=K), RHO)
if __name__ == '__main__':
    pts = E.grid()
    with Pool(3) as p:
        vals = p.map(f, [(k, q[0]) for k in range(len(RUNS)) for q in pts], chunksize=16)
    n = len(pts)
    out = [dict(T=[RUNS[k][0].real, RUNS[k][0].imag], K=RUNS[k][1], vals=[[x.real, x.imag] for x in vals[k * n:(k + 1) * n]]) for k in range(len(RUNS))]
    json.dump(dict(model=LAB, pts=[[q[0].real, q[0].imag, q[1], q[2], q[3]] for q in pts], runs=out), open('cusp0_runs.json', 'w'))
    print('done')
