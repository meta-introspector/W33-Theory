"""Pass 11134 scan: the zero-radius (T-dual) tachyon species -- l^2 = 1 vectors in E8xE8 + V0 + N W_neutral (N mod 3), beta phase -1, for the 21 neutral-exit survivors; frozen as data/w33_pass11134_dual_tachyon_multiplet_21.json"""
import sys, json
import numpy as np
from collections import Counter
sys.path.insert(0, r'C:\Repos\Theory of Everything\analysis')
import w33_pass11107_winding_tachyons as T7
import w33_so16_one_loop as E
R = r'C:\Repos\Theory of Everything\data'
tor = json.load(open(R + r'\w33_pass11119_torus_resolved_tachyons_491.json'))
labs = json.load(open(R + r'\w33_pass11119_neutral_tori_across_491.json'))['passing_gauntlet']
resc = json.load(open(R + r'\w33_pass11108_rescan_levels.json'))
sh = json.load(open(R + r'\w33_pass11106_survivor_shifts.json'))['models']
out = {}
for lab in labs:
    rows = resc[lab]['rows'] if lab in resc else next(v['rows'] for v in sh.values() if v['label'] == lab)
    V0, V, W = E.parse_model(rows)
    wlt = [i for i in (0, 2, 4) if np.any(W[i] != 0)]
    t = 0 if 'neutral' in tor[lab]['torus1']['kinds'] else 1
    Wn = W[wlt[t]]
    states = []
    for N in range(3):
        w = V0 + N * Wn
        for x1 in T7.e8_near(w[:8], 1.0 + 1e-6):
            for x2 in T7.e8_near(w[8:], 1.0 + 1e-6):
                pi = np.concatenate([x1, x2]); l = pi + w
                if abs(l @ l - 1) > 1e-9: continue
                if abs(np.exp(2j * np.pi * (pi @ V0 + V0 @ V0)) + 1) > 1e-6: continue
                states.append((N, tuple(np.round(l, 6)), round(float(l[:8] @ l[:8]), 6)))
    split = Counter((N, s) for N, _, s in states)
    out[lab] = dict(total=len(states), by_N={str(N): sum(1 for n, _, _ in states if n == N) for N in range(3)},
                    split={f"N{N}|l1sq={s}": c for (N, s), c in sorted(split.items())},
                    three_W_in_lattice=bool(np.allclose(((3 * Wn) * 2) % 1, 0)))
    print(lab, out[lab], flush=True)
json.dump(out, open('dualtach.json', 'w'), indent=1)
