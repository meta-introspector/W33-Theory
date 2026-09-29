"""Pass 11140 scan: massless vectors at zero radius of the neutral torus (l^2 = 2, Wilson-line offset integral, three phase rules compared); frozen as data/w33_pass11140_zero_radius_roots_21.json"""
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
def cands(V0, Wn):
    out = []
    for N in range(3):
        w = N * Wn
        for x1 in T7.e8_near(w[:8], 2.0 + 1e-6):
            for x2 in T7.e8_near(w[8:], 2.0 + 1e-6):
                pi = np.concatenate([x1, x2]); l = pi + w
                if abs(l @ l - 2) < 1e-9: out.append((N, pi, l))
    return out
RULES = {'a: pi.V0 in Z': lambda N, pi, l, V0, Wn: pi @ V0,
         'b: pi.V0 + N W.V0/2 in Z': lambda N, pi, l, V0, Wn: pi @ V0 + 0.5 * N * (Wn @ V0),
         'c: l.V0 in Z': lambda N, pi, l, V0, Wn: l @ V0}
def closed(roots):
    S = {tuple(np.round(r, 6)) for r in roots}
    arr = np.array(roots)
    bad = 0
    for a in arr[:: max(1, len(arr) // 60)]:
        for b in arr:
            refl = b - (a @ b) * a
            if tuple(np.round(refl, 6)) not in S: bad += 1
    return bad
def classify(roots):
    """rank and number of roots per irreducible component via the connectivity of the root graph (a.b != 0)"""
    arr = np.array(roots); n = len(arr); G = np.abs(arr @ arr.T) > 1e-9
    seen = [False] * n; comps = []
    for i in range(n):
        if seen[i]: continue
        st = [i]; seen[i] = True; comp = []
        while st:
            a = st.pop(); comp.append(a)
            for b in np.nonzero(G[a])[0]:
                if not seen[b]: seen[b] = True; st.append(b)
        comps.append((int(np.linalg.matrix_rank(arr[comp])), len(comp)))
    return sorted(comps, reverse=True)
out = {}
for lab in labs:
    rows = resc[lab]['rows'] if lab in resc else next(v['rows'] for v in sh.values() if v['label'] == lab)
    V0, V, W = E.parse_model(rows)
    wlt = [i for i in (0, 2, 4) if np.any(W[i] != 0)]
    t = 0 if 'neutral' in tor[lab]['torus1']['kinds'] else 1
    Wn, Wo = W[wlt[t]], W[wlt[1 - t]]
    C = cands(V0, Wn)
    res = {}
    for name, rule in RULES.items():
        for off in (False, True):
            roots = [l for N, pi, l in C if abs(((rule(N, pi, l, V0, Wn)) + 1e-9) % 1) < 1e-6
                     and abs(((pi @ Wn + 0.5 * Wn @ (N * Wn)) + 1e-9) % 1) < 1e-6
                     and (not off or abs(((pi @ Wo + 0.5 * Wo @ (N * Wn)) + 1e-9) % 1) < 1e-6)]
            res[f"{name}|offset_int={off}"] = dict(n=len(roots), not_closed=closed(roots), comps=classify(roots))
    out[lab] = res
    print(lab, {k.split(':')[0] + ('|off2' if k.endswith('True') else ''): (v['n'], v['not_closed'], v['comps']) for k, v in res.items() if k.startswith('a')}, flush=True)

json.dump(out, open('roots3.json', 'w'), indent=1)
