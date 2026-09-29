"""Pass 11142 search: two-doublet mixtures h = Hu_k + eps^a Hu_j (k tree-level, a = 1..3, geometric factor g = 1 or 0); frozen as data/w33_pass11142_mixing_search.json"""
import json, sys, numpy as np
from collections import Counter
sys.path.insert(0, r'C:\Repos\Theory of Everything\analysis')
import w33_pass11102_higher_order_quark_hierarchy as H
d = json.load(open("mix.json"))
def mat(e):
    return ({tuple(map(int, k.split(','))): v for k, v in e['orders'].items()}, {tuple(map(int, k.split(','))): v for k, v in e['geo'].items()})
def comb(e1, e2, rows, cols, a, g, rng):
    n1, g1 = mat(e1); n2, g2 = mat(e2); n, geo = {}, {}
    for key in set(n1) | set(n2):
        c = [x for x in ((n1[key] + (g * g1.get(key, 0) if n1[key] == 0 else 0)) if key in n1 else None,
                         (a + n2[key] + (g * g2.get(key, 0) if n2[key] == 0 else 0)) if key in n2 else None) if x is not None]
        n[key] = min(c); geo[key] = 0
    return H.exponents(n, geo, rows, cols, 1.0, rng)
out = {}
for lab, v in d.items():
    sz = v['sizes']; tree = [k for k in range(9) if any(o == 0 for o in mat(v['up'][str(k)])[0].values())]
    rng = np.random.default_rng(0); pats = Counter(); hits = []
    for g in (1.0, 0.0):
        for k in tree:
            for j in range(9):
                if j == k: continue
                for a in (1, 2, 3):
                    up = comb(v['up'][str(k)], v['up'][str(j)], sz['Q'], sz['U'], a, g, rng)
                    pats[(g, str(up))] += 1
                    if up and len(set(up)) == 3 and up[0] == 0: hits.append((g, k, j, a, up))
    out[lab] = dict(tree=tree, patterns={f"g={g}|{p}": c for (g, p), c in pats.items()}, hierarchical_hits=hits)
    print(lab[14:], dict(pats), 'hits', len(hits), flush=True)
json.dump(out, open('mixsearch.json', 'w'), indent=1)
