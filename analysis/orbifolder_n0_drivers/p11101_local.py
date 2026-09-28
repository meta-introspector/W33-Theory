"""Pass 11101 draft: localized light Higgs.  For every model and every INDIVIDUAL Higgs doublet H_k, the up (down, lepton)
Yukawa matrix Y^(k): untwisted triples eps/delta structure (Pass 11098), twisted triples O(1) * eps^n_diff (n_diff = tori
where the three fixed points differ).  Singular values at eps = 1e-2, 1e-4 -> scaling exponents p_i.
usage: p11101_local.py our_v2.dump their_sm.dump out.json"""
import json
import sys
from collections import Counter
from fractions import Fraction as F

import numpy as np

sys.path.insert(0, r'C:\Repos\Theory of Everything\analysis')
import w33_pass11097_fractional_fermion_masses as P97  # noqa: E402
import w33_pass11098_yukawa_textures as Y  # noqa: E402
import w33_so16_selection_rules as S  # noqa: E402


def exps(entries, n1, n2, rng):
    """entries: (i, j, coeff_kind, s, nd). returns singular-value exponents in eps"""
    coef = {}
    for (i, j, kind, s, nd) in entries:
        c = s if kind == 'U' else (rng.standard_normal() + 1j * rng.standard_normal())
        coef[(i, j, nd)] = coef.get((i, j, nd), 0) + c
    out = []
    W = (1.0, 1.7, 2.9)          # anisotropic tori: eps_t = eps^W_t (generic, distinct)
    for eps in (1e-3, 1e-6):
        M = np.zeros((n1, n2), dtype=complex)
        for (i, j, nd), c in coef.items():
            M[i, j] += c * eps ** sum(w * d for w, d in zip(W, nd))
        sv = np.sort(np.linalg.svd(M, compute_uv=False))[::-1]
        out.append(sv)
    a, b = out
    p = []
    for x, y in zip(a, b):
        if x < 1e-300 or y < 1e-300:
            p.append(None)
        else:
            p.append(round(float(np.log(y / x) / np.log(1e-3)), 2))
    return p


def analyse(us, th, rng):
    algs, hid, nu1 = P97.prepare(us, th)
    c, oc, ow, qcol = S.sm_data(th, us)
    qb = P97.conj(qcol, 'A2')
    ferm = [f for f in us['fields'] if f['m'] == 2]
    scal = [f for f in us['fields'] if f['m'] == 6]
    pick = lambda L, col, w, Yv: [f for f in L if f['col'] == col and f['w'] == w and f['Y'] == Yv]
    sets = dict(up=(pick(ferm, qcol, 2, F(1, 6)), pick(ferm, qb, 1, F(-2, 3)), pick(scal, '1', 2, F(1, 2))),
                down=(pick(ferm, qcol, 2, F(1, 6)), pick(ferm, qb, 1, F(1, 3)), pick(scal, '1', 2, F(-1, 2))),
                lepton=(pick(ferm, '1', 2, F(-1, 2)), pick(ferm, '1', 1, F(1)), pick(scal, '1', 2, F(-1, 2))))
    out = {}
    for name, (A, B, H) in sets.items():
        perH = []
        for k, h in enumerate(H):
            ent = []
            for i, a in enumerate(A):
                for j, b in enumerate(B):
                    if not (P97.cubic_ok([S.charge_vector(x, True) for x in (a, b, h)], nu1) and Y.hidden_ok([a, b, h], hid, algs)):
                        continue
                    if all(x['k'] == 0 and x['l'] == 0 for x in (a, b, h)):
                        s = Y.eps(Y.plane(a), Y.plane(b), Y.plane(h))
                        if s:
                            ent.append((i, j, 'U', s, (0, 0, 0)))
                    else:
                        ca, cb, ch = S.classes(a['n']), S.classes(b['n']), S.classes(h['n'])
                        ent.append((i, j, 'T', None, tuple(int(not (ca[t] == cb[t] == ch[t])) for t in range(3))))
            if ent:
                perH.append(dict(k=k, sector=(h['k'], h['l']), entries=len(ent), exponents=exps(ent, len(A), len(B), rng)))
        out[name] = perH
    return out


def main(ours, theirs, outp):
    rng = np.random.default_rng(3)
    US, TH = S.parse_ours(ours), S.parse_theirs(theirs)
    res = {}
    for i in sorted(US):
        r = analyse(US[i], TH[i], rng)
        res[i] = dict(label=US[i]['label'], **r)
        best = [h['exponents'] for h in r['up']]
        print(i, 'up per-Higgs exponents:', best[:6], flush=True)
    json.dump(res, open(outp, 'w'), indent=1, default=str)


if __name__ == '__main__':
    main(*sys.argv[1:4])
