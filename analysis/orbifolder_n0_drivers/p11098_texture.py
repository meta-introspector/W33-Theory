"""Pass 11098: physical structure of the cubic Yukawas (up, down, lepton) of the 104 SO(16)xSO(16) A8 models.
Untwisted triples: y_ij = g * eps(p_i, p_j, p_k) h_k (the internal gamma-matrix/epsilon structure of the 10D gauge
coupling; p = complex plane of the untwisted field) -> the matrix is ANTISYMMETRIC in the planes.
Twisted triples: O(1) if all three fields sit at the same fixed point in every torus, otherwise suppressed by
exp(-area) per torus where the fixed points differ (encoded as eps^n_diff, eps -> 0 for the leading order).
Reports the O(1) (unsuppressed) singular-value pattern for a random light-Higgs direction.
usage: p11098_texture.py our_v2.dump their_sm.dump out.json"""
import json
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, sys.path[0])
from p11098_yukawa import cubic_ok, hidden_ok  # noqa: E402
from p11097_frac import conj  # noqa: E402
from so16_couplings import F, classes, dimof, parse_ours, parse_theirs, sm_data, sm_info  # noqa: E402


def plane(f):
    """complex plane (1,2,3) of an untwisted field from its right-mover momentum, or 0 if none singled out"""
    q = f['qsh'][1:]
    if f['m'] == 6:
        nz = [a for a in range(3) if q[a] != 0]
        return nz[0] + 1 if len(nz) == 1 else 0
    signs = [x > 0 for x in q]
    if signs.count(True) == 2:
        return signs.index(False) + 1
    if signs.count(False) == 2:
        return signs.index(True) + 1
    return 0


def eps(a, b, c):
    """coefficient of an untwisted psi_a psi_b phi_c Yukawa from the 10D gauge coupling: the internal SU(4) product
    4 x 4 -> 6.  Fermion planes 1..3 are the SU(3) triplet weights, plane 0 the SU(3) singlet (-,-,-); the scalar
    plane is its 6-weight.  triplet.triplet.scalar: epsilon_ijk; singlet.triplet_i.scalar_j: delta_ij."""
    if a == 0 and b != 0:
        return 1 if b == c else 0
    if b == 0 and a != 0:
        return 1 if a == c else 0
    if len({a, b, c}) < 3 or 0 in (a, b, c):
        return 0
    return {(1, 2, 3): 1, (2, 3, 1): 1, (3, 1, 2): 1, (1, 3, 2): -1, (3, 2, 1): -1, (2, 1, 3): -1}[(a, b, c)]


def analyse(us, th, rng):
    c, oc, ow, qcol = sm_data(th, us)
    algs = [a for a, _ in us['factors']]
    hid = [j for j in range(len(algs)) if j not in (oc, ow)]
    nu1 = len(us['u1'])
    ferm = [f for f in us['fields'] if f['m'] == 2]
    scal = [f for f in us['fields'] if f['m'] == 6]
    for f in ferm + scal:
        f['Y'], f['col'], f['w'], f['t'], f['frac'] = sm_info(f, c, oc, ow, qcol)
    qb = conj(qcol, 'A2')
    pick = lambda L, col, w, Y: [f for f in L if f['col'] == col and f['w'] == w and f['Y'] == Y]
    Q, U, D = pick(ferm, qcol, 2, F(1, 6)), pick(ferm, qb, 1, F(-2, 3)), pick(ferm, qb, 1, F(1, 3))
    Lp, E = pick(ferm, '1', 2, F(-1, 2)), pick(ferm, '1', 1, F(1))
    Hu, Hd = pick(scal, '1', 2, F(1, 2)), pick(scal, '1', 2, F(-1, 2))
    out = {}
    for name, A, B, H in (('up', Q, U, Hu), ('down', Q, D, Hd), ('lepton', Lp, E, Hd)):
        trip = []
        for i, a in enumerate(A):
            for j, b in enumerate(B):
                for k, h in enumerate(H):
                    if not (cubic_ok(a, b, h, nu1) and hidden_ok([a, b, h], hid, algs)):
                        continue
                    untw = all(x['k'] == 0 and x['l'] == 0 for x in (a, b, h))
                    if untw:
                        s = eps(plane(a), plane(b), plane(h))
                        trip.append((i, j, k, 'U', s, 0))
                    else:
                        ca, cb, ch = classes(a['n']), classes(b['n']), classes(h['n'])
                        nd = sum(1 for t in range(3) if not (ca[t] == cb[t] == ch[t]))
                        trip.append((i, j, k, 'T', None, nd))
        # O(1) matrix at a random light-Higgs direction
        hvec = rng.standard_normal(len(H)) + 1j * rng.standard_normal(len(H)) if H else []
        M = np.zeros((len(A), len(B)), dtype=complex)
        for (i, j, k, kind, s, nd) in trip:
            if nd:
                continue
            coeff = s if kind == 'U' else (rng.standard_normal() + 1j * rng.standard_normal())
            M[i, j] += coeff * hvec[k]
        sv = np.linalg.svd(M, compute_uv=False) if M.size else np.array([])
        sv = sv / sv.max() if sv.size and sv.max() > 0 else sv
        out[name] = dict(triples=len(trip), untwisted=sum(1 for t in trip if t[3] == 'U'),
                         untwisted_zero_eps=sum(1 for t in trip if t[3] == 'U' and t[4] == 0),
                         twisted_same_fixed_point=sum(1 for t in trip if t[3] == 'T' and t[5] == 0),
                         twisted_suppressed=sum(1 for t in trip if t[3] == 'T' and t[5] > 0),
                         o1_rank=int(np.sum(sv > 1e-9)) if sv.size else 0,
                         o1_singular_values=[round(float(x), 6) for x in sv[:3]])
    return out


def main(ours, theirs, outp):
    rng = np.random.default_rng(20260928)
    US = parse_ours(ours)
    TH = parse_theirs(theirs)
    res = {}
    for i in sorted(US):
        r = analyse(US[i], TH[i], rng)
        res[i] = dict(label=US[i]['label'], **r)
        print(i, 'up', r['up']['o1_rank'], r['up']['o1_singular_values'], 'down', r['down']['o1_rank'], 'lep', r['lepton']['o1_rank'], flush=True)
    json.dump(res, open(outp, 'w'), indent=1, default=str)
    for name in ('up', 'down', 'lepton'):
        print(name, 'O(1) rank distribution', dict(Counter(v[name]['o1_rank'] for v in res.values())))
    deg = sum(1 for v in res.values() if v['up']['o1_rank'] == 2 and abs(v['up']['o1_singular_values'][1] - 1) < 1e-6)
    print('up: O(1) rank 2 with two EQUAL singular values (top = charm degenerate):', deg)
    print('up: O(1) rank exactly 1 (heavy top only):', sum(1 for v in res.values() if v['up']['o1_rank'] == 1))


if __name__ == '__main__':
    main(*sys.argv[1:4])
