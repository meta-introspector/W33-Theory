#!/usr/bin/env python3
"""Pass 11098: renormalizable Yukawa textures of the SO(16)xSO(16) W(3,3) A8 models -- no single heavy top.

Cubic couplings psi psi phi under the EXACT tree-level rules of analysis/w33_so16_selection_rules.py, for
up (q ubar H_u), down (q dbar H_d) and lepton (l ebar H_d) Yukawas, in all 104 models of Pass 11095.

Structure, not just existence:
  * untwisted triples: the coupling is the 10D gauge coupling with the internal epsilon structure,
    y_ij = g eps(p_i, p_j, p_k) h_k (p = complex plane of the untwisted field).  A 3x3 matrix sum_k eps_ijk h_k is
    antisymmetric; its singular values are (|h|, |h|, 0) for EVERY Higgs direction h -- two degenerate heavy states and
    one massless one;
  * twisted triples: O(1) when the three fields sit at the same fixed point in every torus, exponentially suppressed
    otherwise.
A random-coefficient rank (the first estimate) over-states the untwisted case (3 instead of 2) -- the trap recorded in
the project memory "structural rank is only an upper bound".
Frozen results for all 104 in data/w33_pass11098_yukawa_104.json; the sample models are recomputed here.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11097_fractional_fermion_masses as P97  # noqa: E402
import w33_so16_selection_rules as S  # noqa: E402

FROZEN = ROOT / "data" / "w33_pass11098_yukawa_104.json"
OUT = ROOT / "data" / "w33_pass11098_yukawa_textures.json"


def hidden_ok(fields, hid, algs):
    for j in hid:
        nontriv = [f['dim'][j] for f in fields if S.dimof(f['dim'][j]) != 1]
        if not nontriv:
            continue
        if len(nontriv) == 2 and nontriv[1] == P97.conj(nontriv[0], algs[j]):
            continue
        return False
    return True


def plane(f):
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
    return 1 if (a, b, c) in ((1, 2, 3), (2, 3, 1), (3, 1, 2)) else -1


def epsilon_singular_values(h):
    """singular values of M_ij = sum_k eps_ijk h_k (normalised): always (1, 1, 0)"""
    M = np.zeros((3, 3), dtype=complex)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                M[i, j] += eps(i + 1, j + 1, k + 1) * h[k]
    sv = np.linalg.svd(M, compute_uv=False)
    return sv / sv.max()


def textures(us, th, rng):
    algs, hid, nu1 = P97.prepare(us, th)
    ferm = [f for f in us['fields'] if f['m'] == 2]
    scal = [f for f in us['fields'] if f['m'] == 6]
    c, oc, ow, qcol = S.sm_data(th, us)
    qb = P97.conj(qcol, 'A2')
    pick = lambda L, col, w, Y: [f for f in L if f['col'] == col and f['w'] == w and f['Y'] == Y]
    sets = dict(up=(pick(ferm, qcol, 2, F(1, 6)), pick(ferm, qb, 1, F(-2, 3)), pick(scal, '1', 2, F(1, 2))),
                down=(pick(ferm, qcol, 2, F(1, 6)), pick(ferm, qb, 1, F(1, 3)), pick(scal, '1', 2, F(-1, 2))),
                lepton=(pick(ferm, '1', 2, F(-1, 2)), pick(ferm, '1', 1, F(1)), pick(scal, '1', 2, F(-1, 2))))
    out = {}
    for name, (A, B, H) in sets.items():
        trip = []
        for i, a in enumerate(A):
            for j, b in enumerate(B):
                for k, h in enumerate(H):
                    if not (P97.cubic_ok([S.charge_vector(x, True) for x in (a, b, h)], nu1) and hidden_ok([a, b, h], hid, algs)):
                        continue
                    if all(x['k'] == 0 and x['l'] == 0 for x in (a, b, h)):
                        trip.append((i, j, k, 'U', eps(plane(a), plane(b), plane(h)), 0))
                    else:
                        ca, cb, ch = S.classes(a['n']), S.classes(b['n']), S.classes(h['n'])
                        trip.append((i, j, k, 'T', None, sum(1 for t in range(3) if not (ca[t] == cb[t] == ch[t]))))
        hv = rng.standard_normal(len(H)) + 1j * rng.standard_normal(len(H)) if H else []
        M = np.zeros((len(A), len(B)), dtype=complex)
        for (i, j, k, kind, s, nd) in trip:
            if nd == 0:
                M[i, j] += (s if kind == 'U' else rng.standard_normal() + 1j * rng.standard_normal()) * hv[k]
        sv = np.linalg.svd(M, compute_uv=False) if M.size else np.array([])
        sv = sv / sv.max() if sv.size and sv.max() > 0 else sv
        out[name] = dict(triples=len(trip), untwisted=sum(1 for t in trip if t[3] == 'U'),
                         o1_rank=int(np.sum(sv > 1e-9)) if sv.size else 0, o1_singular_values=[round(float(x), 6) for x in sv[:3]])
    return out


def main():
    rng = np.random.default_rng(1)
    for _ in range(5):
        sv = epsilon_singular_values(rng.standard_normal(3) + 1j * rng.standard_normal(3))
        assert np.allclose(sorted(sv), [0, 1, 1])
    US, TH = P97.sample_models()
    frozen = json.loads(FROZEN.read_text())
    recheck = {}
    for i in US:
        t = textures(US[i], TH[i], rng)
        f = frozen[str(i)]
        recheck[str(i)] = dict(recomputed={k: (v['o1_rank'], v['untwisted']) for k, v in t.items()},
                               frozen={k: (f[k]['o1_rank'], f[k]['untwisted']) for k in ('up', 'down', 'lepton')})
        print(i, recheck[str(i)], flush=True)
    summ = {name: dict(Counter(str(v[name]['o1_rank']) for v in frozen.values())) for name in ('up', 'down', 'lepton')}
    deg = sum(1 for v in frozen.values() if v['up']['o1_rank'] == 2 and abs(v['up']['o1_singular_values'][1] - 1) < 1e-6)
    res = dict(pass_id=11098, models=len(frozen), o1_rank_distribution=summ,
               up_top_charm_degenerate_untwisted=deg,
               up_twisted_models=sum(1 for v in frozen.values() if v['up']['untwisted'] == 0),
               up_rank_exactly_one=sum(1 for v in frozen.values() if v['up']['o1_rank'] == 1),
               epsilon_theorem="singular values of sum_k eps_ijk h_k are (1,1,0) for every h (checked on random h)",
               sample_recheck=recheck)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True, default=str))
    print(json.dumps({k: v for k, v in res.items() if k != 'sample_recheck'}, indent=1))


if __name__ == "__main__":
    main()
