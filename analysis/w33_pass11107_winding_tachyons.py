#!/usr/bin/env python3
"""Pass 11107: the winding tachyons of the SO(16)xSO(16) A8 survivors -- explicit spectrum, critical radii and charges.

Pass 11106 saw, in the one-loop integrand, a level-matched state that becomes tachyonic when the Wilson-line tori shrink.
Here it is built explicitly.  Witten-twisted sector, right mover r = -v0 (q_R = 0, the O-class ground state), no
oscillators:
    Delta = -1/2 + pR^2 / 2,   level matching  l^2 + pL^2 - pR^2 = 1,   l = pi + V0 + A n  (pi in E8 x E8),
    beta projection: exp 2 pi i (pi.V0 + V0^2) = -1   (the sign fixed by Z[1,1] = -H[1,1] in the one-loop engine),
with p_{L,R} of all three tori (Narain vectors of analysis/w33_so16_one_loop.py) and the family torus at T* = rho.  The
Z3 image of a state with torus momentum/winding is another state, so one combination per orbit survives.

Checks: at Im T_WL = 1.0 the lowest level Delta = -0.211 matches the growth rate of the one-loop integrand (-0.205);
no tachyon at large radius (T_WL = T* = 5i: the 10D SO(16)^2 string is tachyon-free); none at Im T_WL = 2 in any survivor.
Results (both Wilson-line tori at T_WL = i y, T* = rho):
  * critical radii y_c: sqrt3 = 1.73 (models 2, 14, 53, 102), 1.51 (10, 13, 15, 35, 57, 69, 78); model 77 is
    tachyon-free on the whole imaginary axis and tachyonic at the SU(3) point;
  * CORRECTION to Pass 11106: its onset "Im T_WL ~ 1.4" came from a growth-ratio test that misses small |Delta|; the
    enumerated onset is 1.51-1.73.  Its Lambda value at T_WL = 1.5i (model 2) lies in the tachyonic region and is void;
    the monotone fall of Lambda holds on the tachyon-free segment 1.75 <= Im T_WL <= 3;
  * the tachyons are colour singlets; in 10/12 models every tachyonic state has fractional electric charge
    (+-1/3, +-2/3; e.g. SU(2) doublets with Y = +-1/6), so any condensate breaks U(1)_EM; in models 53 and 57 some
    degenerate tachyons are fully SM-neutral (Y = 0, SU(2) singlet), so a neutral condensation direction exists there.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_so16_one_loop as E  # noqa: E402

D = json.loads((ROOT / "data" / "w33_pass11106_survivor_shifts.json").read_text())['models']
CRIT = ROOT / "data" / "w33_pass11107_critical_radii_12.json"
QN = ROOT / "data" / "w33_pass11107_tachyon_quantum_numbers_12.json"
OUT = ROOT / "data" / "w33_pass11107_winding_tachyons.json"
RHO = 0.5 + 0.8660254037844386j


def e8_near(w, r2):
    """E8 vectors x with |x + w|^2 < r2"""
    out = []
    for half in (0.0, 0.5):
        cands = []
        for c in range(8):
            base = -w[c] - half
            lo, hi = int(np.floor(base - np.sqrt(r2))) - 1, int(np.ceil(base + np.sqrt(r2))) + 1
            cands.append([k + half for k in range(lo, hi + 1) if (k + half + w[c]) ** 2 < r2])
        for x in itertools.product(*cands):
            x = np.array(x)
            if round(x.sum()) % 2 == 0 and np.sum((x + w) ** 2) < r2 - 1e-12:
                out.append(x)
    return out


def torus_table(T, off, N, nmax=3, mmax=6):
    """{rounded (pL^2 - pR^2): min pR^2} over n in Z^2 with n1 + n2 = N (N=None: any n), m in Z^2; offsets off (2,)"""
    G, B = E.GB(T)
    Gi = np.linalg.inv(G)
    tab = {}
    rng = range(-nmax, nmax + 1)
    for n1 in rng:
        for n2 in rng:
            if N is not None and n1 + n2 != N:
                continue
            n = np.array([n1, n2], float)
            for m1 in range(-mmax, mmax + 1):
                for m2 in range(-mmax, mmax + 1):
                    base = np.array([m1, m2]) - off + B @ n
                    wp, wm = base + G @ n, base - G @ n
                    pl, pr = 0.5 * wp @ Gi @ wp, 0.5 * wm @ Gi @ wm
                    k = round(pl - pr, 6)
                    if pr < 1.2 and (k not in tab or pr < tab[k][0]):
                        tab[k] = (pr, (n1, n2, m1, m2))
    return tab


def tachyons(m, Tw, Ts, Nrange=3, rows=None):
    V0, V, W = E.parse_model(rows if rows is not None else D[str(m)]['rows'])
    tori = [0, 2, 4]
    wlt = [i for i in tori if np.any(W[i] != 0)]
    Wt = [W[i] for i in wlt]
    tabs_free = torus_table(Ts, np.zeros(2), None)
    found = []
    for N1 in range(-Nrange, Nrange + 1):
        for N2 in range(-Nrange, Nrange + 1):
            w = V0 + N1 * Wt[0] + N2 * Wt[1]
            A = [e8_near(w[:8], 2.0), e8_near(w[8:], 2.0)]
            for x1 in A[0]:
                for x2 in A[1]:
                    pi = np.concatenate([x1, x2])
                    l = pi + w
                    l2 = l @ l
                    if l2 >= 2 - 1e-9:
                        continue
                    ph = np.exp(2j * np.pi * (pi @ V0 + V0 @ V0))
                    if abs(ph + 1) > 1e-6:
                        continue                       # beta projection keeps phase -1
                    pp = pi + V0
                    offs = []
                    for t, (Nt, Wv) in enumerate(zip((N1, N2), Wt)):
                        o = (pp @ Wv + 0.5 * Wv @ (N1 * Wt[0] + N2 * Wt[1])) % 1.0
                        offs.append(np.array([o, o]))
                    t1 = torus_table(Tw[0], offs[0], N1)
                    t2 = torus_table(Tw[1], offs[1], N2)
                    target = round(1 - l2, 6)
                    best = None
                    for k1, (r1, s1) in t1.items():
                        for k2, (r2, s2) in t2.items():
                            k3 = round(target - k1 - k2, 6)
                            if k3 in tabs_free:
                                r3, s3 = tabs_free[k3]
                                tot = r1 + r2 + r3
                                if tot < 1 - 1e-9 and (best is None or tot < best[0]):
                                    best = (tot, s1, s2, s3)
                    if best:
                        found.append(dict(Delta=-0.5 + best[0] / 2, l=l.tolist(), l2=l2, N=(N1, N2), tori=best[1:]))
    found.sort(key=lambda x: x['Delta'])
    return found




def summarize():
    crit = json.loads(CRIT.read_text())
    qn = json.loads(QN.read_text())
    charges = {}
    for m, T, out in qn:
        Qs = sorted({q for o in out for q in o['Q']})
        neutral = [o for o in out if abs(o['Y']) < 1e-9 and all(abs(x) < 1e-9 for x in o['su2_weight'])]
        charges[m] = dict(T_WL=T, n_states=len(out), min_Delta=min(o['Delta'] for o in out), Q=Qs,
                          colour_singlets=all(all(abs(x) < 1e-9 for x in o['colour_weight']) for o in out),
                          sm_neutral_states=len(neutral),
                          all_fractional=all(abs(q * 3 % 3) > 1e-6 for o in out for q in o['Q']))
    res = dict(pass_id=11107, critical_radius_ImT_WL=crit, charges=charges,
               models_all_fractional=sorted((m for m, v in charges.items() if v['all_fractional']), key=int),
               models_with_neutral_tachyon=sorted((m for m, v in charges.items() if v['sm_neutral_states']), key=int))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps({k: v for k, v in res.items() if k != 'charges'}, indent=1))
    return res


if __name__ == "__main__":
    if '--check' in sys.argv:
        print('T_WL = 1.0i: min Delta', tachyons('2', [1j, 1j], RHO)[0]['Delta'], '| 2i:', len(tachyons('2', [2j, 2j], RHO)),
              '| 5i/5i:', len(tachyons('2', [5j, 5j], 5j)))
    summarize()
