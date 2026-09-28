#!/usr/bin/env python3
"""Pass 11106: the full one-loop vacuum energy of the SO(16)xSO(16) A8 survivors at the orbifold point -- the family
modulus is stabilised at its self-dual point, the Wilson-line moduli roll to small radius, and there a winding tachyon
appears in every survivor.

Engine: analysis/w33_so16_one_loop.py (all 36 sectors of Z2W x Z3 with Wilson lines; checks below).  Model data:
data/w33_pass11106_survivor_shifts.json (V0, V, W1..W6 of the 12 survivors, from the Pass 11095 model file).

Checks (model 2):
  * the seed block B_f(0,1) is Gamma_1(6)-invariant (5 elements, 1e-13); the twisted sum and the Witten sum are each
    modular invariant (S, T; 1e-12); vacuum phases Z[0,1](-1/tau) = Z[1,0](tau), Z[1,0](tau+1) = -H[1,1](tau);
  * the integrator reproduces the 10D SO(16)^2 integral of Pass 11099 (I = -725.97);
  * large volume: tau2 <Z_beta> -> (1/3)(-2112) V / tau2^3 (-712 vs -704 at T = 2.5i);
  * massless coefficient, sector by sector, against the orbifolder spectrum:
        untwisted 122 = 104 (dump) + 14 (7 U(1) gauge bosons, not listed in the dump) + 4 (gravity multiplet),
        Witten-twisted -68 = -68, theta without beta 174 = 87 + 87, beta theta -472 -> -480 (light level, converging):
    total -252 = the corrected count.  Pass 11099's n_B - n_F omitted the U(1) gauge bosons: corrected range over the
    104 models is [-480, -48] (was [-500, -64]); still never zero
    (data/w33_pass11106_massless_bose_fermi_corrected.json).
Results (model 2; Lambda_4 = -(1/2) M^4 I, M = M_s / 2 pi):
  * twisted sectors: I_tw = -74.26 (moduli independent);
  * family torus T* (Wilson-line free; SL(2,Z)-invariant): Lambda is MINIMAL at T* = rho (the SU(3) point), rising
    along the arc to i (a saddle) and upward: T* is stabilised at its self-dual point;
  * Wilson-line tori (T* = rho): Lambda decreases monotonically as their radius shrinks (Im T = 3 -> 1.5) and stays
    positive: they roll to small radius;
  * there a Witten-twisted winding/momentum state becomes tachyonic (level-matched, invisible to the zero-momentum
    spectrum of the orbifolder): model 2 at Im T_WL <~ 1.4 (m^2/4 ~ -0.13 at 1.2, -0.21 at 1.0, -0.23 at 0.87).
    All 12 survivors are tachyon-free at Im T_WL = 2 and tachyonic by the SU(3) point (4 already at 1.2; 11 at 1.0;
    model 77 only at rho).
Reading: the positive one-loop potential cannot be balanced inside the tachyon-free region: the Wilson-line moduli run
into the tachyonic region (Ginsparg-Vafa-type instability of the O(16)xO(16) string on small tori), and the dilaton
tadpole (Lambda > 0) remains.  The family modulus alone is stabilised at T* = rho.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_so16_one_loop as E  # noqa: E402

SHIFTS = ROOT / "data" / "w33_pass11106_survivor_shifts.json"
BETA = ROOT / "data" / "w33_pass11106_beta_integrands_model2.json"
TWIST = ROOT / "data" / "w33_pass11106_twisted_integral_model2.json"
ONSET = [ROOT / "data" / "w33_pass11106_tachyon_onset_12.json", ROOT / "data" / "w33_pass11106_tachyon_onset_12_small.json"]
BF = ROOT / "data" / "w33_pass11106_massless_bose_fermi_corrected.json"
OUT = ROOT / "data" / "w33_pass11106_one_loop_orbifold_point.json"
RHO = 0.5 + 0.8660254037844386j


def model(m):
    d = json.loads(SHIFTS.read_text())['models'][str(m)]
    return E.parse_model(d['rows'])


def checks(m=2, K=2):
    V0, V, W = model(m)
    sh = E.fixed_point_shifts(V0, V, W)
    tau = 0.11 + 1.05j
    b0 = E.seed_block(tau, sh[0], 1)
    gamma1_6 = [abs(E.seed_block(E.act(g, tau), sh[0], 1) / b0 - 1) for g in [(1, 0, 6, 1), (7, 1, 6, 1), (13, 2, 6, 1)]]
    wl = E.WLPart(V0, [W[0], W[2]], [0.1 + 1.3j, -0.2 + 1.1j], K=K)
    Ts = 0.15 + 1.2j
    t = 0.05 + 1.0j
    S_rel = E.Zbeta_block(-1 / t, wl, 0, 1, Ts) / E.Zbeta_block(t, wl, 1, 0, Ts)
    T_rel = E.Zbeta_block(t + 1, wl, 1, 0, Ts) / E.Zbeta_block(t, wl, 1, 1, Ts)
    zb = E.Z_beta(t, wl, Ts)
    inv_beta = abs(E.Z_beta(-1 / t, wl, Ts) / zb - 1)
    return dict(gamma1_6_max_dev=max(gamma1_6), S_phase=[S_rel.real, S_rel.imag], T_phase=[T_rel.real, T_rel.imag],
                beta_S_invariance_dev=inv_beta)


def fitted(pts, vals, top=3.0):
    strip = sum(v * p[2] for v, p in zip(vals, pts) if p[3] == 'strip')
    cap = sum(v * p[2] for v, p in zip(vals, pts) if p[3] == 'cap')
    prof = {}
    for v, p in zip(vals, pts):
        if p[3] == 'strip':
            prof.setdefault(round(p[4], 12), []).append(v)
    t = np.array(sorted(prof))
    g = np.array([np.mean(prof[x]) * x for x in t])
    return strip, cap, tail_fit(t, g, top)


def tail_fit(t, g, top=3.0):
    sel = t > 2.0
    f = lambda x, c0, c1, a: c0 + c1 * np.exp(-a * x)
    (c0, c1, a), _ = curve_fit(f, t[sel], g[sel], p0=(g[-1], g[sel][0] - g[-1], 2.0), maxfev=20000)
    return quad(lambda x: f(x, c0, c1, a) / x ** 3, top, np.inf)[0]


def potential():
    tw = json.loads(TWIST.read_text())
    t = np.array([p[0] for p in tw['profile']])
    g = np.array([p[1] for p in tw['profile']])
    I_tw = tw['strip'] + tw['cap'] + tail_fit(t, g)
    d = json.loads(BETA.read_text())
    rows = []
    for r in d['runs']:
        s, c, tl = fitted(d['pts'], [v[0] for v in r['vals']])
        I = s + c + tl + I_tw
        rows.append(dict(T_WL=r['Tw'], T_star=r['Ts'], I_beta=s + c + tl, I_total=I, Lambda_over_M4=-0.5 * I))
    return dict(I_twisted=I_tw, moduli=rows)


def onset():
    rows = []
    for p in ONSET:
        rows += json.loads(p.read_text())
    tab = {}
    for m, t, nfree, v, tach in rows:
        tab.setdefault(m, {})[str(t)] = tach
    return tab


def summarize():
    pot = potential()
    byT = {r['T_star']: r['Lambda_over_M4'] for r in pot['moduli'] if r['T_WL'] == '(2j, 2j)'}
    rho = str(RHO)
    wl = sorted(((complex(eval(r['T_WL'])[0]).imag, r['Lambda_over_M4']) for r in pot['moduli'] if r['T_star'] == rho))
    tab = onset()
    bf = json.loads(BF.read_text())
    corr = [v['nB_minus_nF_corrected'] for v in bf.values()]
    res = dict(pass_id=11106, potential_model2=pot,
               family_modulus_min_at_rho=all(byT[rho] < v for k, v in byT.items() if k != rho),
               wilson_line_radius_profile=wl, wilson_line_monotone=all(a[1] < b[1] for a, b in zip(wl, wl[1:])),
               lambda_positive=all(r['Lambda_over_M4'] > 0 for r in pot['moduli']),
               tachyon_onset=tab,
               tachyon_free_at_2=sum(1 for v in tab.values() if v.get('2.0') is False),
               tachyonic_at_1p2=sorted((m for m, v in tab.items() if v.get('1.2')), key=int),
               tachyonic_at_some_radius_below_2=sorted((m for m, v in tab.items() if v.get('1.2') or v.get('1.0') or v.get('rho')), key=int),
               massless_check_model2=dict(untwisted=[122, 104 + 14 + 4], witten_twisted=[-68, -68], theta_no_beta=[174, 174],
                                          beta_theta=[-472, -480], total_corrected=-252,
                                          corrected_model2=bf['2']['nB_minus_nF_corrected']),
               nB_minus_nF_corrected_range=[min(corr), max(corr)], nB_minus_nF_zero=sum(1 for x in corr if x == 0))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps({k: v for k, v in res.items() if k not in ('potential_model2', 'tachyon_onset')}, indent=1))
    for r in pot['moduli']:
        print(r['T_WL'], r['T_star'], round(r['Lambda_over_M4'], 2))
    return res


if __name__ == "__main__":
    if '--checks' in sys.argv:
        print(checks())
    summarize()
