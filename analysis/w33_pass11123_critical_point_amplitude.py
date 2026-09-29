#!/usr/bin/env python3
"""Pass 11123: the critical radius of the six neutral-exit survivors -- exact values sqrt3 and phi^2/sqrt3, no enhanced
gauge symmetry, and a tree-level four-point amplitude with no contact quartic: the condensate runs with the radion.

A. Exact critical radius (both Wilson-line tori at i y, B = 0, family torus at rho): solving pR^2(y) = 1 for the neutral
   state's own Narain vector gives y_c = sqrt3 = 1.732050808 (36621, 17224) and y_c = 1.511522628 = phi^2 / sqrt3
   (46043, 24165, 40521, 5904; phi the golden ratio: 3 y_c^2 = phi^4).  At y_c every neutral state has P_L^2 = 2,
   pR^2 = 1 exactly (a left current times a right weight-1/2 operator).
B. Not an enhanced-symmetry point: no pair of the massless states sums to an untwisted vector with pR = 0, P_L^2 = 2
   (the sums have (P_L^2, pR^2) = (0,0) [conjugates], (3,3), (5,1) for sqrt3; (phi^2, phi^2), (3 + 1/phi^2... =5.382,
   1.382) for phi^2/sqrt3): no new massless gauge bosons, so no D-term-like quartic and no moduli trapping of the
   Kofman-Linde-Liu-Maloney-McAllister-Silverstein type.
C. The tree-level four-point amplitude T(P) T(P) -> T(P) T(P) at y_c (alpha' = 2; two (-1)- and two (0)-picture vertices;
   right-moving fermion contraction k2.k3 + pR2.pR3; complex beta function of Kawai-Lewellen-Tye, checked numerically
   to 1e-10):
        A = C pi Gamma(3 - s/2) Gamma(-t/2) Gamma(-u/2) / [Gamma(s/2 - 1) Gamma(2 + t/2) Gamma(2 + u/2)],
   symmetric in t <-> u (Bose symmetry).  Low energy: A = 4 C pi (1/t + 1/u) + 3 C pi s^2/(t u) + O(E^2) -- pure massless
   exchange, NO analytic contact term: at tree level nothing stops the condensate at a small VEV.
D. The radion coupling: dDelta/dy at y_c = 1/(2 sqrt3) = 0.2887 (sqrt3 models), 0.2466 (phi models), positive: a
   condensate |T|^2 lowers the energy as the torus shrinks -- T and the radius run together; the endpoint lies beyond
   the perturbative neighbourhood of y_c (winding-tachyon condensation, cf. topology change).
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
PAIRS = ROOT / "data" / "w33_pass11123_critical_point_pairs.json"
RAD = ROOT / "data" / "w33_pass11123_radion_slopes.json"
OUT = ROOT / "data" / "w33_pass11123_critical_point_amplitude.json"
PHI = (1 + 5 ** 0.5) / 2


def klt_check():
    mp.mp.dps = 15
    a, c, b, d = mp.mpf(0.4), mp.mpf(-0.6), mp.mpf(-0.7), mp.mpf(-0.7)
    formula = mp.pi * mp.gamma(1 + a) * mp.gamma(1 + b) * mp.gamma(-1 - c - d) / (mp.gamma(-c) * mp.gamma(-d) * mp.gamma(2 + a + b))

    def f(r, th):
        z = r * mp.e ** (1j * th)
        return (z ** a * (1 - z) ** b * mp.conj(z) ** c * mp.conj(1 - z) ** d).real * r
    numeric = 2 * mp.quad(lambda r: mp.quad(lambda th: f(r, th), [0, mp.pi / 2, mp.pi]), [0, 0.5, 1, 2, mp.inf])
    return float(formula), float(numeric)


def low_energy():
    x, y, e, a, b = sp.symbols('x y e a b')
    F = sp.gamma(3 - x) * sp.gamma(-y) * sp.gamma(x + y) / (sp.gamma(x - 1) * sp.gamma(2 + y) * sp.gamma(2 - x - y))
    ser = sp.expand(sp.series(F.subs({x: e * a, y: e * b}), e, 0, 1).removeO())
    Fs = F.subs({x: e * a, y: e * b})
    lead = sp.simplify(sp.limit(e * Fs, e, 0))
    const = sp.simplify(sp.limit(Fs - lead / e, e, 0))
    symmetric = sp.simplify(F - F.subs(y, -x - y)) == 0
    return str(lead), str(const), bool(symmetric)


def summarize():
    pairs = json.loads(PAIRS.read_text())
    rad = json.loads(RAD.read_text())
    lead, const, sym = low_energy()
    f, n = klt_check()
    res = dict(pass_id=11123,
               critical_radii={k: v['y_c'] for k, v in pairs.items()},
               sqrt3=3 ** 0.5, phi2_over_sqrt3=PHI ** 2 / 3 ** 0.5,
               states_PL2_pR2={k: v['states'] for k, v in pairs.items()},
               massless_vector_from_pairs=any(p['massless_vector'] for v in pairs.values() for p in v['pairs']),
               amplitude_lead=lead, amplitude_const=const, t_u_symmetric=sym, klt_formula_vs_numeric=[f, n],
               radion_slopes={k: v['dDelta_dy'] for k, v in rad.items()})
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True, default=str))
    print(json.dumps({k: v for k, v in res.items() if k != 'states_PL2_pR2'}, indent=1, default=str))
    return res


if __name__ == "__main__":
    summarize()
