#!/usr/bin/env python3
"""Pass 11121: bridge -- the other track's Hesse / H27 / Clifford / 243 objects, object by object, in the explicit
SO(16)xSO(16) W(3,3) string vacua (Passes 11105, 11109, 11111, 11112, 11114, 11116, 11117).

Computed here: the curve of quark-mass determinants.  With the Delta(54) Higgs triplet h the up-type mass determinant is
the Hesse cubic x^3 + y^3 + z^3 - 3 lambda(T) xyz with lambda(T) = (s^3 + 2 d^3)/(3 s d^2), s, d the A2 theta functions of
the family modulus (Passes 11109, 11114).  Its j-invariant j_curve = 27 L (L + 8)^3 / (L - 1)^3, L = lambda^3:
  * T = rho (fixed by ST, the one-loop minimum): lambda = w^2, L = 1 -- a SINGULAR member (triangle), j_curve = infinity;
  * T = i (fixed by S): lambda = 1 + sqrt3 exactly, L = 10 + 6 sqrt3, j_curve = 1728 = j(i) -- the order-4 automorphism
    the fixed point forces;
  * T -> i infinity: lambda -> infinity (the cusp; the fourth singular member).
Away from the fixed points j_curve(T) is NOT j(T), j(3T) or j(T/3) (checked at 1.2i, 1.5i, 2i, 0.3 + 1.1i): the quark-mass
curve is not a modular copy of the family torus; its special values are those the fixed points force.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11109_instanton_yukawas_at_rho as Y  # noqa: E402

OUT = ROOT / "data" / "w33_pass11121_hesse_gauge_bridge.json"

BRIDGE = [
    dict(other='external-A2 qutrit H27 (2026-09-21_physical_external_a2_h27.md)',
         ours='family Delta(54) of every model: X = fixed-point translation, Z = space-group phase (Pass 11105)',
         status='realised'),
    dict(other='H27 centre = FI Z3 exp(2 pi i Q_psi / 3) (2026-09-21_physical_fi_is_h27_center.md)',
         ours='centre = orbifold twist; a gauge element in 104/104 models; exactly the anomalous-U(1) Z3 in models 4, 25, 29 (Pass 11116)',
         status='realised as a gauge element; literally the FI Z3 in 3/104'),
    dict(other='Clifford-648 = H27:SL(2,3) = W33 point stabiliser (2026-09-21_physical_a2_clifford648_w33_bridge.md; 1296 with similitudes)',
         ours='only H27:<-I> (index 12); the quadratic phases are forbidden by the space-group rule (Pass 11105)',
         status='partially realised (-I only)'),
    dict(other='two-qutrit 3^(1+4) = H27_int o H27_ext, commutation graph W(3,3) (2026-09-21_e8_trinification_two_qutrit_pauli243.md)',
         ours='not on generations (faithful degree 3^k: k = 1), not on generation x colour, not on hidden SU(3) (Passes 11111, 11116)',
         status='not realised on any matter'),
    dict(other='matter 81 = 9 x H9 under 3^(1+4) (2026-09-21_e8_matter81_pauli243_restriction.md)',
         ours='requires the E6 trinification centre; the vacua carry SM x U(1)^n x hidden, no E6',
         status='not realised'),
    dict(other='Hesse configuration / Hesse pencil (Passes 11061-11066; 2026-09-23 Hesse notes; Artebani-Dolgachev)',
         ours='the up-quark mass determinant over the Delta(54) Higgs triplet is a Hesse-pencil cubic (forced by Delta(27)); the family modulus picks the member: singular at the one-loop minimum rho (Pass 11114)',
         status='realised: the Hesse pencil is the geometry of the quark masses'),
    dict(other='AG(3,3) / collinearity (27 fixed points of T6/Z3)',
         ours='theta^3 couplings = collinear triples; theta x theta^2 x untwisted = same point (Pass 11105)',
         status='realised'),
]


def couplings(T):
    return Y.theta_a2(T, 0 * Y.FP), Y.theta_a2(T, Y.FP)


def lam(T):
    s, d = couplings(T)
    return (s ** 3 + 2 * d ** 3) / (3 * s * d * d)


def jcurve(T):
    L = lam(T) ** 3
    return 27 * L * (L + 8) ** 3 / (L - 1) ** 3


def main():
    li = lam(1j)
    res = dict(pass_id=11121, bridge=BRIDGE,
               lambda_rho=[lam(Y.RHO).real, lam(Y.RHO).imag], lambda_rho_cubed=[(lam(Y.RHO) ** 3).real, (lam(Y.RHO) ** 3).imag],
               lambda_i=[li.real, li.imag], lambda_i_is_1_plus_sqrt3=bool(abs(li - (1 + 3 ** 0.5)) < 1e-10),
               jcurve_i=[jcurve(1j).real, jcurve(1j).imag],
               lambda_large_T=[float(abs(lam(4j))), float(abs(lam(6j)))])
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps({k: v for k, v in res.items() if k != 'bridge'}, indent=1))
    return res


if __name__ == "__main__":
    main()
