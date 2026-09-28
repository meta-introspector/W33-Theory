#!/usr/bin/env python3
"""Pass 11109: at the one-loop-stabilised family modulus T* = rho the tree-level Yukawas give m_c = m_u = m_t / 2 -- the
heavy top needs Im T* ~ 3.2, where the potential is higher.

In the 12 survivors (Passes 11101-11104) the three generations sit at the three theta-fixed points of the Wilson-line-free
family torus, the light Higgs at one of them.  Tree-level twisted Yukawas factorise over the tori; in the family torus the
coupling is Y_same when all three fields share the fixed point (the top) and Y_dist when the three points are distinct
(the charm-up block, off-diagonal).  The other tori contribute a common factor.  So

    m_c = m_u,   m_{c,u} / m_t = |Y_dist(T*) / Y_same(T*)|      (tree level; Kahler normalisations cancel).

Worldsheet instantons: the classical solution maps the sphere onto TWO copies of the triangle spanned by the three fixed
points (Schwarz-Christoffel), action = 2 x area / (2 pi alpha').  In the Narain normalisation of the one-loop engine
(self-dual circle G = 1, SU(3) point T = rho; torus area = 4 pi^2 alpha' Im T) this makes the couplings the A2 theta
functions of the Kahler modulus,

    Y_c(T) = sum_{v in A2 + c} exp(i pi T |v|^2)    (|root|^2 = 2;  c = 0: same point, c = f, 2f: distinct points),

the known result (Lauer, Mas, Nilles 1989-91; Kobayashi et al. 2018, arXiv:1804.06644).  Check: the triple (Y_0, Y_f, Y_2f)
closes under T -> -1/T with weight (-iT) and the finite Fourier matrix (1/sqrt3)[w^(jk)] -- the single-triangle
normalisation does not.  The minimal distinct-point triangle is 1/6 of the torus cell (checked).

Results:
  * at T* = rho (the minimum of the one-loop potential in the family direction, Pass 11106): |Y_dist / Y_same| = 1/2
    exactly (rho is the fixed point of S T; the coupling vector is an eigenvector of the Fourier matrix);
  * the observed m_c/m_t = 0.0036 (M_Z) needs Im T* = 3.21 (0.0027 at a high scale: 3.35);
  * |Y_dist/Y_same| at i, 1.5i, 2i, 3i: 0.366, 0.130, 0.045, 0.0056.
The potential at Im T* = 3.2 (Pass 11110) is higher than at rho: the stabilised vacuum has no family hierarchy.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11109_instanton_yukawas_at_rho.json"
E1 = np.array([np.sqrt(2), 0.0])
E2 = np.array([-1 / np.sqrt(2), np.sqrt(1.5)])
FP = (2 * E1 + E2) / 3
W3 = np.exp(2j * np.pi / 3)
RHO = 0.5 + 0.8660254037844386j


def theta_a2(T, shift, k=2, N=14):
    """sum_{v in A2 + shift} exp(2 pi i T |v|^2 / k); k = 2 is the double-triangle (physical) normalisation"""
    a = np.arange(-N, N + 1)
    A, B = np.meshgrid(a, a, indexing='ij')
    v = shift[None, None, :] + A[..., None] * E1 + B[..., None] * E2
    return np.sum(np.exp(2j * np.pi * T * np.einsum('ijk,ijk->ij', v, v) / k))


def duality_closure(k, T=0.13 + 1.21j):
    vec = lambda t: np.array([theta_a2(t, c * FP, k) for c in range(3)])
    M = np.array([[W3 ** (i * j) for j in range(3)] for i in range(3)]) / np.sqrt(3)
    return bool(np.allclose(vec(-1 / T), (-1j * T) * (M @ vec(T)), rtol=1e-8))


def ratio(T):
    return abs(theta_a2(T, FP) / theta_a2(T, 0 * FP))


def main():
    res = dict(pass_id=11109,
               minimal_triangle_area_fraction=float(FP @ FP / 4),
               duality_closes_double_triangle=duality_closure(2), duality_closes_single_triangle=duality_closure(4),
               ratio_at_rho=ratio(RHO), ratio_at_rho_minus_1=ratio(RHO - 1),
               ratio_on_axis={str(y): ratio(1j * y) for y in (1.0, 1.5, 2.0, 3.0, 4.0)},
               ImT_for_mc_over_mt={str(t): brentq(lambda y: ratio(1j * y) - t, 1.0, 20.0) for t in (0.0036, 0.0027)})
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    main()
