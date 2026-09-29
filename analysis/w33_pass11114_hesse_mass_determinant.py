#!/usr/bin/env python3
"""Pass 11114: the up-quark mass determinant is a Hesse cubic in the Higgs direction; at the stabilised family modulus
it is the singular (triangle) member -- and no Higgs alignment beats the localized Higgs: m_2/m_top >= 1/2 at T = rho.

Setup (the 12 survivors; Passes 11101-11109).  The three Higgs doublets H_k sit at the three fixed points of the
family torus, a Delta(54) triplet.  With a general light-Higgs direction h = (h1, h2, h3) the tree-level up-type matrix is

    M(h) = [[s h1, d h3, d h2], [d h3, s h2, d h1], [d h2, d h1, s h3]],   s = Y_same(T*), d = Y_dist(T*)

(Y = A2 theta functions of the family modulus, Pass 11109; same-point and all-distinct couplings only -- two equal
points and a third different is not collinear, hence forbidden).  Then

    det M(h) = h1 h2 h3 (s^3 + 2 d^3) - s d^2 (h1^3 + h2^3 + h3^3),

a member of the HESSE PENCIL x^3 + y^3 + z^3 - 3 lambda x y z with lambda(T) = (s^3 + 2 d^3) / (3 s d^2): the other
track's Hesse configuration (Passes 11061-11066, 2026-09-23 Hesse notes) is the geometry of the quark-mass determinant.
Membership is forced by Delta(27) (its invariant cubics ARE the Hesse pencil; Artebani-Dolgachev); the content is lambda(T).

Results:
  * at T* = rho (the one-loop minimum, Pass 11106/11110): d/s = e^(i pi/3)/2 and lambda = w^2 EXACTLY, so lambda^3 = 1:
    the determinant is a SINGULAR Hesse member, a triangle of three lines;
  * on the standard Delta(27) alignments the spectrum at rho is (1, 1/2, 1/2) [(1,0,0), (1,1,1), (1,w,w^2), ...] or
    (1, 1, 0) [(1,1,w), (1,w,1), (1,w^2,w^2)]: a massless up but top = charm;
  * rank 1 (one heavy, two massless) needs M_ij^2 = M_ii M_jj, i.e. (d/s)^6 = 1: impossible for |d/s| != 1;
  * min over all Higgs directions of sigma_2/sigma_1 = |d/s| (60 random starts x Nelder-Mead at each of rho, i, 1.5i),
    attained by the LOCALIZED Higgs: 1/2 at rho, (sqrt3-1)/2 at i.  No alignment can produce the hierarchy at the
    stabilised point.
  * T-dependence beyond one loop: T6/Z3 has no N = 2 subsectors, so the (DKL) gauge threshold corrections are
    independent of the Kahler moduli; hidden-sector condensates cannot move T* away from rho at this order.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11109_instanton_yukawas_at_rho as Y  # noqa: E402

OUT = ROOT / "data" / "w33_pass11114_hesse_mass_determinant.json"
W = np.exp(2j * np.pi / 3)


def mass_matrix(h, s, d):
    h1, h2, h3 = h
    return np.array([[s * h1, d * h3, d * h2], [d * h3, s * h2, d * h1], [d * h2, d * h1, s * h3]])


def couplings(T):
    return Y.theta_a2(T, 0 * Y.FP), Y.theta_a2(T, Y.FP)


def hesse_lambda(T):
    s, d = couplings(T)
    return (s ** 3 + 2 * d ** 3) / (3 * s * d * d)


def determinant_identity(T, trials=20, seed=1):
    s, d = couplings(T)
    rng = np.random.default_rng(seed)
    err = 0
    for _ in range(trials):
        h = rng.normal(size=3) + 1j * rng.normal(size=3)
        lhs = np.linalg.det(mass_matrix(h, s, d))
        rhs = h[0] * h[1] * h[2] * (s ** 3 + 2 * d ** 3) - s * d * d * (h ** 3).sum()
        err = max(err, abs(lhs - rhs) / max(1, abs(lhs)))
    return err


def alignments(T):
    s, d = couplings(T)
    out = {}
    for name, h in [('(1,0,0)', (1, 0, 0)), ('(1,1,1)', (1, 1, 1))] + \
            [(f'(1,w^{a},w^{b})', (1, W ** a, W ** b)) for a, b in itertools.product(range(3), repeat=2)]:
        h = np.array(h, complex)
        v = np.linalg.svd(mass_matrix(h / np.linalg.norm(h), s, d), compute_uv=False)
        out[name] = [round(float(x), 9) for x in v / v[0]]
    return out


def min_ratio(T, starts=60, seed=0):
    s, d = couplings(T)

    def f(x):
        h = np.array([1, x[0] + 1j * x[1], x[2] + 1j * x[3]])
        v = np.linalg.svd(mass_matrix(h, s, d), compute_uv=False)
        return v[1] / v[0]
    rng = np.random.default_rng(seed)
    best = min((minimize(f, rng.normal(size=4) * 2, method='Nelder-Mead',
                         options=dict(xatol=1e-12, fatol=1e-14, maxiter=20000)) for _ in range(starts)), key=lambda r: r.fun)
    return float(best.fun), float(abs(d / s))


def main(starts=60):
    lam = hesse_lambda(Y.RHO)
    res = dict(pass_id=11114, lambda_at_rho=[lam.real, lam.imag], lambda_cubed_at_rho=[(lam ** 3).real, (lam ** 3).imag],
               d_over_s_at_rho=[complex(couplings(Y.RHO)[1] / couplings(Y.RHO)[0]).real, complex(couplings(Y.RHO)[1] / couplings(Y.RHO)[0]).imag],
               determinant_identity_max_rel_err=determinant_identity(Y.RHO),
               alignments_at_rho=alignments(Y.RHO),
               min_sigma2_over_sigma1={str(T): min_ratio(T, starts) for T in (Y.RHO, 1j, 1.5j)})
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    main()
