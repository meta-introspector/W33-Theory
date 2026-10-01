#!/usr/bin/env python3
"""Shared exact algebra for Passes 11255--11259.

The basic G26 invariants use the standard coordinates of the crystallographic
reflection representation.  They are not silently identified with the three
vectors in Pass 11218's regular abelian centralizer; Pass 11255 proves that
centralizer contains a nonzero nilpotent element and is not a Cartan subspace.
"""
from __future__ import annotations

import sympy as sp

X, Y, Z = sp.symbols("x y z")
VARS = (X, Y, Z)


def invariants():
    x, y, z = VARS
    u6 = x**6 + y**6 + z**6 - 10 * (x**3*y**3 + x**3*z**3 + y**3*z**3)
    u9 = (x**3-y**3) * (x**3-z**3) * (y**3-z**3)
    u12 = (
        x**12 + y**12 + z**12
        - 110 * sum(a**9*b**3 for a in VARS for b in VARS if a != b)
        + 462 * sum(VARS[i]**6 * VARS[j]**6 for i in range(3) for j in range(i+1, 3))
    )
    return sp.expand(u6), sp.expand(u12), sp.expand(u9**2)


def invariant_jacobian():
    fs = invariants()
    J = sp.Matrix([[sp.diff(f, v) for v in VARS] for f in fs])
    return J, sp.factor(J.det())


def target_potential_data(point=(1, 2, 3)):
    """Dimensionless orbit-distance potential and its exact Hessian at point.

    V=sum_i (u_i/u_i(p)-1)^2.  This is a controlled invariant potential,
    not a fitted physical scalar potential.  At a regular point its zero set
    is one closed G26 orbit and its Hessian is positive definite.
    """
    fs = invariants()
    sub = dict(zip(VARS, point))
    vals = tuple(sp.Integer(f.subs(sub)) for f in fs)
    assert all(vals)
    Jn = sp.Matrix([
        [sp.diff(f, v).subs(sub) / vals[i] for v in VARS]
        for i, f in enumerate(fs)
    ])
    H = sp.simplify(2 * Jn.T * Jn)
    return fs, vals, Jn, H


def rational_matrix(M):
    return [[str(sp.factor(x)) for x in row] for row in M.tolist()]

