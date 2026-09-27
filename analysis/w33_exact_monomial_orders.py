#!/usr/bin/env python3
"""Exact lowest order of a monomial with prescribed charges (all orders, no cap, no tolerance).

An independent, tolerance-free certificate for the mixed-integer programs of Passes 10978-10980 (which carry the
Z_N congruences in integer slack variables); Pass 11024 checks those programs against it.

SIGN CONVENTION (the one trap): a coupling  X_1 ... X_m  S^e  is allowed iff  sum_i e_i v_i = w - sum_j v(X_j),
where w is the superpotential's charge vector (R-charges; zero for non-R components).  Use `coupling(...)`, which
builds that target; `order(target)` takes the raw right-hand side.

Problem.  Fields with charge vectors v_1..v_n in Q^d, the first nq components continuous (U(1)s) and the rest
discrete (defined mod 1).  Given a target t, find e in Z_{>=0}^n with sum e_i v_i = t on the continuous
components and = t mod 1 on the discrete ones, minimising sum e_i.

Method.  With den the common denominator, the integer system [den*Vc 0; den*Vd -den*I] (e, m) = den*t is solved
by Smith normal form: all integer solutions are e = e0 + N z, z in Z^k, k = dim of the kernel of the continuous
charges.  The kernel basis is LLL-reduced and e0 Babai-recentred in exact integer arithmetic, so the remaining
problem -- min sum(e0 + N z) subject to e0 + N z >= 0 -- is a k-dimensional integer program with small entries and
no congruences.  Every returned solution is re-verified in exact rational arithmetic; 'None' means the integer
program is proven infeasible (or the lattice system has no integer solution at all); a solver that does not
finish raises instead of returning None.
"""
from __future__ import annotations

from fractions import Fraction
from math import lcm

import numpy as np
import sympy as sp
from scipy.optimize import Bounds, LinearConstraint, milp
from sympy.matrices.normalforms import smith_normal_decomp
from sympy.polys.matrices import DomainMatrix


def _dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def lll_reduce(basis, delta=Fraction(3, 4)):
    """textbook LLL in exact rational arithmetic (integer basis vectors as rows; small rank)."""
    b = [list(v) for v in basis]
    k = len(b)
    if k <= 1:
        return b

    def gso():
        bs, mu = [], [[Fraction(0)] * k for _ in range(k)]
        for i in range(k):
            v = [Fraction(x) for x in b[i]]
            for j in range(i):
                mu[i][j] = Fraction(_dot(b[i], bs[j])) / _dot(bs[j], bs[j])
                v = [x - mu[i][j] * y for x, y in zip(v, bs[j])]
            bs.append(v)
        return bs, mu
    bs, mu = gso()
    i = 1
    while i < k:
        for j in range(i - 1, -1, -1):
            q = round(mu[i][j])
            if q:
                b[i] = [x - q * y for x, y in zip(b[i], b[j])]
                bs, mu = gso()
        if _dot(bs[i], bs[i]) >= (delta - mu[i][i - 1] ** 2) * _dot(bs[i - 1], bs[i - 1]):
            i += 1
        else:
            b[i], b[i - 1] = b[i - 1], b[i]
            bs, mu = gso()
            i = max(i - 1, 1)
    return b


def babai(e0, basis):
    """exact nearest-plane-style recentring: subtract the rounded exact least-squares combination of the basis."""
    k = len(basis)
    if not k:
        return e0
    G = sp.Matrix(k, k, lambda i, j: _dot(basis[i], basis[j]))
    rhs = sp.Matrix([-_dot(basis[i], e0) for i in range(k)])
    z = G.LUsolve(rhs)
    zr = [int(sp.floor(z[i] + sp.Rational(1, 2))) for i in range(k)]
    return [e0[i] + sum(zr[j] * basis[j][i] for j in range(k)) for i in range(len(e0))]


class ExactOrders:
    def __init__(self, vecs, nq, extra_denominators=()):
        self.vecs = [[Fraction(x) for x in v] for v in vecs]
        self.nq, self.dim, self.n = nq, len(self.vecs[0]), len(self.vecs)
        self.den = lcm(*[x.denominator for v in self.vecs for x in v], *extra_denominators)
        nd = self.dim - nq
        A = sp.zeros(self.dim, self.n + nd)
        for i, v in enumerate(self.vecs):
            for c in range(self.dim):
                A[c, i] = int(v[c] * self.den)
        for j in range(nd):
            A[nq + j, self.n + j] = -self.den
        D, U, V = smith_normal_decomp(A, domain=sp.ZZ)
        self.U = [[int(x) for x in U.row(i)] for i in range(U.rows)]
        self.Vm = V
        self.Vi = [[int(x) for x in V.row(i)] for i in range(V.rows)]
        self.diag = [int(D[i, i]) for i in range(min(D.shape))]
        self.r = sum(1 for d in self.diag if d != 0)
        N = V[:, self.r:]
        if N.shape[1]:
            rows = [[int(N[i, j]) for i in range(self.n)] for j in range(N.shape[1])]
            self.B = lll_reduce(rows)
            self.B = lll_reduce(self.B)
        else:
            self.B = []
        self.k = len(self.B)
        self.Nf = np.array(self.B, dtype=float).T if self.B else np.zeros((self.n, 0))
        self.cache = {}

    def solve(self, target):
        key = tuple(Fraction(x) for x in target)
        if key in self.cache:
            return self.cache[key]
        res = self._solve(key)
        self.cache[key] = res
        return res

    def _solve(self, target):
        t = []
        for x in target:
            y = x * self.den
            if y.denominator != 1:
                return None
            t.append(int(y))
        Ut = [sum(a * b for a, b in zip(row, t)) for row in self.U]
        ncol = len(self.Vi[0])
        y0 = []
        for i in range(ncol):
            if i < self.r:
                if Ut[i] % self.diag[i] != 0:
                    return None
                y0.append(Ut[i] // self.diag[i])
            else:
                y0.append(0)
        if any(Ut[i] != 0 for i in range(self.r, len(self.U))):
            return None
        e0 = [sum(a * b for a, b in zip(self.Vi[i], y0)) for i in range(self.n)]
        for _ in range(3):  # exact recentring
            e1 = babai(e0, self.B)
            if e1 == e0:
                break
            e0 = e1
        if not self.k:
            e = e0
        else:
            assert max(abs(x) for x in e0) < 2 ** 40
            r = milp(self.Nf.sum(axis=0), constraints=[LinearConstraint(self.Nf, -np.array(e0, dtype=float), np.inf)],
                     integrality=np.ones(self.k), bounds=Bounds(-1e6 * np.ones(self.k), 1e6 * np.ones(self.k)),
                     options={"time_limit": 300})
            if r.status == 2:
                return None
            if r.status != 0:
                raise RuntimeError(f"exact orders: milp status {r.status} ({r.message})")
            zr = [int(round(v)) for v in r.x]
            e = [e0[i] + sum(zr[j] * self.B[j][i] for j in range(self.k)) for i in range(self.n)]
        if any(x < 0 for x in e):
            return None
        tot = [sum(Fraction(e[i]) * self.vecs[i][c] for i in range(self.n)) for c in range(self.dim)]
        assert all(tot[c] == target[c] for c in range(self.nq))
        assert all((tot[c] - target[c]).denominator == 1 for c in range(self.nq, self.dim))
        return sum(e), e

    def lattice_solvable(self, target):
        """does an INTEGER (any-sign) exponent vector exist?  False = an exact discrete/U(1) charge of the vacuum forbids
        the operator at every order (protection by an unbroken symmetry); True with order()=None = holomorphy/cone
        obstruction only."""
        t_ = []
        for x in target:
            y = Fraction(x) * self.den
            if y.denominator != 1:
                return False
            t_.append(int(y))
        Ut = [sum(a * b for a, b in zip(row, t_)) for row in self.U]
        if any(Ut[i] % self.diag[i] for i in range(self.r)):
            return False
        return not any(Ut[i] for i in range(self.r, len(self.U)))

    def order(self, target):
        """lowest degree of S^e with sum e_i v_i = target (raw right-hand side), or None."""
        r = self.solve(target)
        return None if r is None else r[0]

    def coupling_lattice(self, field_vecs, wvec):
        target = [Fraction(w) - sum(Fraction(v[c]) for v in field_vecs) for c, w in enumerate(wvec)]
        return self.lattice_solvable(target)

    def coupling(self, field_vecs, wvec):
        """lowest number of vacuum fields that complete the coupling prod(fields) into a superpotential term."""
        target = [Fraction(w) - sum(Fraction(v[c]) for v in field_vecs) for c, w in enumerate(wvec)]
        return self.order(target)


def brute_force_order(vecs, nq, target, dmax):
    """independent check: smallest degree <= dmax, or None (by enumeration)."""
    import itertools
    den = lcm(*[Fraction(x).denominator for v in vecs for x in v], *[Fraction(x).denominator for x in target])
    A = [[int(Fraction(x) * den) for x in v] for v in vecs]
    t = [int(Fraction(x) * den) for x in target]
    for deg in range(dmax + 1):
        for combo in itertools.combinations_with_replacement(range(len(vecs)), deg):
            tot = [sum(A[i][c] for i in combo) for c in range(len(t))]
            if all(tot[c] == t[c] for c in range(nq)) and all((tot[c] - t[c]) % den == 0 for c in range(nq, len(t))):
                return deg
    return None
