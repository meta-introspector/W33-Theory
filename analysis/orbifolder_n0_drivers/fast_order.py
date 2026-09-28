"""Minimal order of a monomial of listed fields with prescribed charge (continuous U(1)s exact, discrete coordinates mod
their periods) by a direct integer program (HiGHS).  Validated against w33_exact_monomial_orders.ExactOrders."""
from fractions import Fraction as F
from math import lcm

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp


class FastOrders:
    def __init__(self, vecs, nq, periods, emax=60):
        """vecs: charge vectors (continuous first nq coords, then discrete coords with the given periods)"""
        keys = {}
        for v in vecs:
            keys.setdefault(tuple(F(x) for x in v), 0)
        self.vecs = [list(k) for k in keys]            # identical charge vectors merged (orders unchanged)
        self.nq, self.periods, self.emax = nq, list(periods), emax
        self.dim = nq + len(periods)
        self.den = lcm(*[x.denominator for v in self.vecs for x in v], *[F(p).denominator for p in periods], 1)
        n, nd = len(self.vecs), len(periods)
        A = np.zeros((self.dim, n + nd))
        for i, v in enumerate(self.vecs):
            for c in range(self.dim):
                A[c, i] = int(v[c] * self.den)
        for j, p in enumerate(periods):
            A[nq + j, n + j] = int(F(p) * self.den)
        self.A, self.n, self.nd = A, n, nd
        self.cobj = np.concatenate([np.ones(n), np.zeros(nd)])
        self.bounds = Bounds(np.concatenate([np.zeros(n), -1e4 * np.ones(nd)]),
                             np.concatenate([emax * np.ones(n), 1e4 * np.ones(nd)]))

    def order(self, target, min_total=0):
        t = [F(x) * self.den for x in target]
        if any(x.denominator != 1 for x in t):
            return None
        b = np.array([float(int(x)) for x in t])
        cons = [LinearConstraint(self.A, b, b)]
        if min_total:
            cons.append(LinearConstraint(self.cobj.reshape(1, -1), min_total, np.inf))
        r = milp(self.cobj, constraints=cons, integrality=np.ones(self.n + self.nd),
                 bounds=self.bounds, options={"time_limit": 30})
        if r.status == 2:
            return None
        if r.status != 0:
            return 'timeout'
        e = [int(round(x)) for x in r.x[:self.n]]
        m = [int(round(x)) for x in r.x[self.n:]]
        tot = [sum(F(e[i]) * self.vecs[i][c] for i in range(self.n)) for c in range(self.dim)]
        assert all(tot[c] == F(target[c]) for c in range(self.nq))
        assert all(tot[self.nq + j] + m[j] * F(self.periods[j]) == F(target[self.nq + j]) for j in range(self.nd))
        return sum(e)

    def feasible(self, target, max_total):
        """is there a monomial of total degree <= max_total with the target charge?  True / False / 'timeout'"""
        t = [F(x) * self.den for x in target]
        if any(x.denominator != 1 for x in t):
            return False
        b = np.array([float(int(x)) for x in t])
        cons = [LinearConstraint(self.A, b, b), LinearConstraint(self.cobj.reshape(1, -1), -np.inf, max_total)]
        ub = np.concatenate([min(self.emax, max_total) * np.ones(self.n), 1e4 * np.ones(self.nd)])
        r = milp(self.cobj, constraints=cons, integrality=np.ones(self.n + self.nd),
                 bounds=Bounds(np.concatenate([np.zeros(self.n), -1e4 * np.ones(self.nd)]), ub), options={"time_limit": 20})
        if r.status == 2:
            return False
        if r.status != 0:
            return 'timeout'
        return True

    def coupling(self, field_vecs, wvec, min_total=0):
        target = [F(w) - sum(F(v[c]) for v in field_vecs) for c, w in enumerate(wvec)]
        return self.order(target, min_total)
