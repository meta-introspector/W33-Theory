#!/usr/bin/env python3
"""Fast exact realizability of a matter parity: an independent engine for `realizable(vecs, odd)` of Pass 10967.

Question.  Charge vectors v = (q_v | e_v): q_v in Q^nq on the U(1) directions, e_v on the discrete directions
(space group, R), each discrete coordinate j of order o_j.  Is there a group element x = (x_c | x_d), x_c real and
x_d integer (the discrete periodicities force x_d in Z^nd, and only x_d mod o matters), with
    v.x in 1/2 + Z for every odd v,   v.x in Z for every even v ?

Method.  For fixed x_d a real x_c exists iff  sum_v c_v (tau_v - e_v.x_d) in Z  for every integer relation c among the
U(1) parts (sum_v c_v q_v = 0), tau = 1/2 (odd) or 0 (even).  The saturated relation lattice depends only on the
U(1) parts, and x_d runs over a finite grid, so the test is a vectorised check over the grid.  Exact: all arithmetic
is in integers scaled by the common denominator.
"""
from __future__ import annotations

import itertools
from fractions import Fraction
from math import lcm

import numpy as np
import sympy as sp


def _relations(Q):
    """saturated integer basis of {c : c^T Q = 0} for a rational matrix Q (rows = vectors)."""
    if not Q:
        return []
    M = sp.Matrix(Q).T  # nq x V
    ker = M.nullspace()
    if not ker:
        return []
    B = sp.Matrix.hstack(*[v * sp.ilcm(*[sp.fraction(x)[1] for x in v]) for v in ker]).T  # k x V rational->int
    B = sp.Matrix([[int(x) for x in B.row(i)] for i in range(B.rows)])
    # saturate: the lattice ker ∩ Z^V is spanned by the rows of U^{-1}[:k] from the Smith form of B
    from sympy.matrices.normalforms import smith_normal_decomp
    D, U, V = smith_normal_decomp(B, domain=sp.ZZ)
    Vi = V.inv()
    k = sum(1 for i in range(min(D.shape)) if D[i, i] != 0)
    return [[int(x) for x in Vi.row(i)] for i in range(k)]


class FastParity:
    def __init__(self, nq, orders, odd_vecs, even_vecs):
        self.nq, self.orders = nq, list(orders)
        self.odd = [[Fraction(x) for x in v] for v in odd_vecs]
        self.even = [[Fraction(x) for x in v] for v in even_vecs]
        allv = self.odd + self.even
        self.den = lcm(*[x.denominator for v in allv for x in v], 2) if allv else 2
        nd = len(self.orders)
        self.grid = np.array(list(itertools.product(*[range(o) for o in self.orders])), dtype=np.int64).reshape(-1, nd)
        self._rel_cache = {}

    def _phase_int(self, v):
        """den * e_v . x_d over the grid (as integers)"""
        e = np.array([int(Fraction(x) * self.den) for x in v[self.nq:]], dtype=np.int64)
        return self.grid @ e

    def realizable(self, extra_even=()):
        vecs = self.odd + self.even + [[Fraction(x) for x in v] for v in extra_even]
        tau = [Fraction(1, 2)] * len(self.odd) + [Fraction(0)] * (len(vecs) - len(self.odd))
        key = tuple(tuple(v[:self.nq]) for v in vecs)
        rel = self._rel_cache.get(key)
        if rel is None:
            rel = _relations([[sp.Rational(x.numerator, x.denominator) for x in v[:self.nq]] for v in vecs])
            self._rel_cache[key] = rel
        ok = np.ones(len(self.grid), dtype=bool)
        ph = np.stack([self._phase_int(v) for v in vecs], axis=1)  # grid x V   (den * e.x_d)
        t = np.array([int(x * self.den) for x in tau], dtype=np.int64)
        for c in rel:
            cv = np.array(c, dtype=np.int64)
            val = (t @ cv) - ph @ cv          # den * sum c (tau - e.x_d)
            ok &= (val % self.den) == 0
            if not ok.any():
                return False
        return bool(ok.any())


def realizable_choice(nq, orders, odd_vecs, even_vecs, types):
    """Exact: is there one variant per type (types = list of lists of charge vectors; all variants of a type share
    their U(1) part) such that odd_vecs are odd and even_vecs + the chosen variants are even?  Returns the list of
    chosen variant indices, or None.

    For fixed discrete element x_d the condition is: for every relation c among the U(1) parts,
        sum_base c_v (tau_v - e_v.x_d) - sum_u c_u (e_{k_u}.x_d) in Z.
    Relations not touching the types filter x_d; for the rest a dynamic programme over the finite group of residue
    vectors decides the existence of a choice (subset-sum), with backtracking for the witness."""
    odd = [[Fraction(x) for x in v] for v in odd_vecs]
    even = [[Fraction(x) for x in v] for v in even_vecs]
    typs = [[[Fraction(x) for x in v] for v in vs] for vs in types]
    for vs in typs:
        assert all(v[:nq] == vs[0][:nq] for v in vs)
    reps = [vs[0] for vs in typs]
    allv = odd + even + reps
    den = lcm(*[x.denominator for v in odd + even + [w for vs in typs for w in vs] for x in v], 2)
    rel = _relations([[sp.Rational(x.numerator, x.denominator) for x in v[:nq]] for v in allv])
    nb = len(odd) + len(even)
    grid = np.array(list(itertools.product(*[range(o) for o in orders])), dtype=np.int64).reshape(-1, len(orders))
    ph = lambda v: grid @ np.array([int(Fraction(x) * den) for x in v[nq:]], dtype=np.int64)
    tau = np.array([den // 2] * len(odd) + [0] * len(even), dtype=np.int64)
    base_ph = np.stack([ph(v) for v in odd + even], axis=1) if nb else np.zeros((len(grid), 0), dtype=np.int64)
    R = np.array(rel, dtype=np.int64).reshape(-1, len(allv))
    Rb, Ru = R[:, :nb], R[:, nb:]
    touching = np.any(Ru != 0, axis=1)
    bval = (tau @ Rb.T)[None, :] - base_ph @ Rb.T            # grid x g   (den * base part)
    ok = np.all((bval[:, ~touching] % den) == 0, axis=1)
    idx = np.nonzero(ok)[0]
    Rt = Ru[touching]                                        # g' x m
    if Rt.shape[0] == 0:
        return [0] * len(typs) if len(idx) else None
    var_ph = [np.stack([ph(v) for v in vs], axis=1) for vs in typs]   # per type: grid x variants
    for xi in idx:
        target = tuple(int(x) % den for x in bval[xi, touching])       # need sum_u c_u e_k.x_d == base part
        reach = {tuple([0] * Rt.shape[0]): None}
        layers = []
        for u in range(len(typs)):
            contrib = [tuple(int(Rt[g, u] * var_ph[u][xi, k]) % den for g in range(Rt.shape[0])) for k in range(len(typs[u]))]
            nxt = {}
            for s in reach:
                for k, c in enumerate(contrib):
                    t = tuple((a + b) % den for a, b in zip(s, c))
                    if t not in nxt:
                        nxt[t] = (s, k)
            layers.append(nxt)
            reach = nxt
        if target in reach:
            choice, s = [], target
            for u in range(len(typs) - 1, -1, -1):
                s, k = layers[u][s]
                choice.append(k)
            return choice[::-1]
    return None
