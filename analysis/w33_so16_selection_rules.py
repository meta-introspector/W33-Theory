#!/usr/bin/env python3
"""Coupling selection rules for the SO(16)xSO(16) (Z2W x Z3) W(3,3) A8 models of Pass 11095 (shared by Passes 11097, 11098).

Fields come from the patched orbifolder (analysis/orbifolder_n0_drivers/levdump2.cpp): sector (k = Witten index,
l = theta power), lattice vector n, U(1) charges, non-Abelian dims, right-mover q_sh and left oscillator contribution.
The Standard-Model embedding (hypercharge, colour, SU(2)_L) comes from the non-SUSY orbifolder's SM configuration and is
carried into our basis by exact linear algebra (factors matched by root span).

Rules for a coupling of listed vertex operators (conjugate states are separate listed fields):
  gauge : U(1) charges sum to 0; non-Abelian invariance (SM exactly; hidden: conjugate pair or singlets)
  Z2W   : sum k = 0 mod 2 (the Witten element acts trivially on space)
  PG    : sum l = 0 mod 3
  SG    : sum of fixed-point classes (n_{2a-1} + n_{2a}) mod 3 = 0 in each torus (Lambda/(1-theta)Lambda = Z3^3)
  R     : H-momentum R^a = q_sh^a + osc^a: sum EXACTLY 0 for a cubic psi psi phi (pictures -1/2,-1/2,-1), 0 mod 3 at
          higher order (picture changing on a prime plane)
The cubic rule was checked on the untwisted U1 U2 U3 coupling of a SUSY orbifold (internal q_sh sum 0).
"""
from __future__ import annotations


import json
import re
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from math import gcd

import networkx as nx
import sympy as sp



def ex(x):
    f = F(float(x)).limit_denominator(10 ** 6)
    assert abs(float(f) - float(x)) < 1e-9, x
    assert 5184 % f.denominator == 0, (x, f)
    return f


def parse_ours(path):
    blocks, cur = {}, None
    for ln in open(path, errors='ignore'):
        m = re.match(r'MODEL (\d+) (\S+)', ln)
        if m:
            cur = int(m.group(1))
            blocks[cur] = dict(label=m.group(2), factors=[], u1=[], fields=[])
            continue
        if cur is None:
            continue
        b = blocks[cur]
        m = re.match(r'FACTOR (\d+) (\S+) \|(.*)', ln)
        if m:
            b['factors'].append((m.group(2), [[ex(x) for x in r.split()] for r in m.group(3).split('|')]))
            continue
        m = re.match(r'U1STD (\d+) \|(.*)', ln)
        if m:
            b['u1'].append([ex(x) for x in m.group(2).split()])
            continue
        m = re.match(r'S k=(\d+) l=(\d+) m=(\d+) dim=(\S+) q=(\S+) osc=(\S+) n=(\S+) qsh=(\S+)', ln)
        if m:
            b['fields'].append(dict(k=int(m.group(1)), l=int(m.group(2)), m=int(m.group(3)), dim=m.group(4).split(','),
                                    q=[ex(x) for x in m.group(5).split(',')],
                                    osc=[ex(x) for x in m.group(6).split(',')],
                                    n=[F(x) for x in m.group(7).split(',')],
                                    qsh=[ex(x) for x in m.group(8).split(',')]))
    return blocks




def sm_data(th, us):
    """hypercharge coefficients in our U(1) basis, colour and SU(2) factor indices in our factor list"""
    lab = [f for f in th['fields'] if f['lab'] in SMY and f['lab'] not in ('n', 'bn')]
    c3 = next(j for j, x in enumerate(next(f for f in lab if f['lab'] == 'q')['dim']) if dimof(x) == 3)
    w2 = next(j for j, x in enumerate(next(f for f in lab if f['lab'] == 'l')['dim']) if dimof(x) == 2)
    qcol = next(f for f in lab if f['lab'] == 'q')['dim'][c3]
    A = sp.Matrix([[Rsym(x) for x in f['q']] for f in lab])
    bY = sp.Matrix([Rsym(SMY[f['lab']]) for f in lab])
    y, params = A.gauss_jordan_solve(bY)
    y = y.subs({p: 0 for p in params})
    assert A * y == bY
    tY = sp.Matrix([[Rsym(x) for x in u] for u in th['u1']]).T * y
    oc = [j for j, (a, roots) in enumerate(us['factors']) if span_equal(roots, th['factors'][c3][1])]
    ow = [j for j, (a, roots) in enumerate(us['factors']) if span_equal(roots, th['factors'][w2][1])]
    assert len(oc) == 1 and len(ow) == 1
    TST = sp.Matrix([[Rsym(x) for x in u] for u in us['u1']])
    c, params = TST.T.gauss_jordan_solve(tY)
    c = c.subs({p: 0 for p in params})
    assert TST.T * c == tY
    # which colour token is the quark triplet in OUR convention: match by the q field's colour token in THEIR config and
    # the span match (same roots, same weights) -> tokens agree
    return [F(int(sp.fraction(x)[0]), int(sp.fraction(x)[1])) for x in c], oc[0], ow[0], qcol


def classes(n):
    return [int((n[2 * a] + n[2 * a + 1]) % 3) for a in range(3)]


def charge_vector(f, rule_R=True):
    """[U(1)s..., k (mod 2), l (mod 3), SG classes (mod 3) x3, R^1..R^3 (mod 3)]"""
    v = list(f['q']) + [F(f['k']), F(f['l'])] + [F(x) for x in classes(f['n'])]
    if rule_R:
        v += [f['qsh'][a] + f['osc'][a - 1] for a in (1, 2, 3)]
    return v


def periods(nu1, rule_R=True, exact_R=False):
    """periodicity vectors of the discrete coordinates"""
    dim = nu1 + 5 + (3 if rule_R else 0)
    P = []
    for i, p in [(nu1, 2), (nu1 + 1, 3), (nu1 + 2, 3), (nu1 + 3, 3), (nu1 + 4, 3)]:
        e = [F(0)] * dim
        e[i] = F(p)
        P.append(e)
    if rule_R and not exact_R:
        for a in range(3):
            e = [F(0)] * dim
            e[nu1 + 5 + a] = F(3)
            P.append(e)
    return P


class Lattice:
    """integer row-echelon (Hermite) basis of the Z-span of rational vectors; membership test"""

    def __init__(self, vecs):
        self.den = 1
        for v in vecs:
            for x in v:
                self.den = self.den * x.denominator // gcd(self.den, x.denominator)
        rows = [[int(x * self.den) for x in v] for v in vecs]
        self.dim = len(vecs[0]) if vecs else 0
        basis = []
        col = 0
        rows = [r for r in rows if any(r)]
        while rows and col < self.dim:
            nz = [r for r in rows if r[col] != 0]
            zr = [r for r in rows if r[col] == 0]
            while len(nz) > 1:
                nz.sort(key=lambda r: abs(r[col]))
                p = nz[0]
                new = [p]
                for r in nz[1:]:
                    q = r[col] // p[col]
                    rr = [a - q * b for a, b in zip(r, p)]
                    if rr[col] != 0:
                        new.append(rr)
                    elif any(rr):
                        zr.append(rr)
                nz = new
            if nz:
                p = nz[0]
                if p[col] < 0:
                    p = [-a for a in p]
                basis.append((col, p))
            rows = zr
            col += 1
        self.basis = basis

    def contains(self, v):
        t = [x * self.den for x in v]
        if any(x.denominator != 1 for x in t):
            return False
        t = [int(x) for x in t]
        for col, p in self.basis:
            if t[col] % p[col] != 0:
                return False
            q = t[col] // p[col]
            t = [a - q * b for a, b in zip(t, p)]
        return not any(t)


def sm_info(f, c, oc, ow, qcol):
    Y = sum(q * ci for q, ci in zip(f['q'], c))
    col = f['dim'][oc]
    w = dimof(f['dim'][ow])
    t = 1 if col == qcol else (-1 if dimof(col) == 3 else 0)
    charges = [Y + F(tt, 2) for tt in range(-(w - 1), w, 2)]
    frac = any(((3 * Q + t) % 3 != 0) or (3 * Q).denominator != 1 for Q in charges) if dimof(col) in (1, 3) else False
    return Y, col, w, t, frac


def conj_token(tok):
    if tok.endswith('adj') or dimof(tok) == 1:
        return tok
    return tok[1:] if tok.startswith('-') else '-' + tok

SMY = {"q": F(1, 6), "bu": F(-2, 3), "bd": F(1, 3), "l": F(-1, 2), "be": F(1), "bl": F(1, 2), "d": F(-1, 3),
       "u": F(2, 3), "bq": F(-1, 6), "e": F(-1), "bn": F(0), "n": F(0)}


def Rsym(x):
    return sp.Rational(x.numerator, x.denominator)


def dimof(t):
    return int(re.match(r'-?(\d+)', t).group(1))


def span_equal(A, B):
    MA, MB = sp.Matrix([[Rsym(x) for x in r] for r in A]), sp.Matrix([[Rsym(x) for x in r] for r in B])
    return MA.rank() == MB.rank() == sp.Matrix.vstack(MA, MB).rank()


def parse_theirs(path):
    """the non-SUSY orbifolder SM configuration dump (nsosm.cpp): FACTOR / U1 / L (labelled left fermions)"""
    blocks, cur = {}, None
    for ln in open(path, errors='ignore'):
        m = re.match(r'MODEL (\d+) (\S+)', ln)
        if m:
            cur = int(m.group(1))
            blocks[cur] = dict(label=m.group(2), factors=[], u1=[], fields=[])
            continue
        if cur is None:
            continue
        b = blocks[cur]
        m = re.match(r'FACTOR (\d+) (\S+) \|(.*)', ln)
        if m:
            b['factors'].append((m.group(2), [[ex(x) for x in r.split()] for r in m.group(3).split('|')]))
            continue
        m = re.match(r'U1 (\d+) \|(.*)', ln)
        if m:
            b['u1'].append([ex(x) for x in m.group(2).split()])
            continue
        m = re.match(r'L (\S+) dim=(\S+) q=(\S+)', ln)
        if m:
            b['fields'].append(dict(lab=m.group(1), dim=m.group(2).split(','), q=[ex(x) for x in m.group(3).split(',')]))
    return blocks


# ---------------------------------------------------------------- minimal orders (direct integer program, HiGHS)
from math import lcm  # noqa: E402

import numpy as np  # noqa: E402
from scipy.optimize import Bounds, LinearConstraint, milp  # noqa: E402



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

    def coupling(self, field_vecs, wvec, min_total=0):
        target = [F(w) - sum(F(v[c]) for v in field_vecs) for c, w in enumerate(wvec)]
        return self.order(target, min_total)
