"""Selection-rule engine for the SO(16)xSO(16) (Z2W x Z3) A8 models (Passes 11097, 11098).

Fields from our patched orbifolder (levdump2): sector (k = Witten index, l = theta power), lattice vector n (fixed-point
labels), U(1) charges (standard basis), non-Abelian dims, right-mover q_sh, left oscillator contribution.
Hypercharge, colour and SU(2)_L from the non-SUSY orbifolder's SM configuration (applied to our basis by exact linear
algebra, as in Pass 11095).

Rules for a coupling (vertex operators of the listed fields; conjugates are separate listed fields):
  gauge   : sum of U(1) charges = 0; non-Abelian invariance (SM exactly; hidden: see VARIANT)
  Z2W     : sum k = 0 mod 2           (Witten element acts trivially on space)
  PG      : sum l = 0 mod 3
  SG      : sum of fixed-point classes (n_{2a-1} + n_{2a}) = 0 mod 3, a = 1,2,3   (Lambda / (1-theta) Lambda = Z3^3)
  R       : sum of R^a = q_sh^a + osc^a (a = 1,2,3) = 0 exactly at cubic order (tree-level H-momentum conservation,
            pictures (-1/2,-1/2,-1)), = 0 mod 3 at higher order (picture changing, prime plane)
A fermion bilinear psi_i psi_j times any monomial of SM-neutral scalars with VEVs is allowed iff the total charge
vector lies in the Z-span of the VEV scalars' vectors plus the periodicities (non-SUSY: no holomorphy, phi and phi* are
both listed).  This is an existence test at some order, not a coupling strength."""
import json
import re
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from math import gcd

import networkx as nx
import sympy as sp

sys.path.insert(0, sys.path[0])
from tach_charges import R as Rsym, dimof  # noqa: E402
from verify_sm import SMY, span_equal  # noqa: E402


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


def parse_theirs(path):
    from tach_charges import parse_blocks
    return parse_blocks(path)


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
