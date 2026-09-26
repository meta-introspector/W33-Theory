#!/usr/bin/env python3
"""Pass 10978: the E6 lattice fixes the admissible R-rule of Z12-I; the last 9 models are decided; audit.

A. Lattice symmetries.  The Z12-I twist on the E6 root lattice is the Coxeter element c.  Generating W(E6)
   (51840 elements) exactly: C_W(c) = <c> (order 12) and -1 is not in W(E6).  Aut(E6) = W x {+-1}, so the
   symmetries commuting with the twist are <c> x <-1>.  -c^6 acts as -1 exactly on the order-3 plane
   (eigenphases 4, 8 of c) and trivially elsewhere.  Hence the ONLY R-symmetry beyond the orbifold projection
   is a Z2^R on the order-3 plane: sum_alpha R^3_alpha = -1 mod 2.  The per-plane Z12 x Z12 x Z3 rules used
   tentatively in Pass 10968 are not supported by the lattice.
B. Bracketing.  For the 14 Z12-I models that reach a parity-viable FI-cancelling vacuum, every verdict is
   computed with the weakest rules (gauge + point group) and with the strongest admissible ones (+ Z2^R).
   Exotic d-triplet rank: all-orders lattice rank (an upper bound: 'deficient' is robust).  Higgs rank:
   explicit monomials to order 9 with random O(1) coefficients (a lower bound: 'no light pair' is robust).
C. Z12I_1063_303: with Z2^R its 5-direction vacuum has exotics full rank and exactly one Higgs pair
   protected at ALL orders, but an outside singlet n_36 appears linearly at order 4 (F ~ eps^3 M_s^2), and
   absorbing any low-order linear singlet lifts the protected pair.  Flat-or-light, again.
D. Audit of Pass 10967's Z6II_23 statement 'no u^c d^c d^c at orders 3, 4': it holds with orbifolder rules,
   with the established prime-plane rules only, and with the gamma-corrected G2 charge.
"""
from __future__ import annotations

import gzip
import importlib.util
import itertools
import json
import re
from collections import Counter, deque
from fractions import Fraction
from math import lcm
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form
from scipy.optimize import linear_sum_assignment, milp, LinearConstraint, Bounds

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("p10967", ROOT / "analysis" / "w33_pass10967_fi_vacuum_regenerates_rpv.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
Z12 = ROOT / "data" / "w33_pass10968_z12_left_chiral_ledger.json.gz"
P10968 = ROOT / "data" / "w33_pass10968_r_symmetry_parity_and_massless_exotics.json"
GAMMA = ROOT / "data" / "w33_pass10974_z6ii23_gamma_phases.json.gz"
OUT = ROOT / "data" / "w33_pass10978_z12_e6_admissible_r_rule.json"
ab = lambda s: abs(int(re.match(r"(-?\d+)", s).group(1)))


# ------------------------------------------------------------------ A
def e6_centralizer():
    Cm = np.array([[2, 0, -1, 0, 0, 0], [0, 2, 0, -1, 0, 0], [-1, 0, 2, -1, 0, 0],
                   [0, -1, -1, 2, -1, 0], [0, 0, 0, -1, 2, -1], [0, 0, 0, 0, -1, 2]])
    S = []
    for i in range(6):
        Mx = np.eye(6, dtype=int)
        for j in range(6):
            Mx[i, j] -= Cm[j, i]
        S.append(Mx)
    I = np.eye(6, dtype=int)
    seen = {I.tobytes(): I}
    dq = deque([I])
    while dq:
        g = dq.popleft()
        for s in S:
            h = s @ g
            if h.tobytes() not in seen:
                seen[h.tobytes()] = h
                dq.append(h)
    W = list(seen.values())
    c = I.copy()
    for s in S:
        c = c @ s
    order = next(k for k in range(1, 50) if np.array_equal(np.linalg.matrix_power(c, k), I))
    cent = [g for g in W if np.array_equal(g @ c, c @ g)]
    powers = [np.linalg.matrix_power(c, k) for k in range(order)]
    ev = sorted(np.round(np.linalg.eigvals((-np.linalg.matrix_power(c, 6)).astype(float)).real, 9).tolist())
    phases = sorted(int(round(np.angle(z) / (2 * np.pi) * 12)) % 12 for z in np.linalg.eigvals(c.astype(float)))
    return dict(W_order=len(W), coxeter_order=order, centralizer_order=len(cent),
                centralizer_is_cyclic_c=all(any(np.array_equal(g, p) for p in powers) for g in cent),
                minus_one_in_W=any(np.array_equal(g, -I) for g in W), eigenvalues_minus_c6=ev, eigenphases_c=phases)


# ------------------------------------------------------------------ helpers for Z12-I
def compositions(n, k):
    for tot in range(1, k + 1):
        for c in itertools.combinations(range(tot + n - 1), n - 1):
            prev, out = -1, []
            for x in c:
                out.append(x - prev - 1)
                prev = x
            out.append(tot + n - 2 - prev)
            yield out


class Z12Model:
    def __init__(self, m, with_z2):
        fl = m["left"]
        self.fl = fl
        self.nq = len(fl[0]["q"])
        order = m["pg_order"]
        self.V = {f["name"]: [Fraction(x) for x in f["q"]] + [Fraction(f["pg"]) / order] for f in fl}
        self.W = [Fraction(0)] * (self.nq + 1)
        self.integ = [[Fraction(0)] * self.nq + [Fraction(1)]]
        if with_z2:
            for f in fl:
                self.V[f["name"]] = self.V[f["name"]] + [Fraction(f["rq"].split(",")[2]) / 2]
            self.W = self.W + [Fraction(-1, 2)]
            self.integ = [v + [Fraction(0)] for v in self.integ] + [[Fraction(0)] * (self.nq + 1) + [Fraction(1)]]
        self.dims = {f["name"]: f["dim"] for f in fl}
        self.den = lcm(*[x.denominator for v in list(self.V.values()) + [self.W] for x in v])

    def lab(self, p):
        return [f["name"] for f in self.fl if M.P.base_of(f["name"]) == p]

    def sing(self):
        return [n for n in self.V if M.P.base_of(n) == "n" and all(ab(x) == 1 and "adj" not in x for x in self.dims[n].split(","))]

    def lattice(self, S):
        gens = [self.V[s] for s in S] + self.integ
        G = sp.Matrix([[int(x * self.den) for x in g] for g in gens]).T
        H = hermite_normal_form(G)
        B = sp.Matrix.hstack(*[H[:, j] for j in range(H.cols) if any(H[:, j])])

        def inl(u):
            try:
                c, prm = B.gauss_jordan_solve(sp.Matrix([int(x * self.den) for x in u]))
            except ValueError:
                return False
            c = c.subs({p: 0 for p in prm})
            return all(sp.fraction(x)[1] == 1 for x in c)
        return inl

    def allowed(self, fields):
        t = [sum(self.V[f][k] for f in fields) - self.W[k] for k in range(len(self.W))]
        return all(v == 0 for v in t[:self.nq]) and all(v.denominator == 1 for v in t[self.nq:])

    def lattice_rank(self, S, Xb, X):
        inl = self.lattice(S)
        R_, C_ = self.lab(Xb), self.lab(X)
        if not R_:
            return (0, 0)
        B = np.array([[int(inl([p + q - w for p, q, w in zip(self.V[a], self.V[b], self.W)])) for b in C_] for a in R_])
        rr, cc = linear_sum_assignment(-B)
        return int(B[rr, cc].sum()), len(R_)

    def numeric_rank(self, S, Xb, X, K=9):
        rng = np.random.default_rng(7)
        vev = {s: 0.3 * np.exp(2j * np.pi * rng.random()) * (0.5 + rng.random()) for s in S}
        monos = [[0] * len(S)] + list(compositions(len(S), K - 2))
        R_, C_ = self.lab(Xb), self.lab(X)
        Mx = np.zeros((len(R_), len(C_)), complex)
        for i, a in enumerate(R_):
            for j, b in enumerate(C_):
                for mo in monos:
                    flds = [a, b] + [s for s, e in zip(S, mo) for _ in range(e)]
                    if self.allowed(flds):
                        Mx[i, j] += (rng.normal() + 1j * rng.normal()) * np.prod([vev[s] ** e for s, e in zip(S, mo)])
        sv = np.linalg.svd(Mx, compute_uv=False)
        return int((sv > 1e-12 * max(sv.max(), 1e-300)).sum()), len(R_)

    def ilp_rank(self, S, Xb, X):
        """structural rank over entries realised by an explicit monomial with NON-NEGATIVE exponents (<= 200 each)."""
        R_, C_ = self.lab(Xb), self.lab(X)
        if not R_:
            return (0, 0)
        B = np.array([[int(self.min_order(S, [p + q - w for p, q, w in zip(self.V[a], self.V[b], self.W)]) is not None)
                       for b in C_] for a in R_])
        rr, cc = linear_sum_assignment(-B)
        return int(B[rr, cc].sum()), len(R_)

    def min_order(self, S, target, cap=200):
        n, dim = len(S), len(self.W)
        A = np.array([[int(x * self.den) for x in self.V[s]] for s in S], dtype=float).T
        ncong = dim - self.nq
        t = np.array([int(x * self.den) for x in target], dtype=float)
        Aeq = np.zeros((dim, n + ncong)); Aeq[:, :n] = A
        for j in range(ncong):
            Aeq[self.nq + j, n + j] = -self.den
        r = milp(np.concatenate([np.ones(n), np.zeros(ncong)]), constraints=[LinearConstraint(Aeq, -t, -t)],
                 integrality=np.ones(n + ncong), bounds=Bounds(np.concatenate([np.zeros(n), -1e7 * np.ones(ncong)]),
                                                             np.concatenate([cap * np.ones(n), 1e7 * np.ones(ncong)])),
                 options={"time_limit": 60})
        return int(round(r.fun)) if r.status == 0 else None


def verdict(dW, hW, dZ, hZ, hLW=None, iW=None, hIW=None):
    """dW/dZ: exotic-d rank (weak / Z2); hW/hZ: Higgs rank, explicit monomials to order 9 (lower bounds);
    hLW: Higgs lattice rank (weak); iW / hIW: exotic-d / Higgs rank from non-negative monomials at any order (weak)."""
    if iW is not None and iW[0] < iW[1]:
        return "dead: exotic d massless at all orders even with the weakest rules"
    if hIW is not None and hIW[0] < hIW[1] - 1:
        return "dead: extra Higgs doublets massless at all orders even with the weakest rules"
    exotic_dead_weak = dW[0] < dW[1]
    exotic_dead_z2 = dZ[0] < dZ[1]
    need = hW[1] - 1
    higgs_heavy_weak, higgs_heavy_z2 = hW[0] > need, hZ[0] > need
    if exotic_dead_weak:
        return "dead: exotic d massless at all orders even with the weakest rules"
    if hLW is not None and hLW[0] < need:
        return "dead: extra Higgs doublets massless at all orders even with the weakest rules"
    if higgs_heavy_z2:
        return "dead: no light Higgs pair even with the strongest admissible rules (mu-problem)"
    if exotic_dead_z2 and higgs_heavy_weak:
        return "dead at both ends: weak rules -> no light Higgs; Z2^R -> massless exotic d"
    if not exotic_dead_z2 and hZ[0] == need:
        return "candidate under Z2^R: exotics massive, exactly one light Higgs pair"
    return "other"


def main():
    A = e6_centralizer()
    print("A", A, flush=True)
    with gzip.open(Z12, "rt") as fh:
        z12 = json.load(fh)
    prev = json.loads(P10968.read_text())["z12_models"]
    B, Bsum = {}, Counter()
    for name, r in sorted(prev.items()):
        if "support" not in r:
            continue
        S = r["support"]
        weak, z2 = Z12Model(z12[name], False), Z12Model(z12[name], True)
        dW, dZ = weak.lattice_rank(S, "d", "bd"), z2.lattice_rank(S, "d", "bd")
        hW, hZ = weak.numeric_rank(S, "bl", "l"), z2.numeric_rank(S, "bl", "l")
        hLW, hLZ = weak.lattice_rank(S, "bl", "l"), z2.lattice_rank(S, "bl", "l")
        iW, iZ = weak.ilp_rank(S, "d", "bd"), z2.ilp_rank(S, "d", "bd")
        hIW, hIZ = weak.ilp_rank(S, "bl", "l"), z2.ilp_rank(S, "bl", "l")
        v = verdict(dW, hW, iZ, hZ, hLW, iW, hIW)
        B[name] = dict(support=S, exotic_d_lattice=dict(weak=dW, z2=dZ), exotic_d_ilp=dict(weak=iW, z2=iZ),
                       higgs_numeric=dict(weak=hW, z2=hZ), higgs_lattice=dict(weak=hLW, z2=hLZ),
                       higgs_ilp=dict(weak=hIW, z2=hIZ), verdict=v)
        Bsum[v.split(":")[0]] += 1
        print(name.split("|")[1], "exotic ilp", iW, iZ, " higgs ilp", hIW, hIZ, " higgs num", hW, hZ, "->", v, flush=True)
    # ---- C: the candidate
    name = "Z12-I|Z12I_1063__SM_20260926_303"
    z2 = Z12Model(z12[name], True)
    S = B[name]["support"]
    inl = z2.lattice(S)
    lin = []
    for o in z2.sing():
        if o in S or not inl([p - w for p, w in zip(z2.V[o], z2.W)]):
            continue
        k = z2.min_order(S, [p - w for p, w in zip(z2.V[o], z2.W)])
        if k is not None:
            lin.append((k + 1, o))
    lin.sort()
    absorb = {}
    for k, o in lin[:5]:
        S2 = S + [o]
        absorb[o] = dict(order=k, higgs_lattice=z2.lattice_rank(S2, "bl", "l"), exotic_d_lattice=z2.lattice_rank(S2, "d", "bd"))
    Cres = dict(higgs_lattice=z2.lattice_rank(S, "bl", "l"), exotic_d_lattice=z2.lattice_rank(S, "d", "bd"),
                W_S_possible=inl([-w for w in z2.W]), linear_terms=lin, single_absorptions=absorb)
    print("C", Cres, flush=True)
    # ---- D: audit of Pass 10967 udd statement in Z6II_23
    ledger, sha = M.P.load_ledger()
    disc = M.load(M.DISC)
    zname = "Z6-II|Z6II_23__SM_20260917_2913"
    d = disc[zname.split("|")[1]]
    with gzip.open(GAMMA, "rt") as fh:
        G = json.load(fh)
    Q = {f["name"]: [Fraction(x) for x in f["q"]] for f in ledger[zname]["left"]}
    NR = {f: [Fraction(x) for x in v["nonR"]] for f, v in d["fields"].items()}
    RR = {f: [Fraction(x) for x in v["R"]] for f, v in d["fields"].items()}
    gam = {f: Fraction(G[f]["gamma"].split(",")[0]) for f in G}
    RN = [(int(x.split(":")[0]), Fraction(x.split(":")[1])) for x in d["R"]]
    nq = len(next(iter(Q.values())))

    def ok(fields, mode):
        if any(sum(Q[f][i] for f in fields) != 0 for i in range(nq)):
            return False
        for j, o in enumerate(d["nonR_orders"]):
            if sum(NR[f][j] for f in fields) % o != 0:
                return False
        for j in {"orbifolder": [0, 1, 2], "prime_planes": [1, 2], "gamma_corrected": [0, 1, 2]}[mode]:
            n, w = RN[j]
            tot = sum(RR[f][j] + (6 * gam[f] if (mode == "gamma_corrected" and j == 0) else 0) for f in fields)
            if (tot - w) % n != 0:
                return False
        return True
    lab = lambda p: sorted(f for f in Q if M.P.base_of(f) == p)
    Dres = {}
    for mode in ("orbifolder", "prime_planes", "gamma_corrected"):
        c3 = sum(ok((u, a, b), mode) for u in lab("bu") for a, b in itertools.combinations_with_replacement(lab("bd"), 2))
        c4 = sum(ok((u, a, b, s), mode) for u in lab("bu") for a, b in itertools.combinations_with_replacement(lab("bd"), 2) for s in lab("n"))
        Dres[mode] = dict(order3=c3, order4=c4)
    print("D", Dres, flush=True)
    OUT.write_text(json.dumps(dict(pass_id=10978, A=A, B=dict(summary=dict(Bsum), models=B), C=Cres, D=Dres),
                              indent=1, sort_keys=True, default=str))


if __name__ == "__main__":
    main()
