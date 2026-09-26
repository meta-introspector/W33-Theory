#!/usr/bin/env python3
"""Pass 10979: the W(3,3) twist on the prime orbifolds T6/Z3 and T6/Z3xZ3.

Every plane of T6/Z3 and T6/Z3xZ3 on SU(3)^3 is prime (order 3), so by Bizet et al. (arXiv:1301.2322) every
R-symmetry orbifolder reports is established -- these are the first families of the programme where the
matter-parity question is decided with NO rule in doubt.  The W(3,3) twist is the A8 Kac class (the 27 fixed
points of T6/Z3 are the 27 far points of a W(3,3) point).

A. Scans (frozen): T6/Z3 with the unique A8 modular-invariant Kac pair: 0 Standard Models in 64,000 Wilson-line
   draws.  T6/Z3xZ3 (V1 = A8 pair): 13 SMs in 40,000 draws, 12 distinct spectra.
B. Census, torus x space group x all (established) R-combinations: matter parity exists in 4 of 13; one closes by
   Farkas on the even-able singlets; 3 have parity-preserving FI-cancelling D-flat rays.
C. MSSM-viable parity (exact matcher) on every realizable FI support: c1 none; c3 24 and c4 12 supports viable.
   On every one the single exotic d-triplet pair is massless to all orders and W|_S = 0 identically (the
   flat-or-massive dichotomy again), but outside singlets enter W linearly at order 3-4.
D. Absorbing those singlets (the vacuum they force): 28 distinct extended vacua with the joint parity still
   MSSM-viable, exotic d massive, FI-cancelling D-flat with every field nonzero.  The Higgs pair stays massless
   to all orders in 22 (mu from supersymmetry breaking) and gets mass at order 11-12 in 6.
E. F-flatness, exactly.  With every vacuum field nonzero, F_i = phi_i^-1 sum_k c_k e_ki m_k(phi), so F = 0
   forces sum_k (c_k m_k) e_k = 0 with every c_k m_k != 0: the exponent vectors of the W-monomials on S must be
   linearly dependent.  Enumerating all W-allowed monomials to degree 14 (exact integer charge test under every
   selection rule): the vectors are independent through degree 11 in all 28 vacua and through degree 14 in 10.
   So no supersymmetric vacuum exists from the superpotential through degree 11, for ANY nonzero couplings; the
   first possible balance is cubic/quartic against degree >= 12, i.e. |<phi>|/M_s ~ (c_low/c_12)^(1/9), a VEV set
   by a coupling ratio rather than by the FI term.  Whether string couplings supply that hierarchy is a
   worldsheet-instanton question, named here as open.
"""
from __future__ import annotations

import gzip
import importlib.util
import json
from collections import Counter
from fractions import Fraction
from math import lcm
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.optimize import Bounds, LinearConstraint, linear_sum_assignment, linprog, milp

ROOT = Path(__file__).resolve().parents[1]


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / "analysis" / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


P8 = _load("p10968", "w33_pass10968_r_symmetry_parity_and_massless_exotics.py")
P74 = _load("p10974", "w33_pass10974_valid_r_rules_close_the_rotation_door.py")
M = P8.M
BLOB = ROOT / "data" / "w33_pass10979_z3xz3_ledger.json.gz"
OUT = ROOT / "data" / "w33_pass10979_z3xz3_parity_vacua_f_obstruction.json"
DMAX = 14


def load():
    with gzip.open(BLOB, "rt") as fh:
        return json.load(fh)


def census(L, D):
    out, cs = {}, Counter()
    for name in sorted(L):
        b = name.split("|")[1]
        RN = [(int(x.split(":")[0]), Fraction(x.split(":")[1])) for x in D[b]["R"]]
        assert all(n == 3 for n, _ in RN)  # every plane prime: every R-symmetry established
        Rc = {f: [Fraction(x) for x in v["R"]] for f, v in D[b]["fields"].items()}
        r = P74.summarize(P74.with_R(name, L[name], D[b], Rc, [n for n, _ in RN], [w for _, w in RN]))
        out[name] = r
        cs["models"] += 1
        cs["parity_exists"] += bool(r["parity_exists"])
        v = r.get("verdict") or ""
        cs["closed_farkas"] += v.startswith("closed")
        cs["fi_rays_realizable"] += v == "COUNTEREXAMPLE"
    return out, dict(cs)


class Tools:
    """exact selection rules for one model: minimal non-negative monomial orders (ILP), mass ranks,
    MSSM-viable parity with a given even support."""

    def __init__(self, Z):
        self.Z = Z
        self.den = lcm(*[x.denominator for v in list(Z.vec.values()) + [Z.wvec] + Z.virtual for x in v])
        self.dim = len(Z.wvec)
        ent = lambda p: [(n_, Z.vec[n_], Z.mult(n_)) for n_ in Z.lab(p)]
        self.species = [(ent("q"), ent("bq")), (ent("bu"), ent("u")), (ent("bd"), ent("d")), (ent("be"), ent("e"))]
        self.leptons = (ent("l"), ent("bl"))

    def min_order(self, S, target, cap=150):
        Z, den, n = self.Z, self.den, len(S)
        A = np.array([[int(x * den) for x in Z.vec[s]] for s in S], dtype=float).T
        ncong = self.dim - Z.nq
        t = np.array([int(x * den) for x in target], dtype=float)
        Aeq = np.zeros((self.dim, n + ncong))
        Aeq[:, :n] = A
        for j in range(ncong):
            Aeq[Z.nq + j, n + j] = -den
        r = milp(np.concatenate([np.ones(n), np.zeros(ncong)]), constraints=[LinearConstraint(Aeq, -t, -t)],
                 integrality=np.ones(n + ncong),
                 bounds=Bounds(np.concatenate([np.zeros(n), -1e7 * np.ones(ncong)]),
                               np.concatenate([cap * np.ones(n), 1e7 * np.ones(ncong)])),
                 options={"time_limit": 30})
        return int(round(r.fun)) if r.status == 0 else None

    def rank(self, S, Xb, X):
        """(structural rank, number of Xb, lowest order reaching that rank) of the Xb.X mass matrix."""
        Z = self.Z
        R_, C_ = Z.lab(Xb), Z.lab(X)
        if not R_:
            return (0, 0, None)
        O = np.full((len(R_), len(C_)), 999)
        for i, a in enumerate(R_):
            for j, c in enumerate(C_):
                o = self.min_order(S, [p + q - w for p, q, w in zip(Z.vec[a], Z.vec[c], Z.wvec)])
                if o is not None:
                    O[i, j] = o + 2
        B = (O < 999).astype(int)
        rr, cc = linear_sum_assignment(-B)
        rk = int(B[rr, cc].sum())
        reach = None
        for k in sorted(set(O.flatten())):
            if k == 999:
                break
            Bk = (O <= k).astype(int)
            r2, c2 = linear_sum_assignment(-Bk)
            if Bk[r2, c2].sum() == rk:
                reach = int(k)
                break
        return (rk, len(R_), reach)

    def viable(self, S):
        Z = self.Z
        realz = lambda ev, od: M.realizable([list(v) for v in od] + [list(v) for v in ev] + [list(v) for v in Z.virtual],
                                            list(range(len(od))))
        return P8.match_search(realz, [Z.vec[s] for s in S], self.species, self.leptons)[0]

    def linear_singlets(self, S):
        Z = self.Z
        lin = []
        for f in Z.sing:
            o = f["name"]
            if o in S:
                continue
            k = self.min_order(S, [p - w for p, w in zip(Z.vec[o], Z.wvec)])
            if k is not None:
                lin.append((k + 1, o))
        return sorted(lin)

    def dflat_full(self, S):
        """FI-cancelling D-flat point with EVERY field of S nonzero (LP per field, then average)."""
        Z = self.Z
        A = np.array([[float(Z.vec[x][k]) for x in S] for k in range(1, Z.nq)])
        sgn = 1 if Z.tr0 > 0 else -1
        c0 = np.array([float(sgn * Z.vec[x][0]) for x in S])
        for i in range(len(S)):
            o = np.zeros(len(S))
            o[i] = -1
            rr = linprog(o, A_eq=np.vstack([A, c0]), b_eq=np.concatenate([np.zeros(A.shape[0]), [-1.0]]),
                         bounds=[(0, 10)] * len(S), method="highs")
            if not (rr.status == 0 and -rr.fun > 1e-9):
                return False
        return True

    def monomials(self, S, dmax=DMAX):
        """every W-allowed monomial in the fields of S with degree 2..dmax (exact integer test, all rules)."""
        Z = self.Z
        den = lcm(*[x.denominator for v in list(Z.vec.values()) + [Z.wvec] for x in v])
        A = np.array([[int(x * den) for x in Z.vec[s]] for s in S], dtype=np.int64)
        w = np.array([int(x * den) for x in Z.wvec], dtype=np.int64)
        nq, n, out = Z.nq, len(S), []

        def rec(i, left, e, tot):
            if i == n:
                d = dmax - left
                t = tot - w
                if d >= 2 and np.all(t[:nq] == 0) and np.all(t[nq:] % den == 0):
                    out.append((d, tuple(e)))
                return
            for k in range(left + 1):
                e.append(k)
                rec(i + 1, left - k, e, tot + k * A[i])
                e.pop()
        rec(0, dmax, [], np.zeros(A.shape[1], dtype=np.int64))
        return sorted(out)


def f_obstruction(mons):
    """lowest degree at which the exponent vectors become linearly dependent (None: independent to DMAX)."""
    rows = []
    for d, e in mons:
        rows.append(list(e))
        if sp.Matrix(rows).rank() < len(rows):
            return d
    return None


def main():
    blob = load()
    L, D = blob["ledger"], blob["disc"]
    sig = {(v["fi_header"][:8], tuple(sorted(Counter((f["dim"], f["k"]) for f in v["left"]).items())))
           for v in L.values()}
    A = dict(scan=blob["scan"], distinct_spectra=len(sig))
    print("A", {k: (v["tries"], v["standard_models"]) for k, v in blob["scan"].items()}, "distinct", len(sig), flush=True)
    B, Bsum = census(L, D)
    print("B", Bsum, flush=True)
    C, E = {}, {}
    Esum = Counter()
    for name in sorted(n for n, r in B.items() if r.get("verdict") == "COUNTEREXAMPLE"):
        b = name.split("|")[1]
        Z = P8.Z6IIR(name, {name: L[name]}, {b: D[b]}, {})
        T = Tools(Z)
        vac, supports = P8.z6ii23_vacua(Z)
        rows, seen = [], set()
        for S in supports:
            v = T.viable(S)
            rec = dict(support=S, mssm_parity=v)
            if v == "viable":
                rec.update(exotic_d=T.rank(S, "d", "bd"), higgs=T.rank(S, "bl", "l"),
                           W_S_lowest_order=T.min_order(S, [-w for w in Z.wvec]), linear_terms=T.linear_singlets(S))
                X = list(dict.fromkeys(list(S) + [o for _, o in rec["linear_terms"]]))
                if frozenset(X) not in seen:
                    seen.add(frozenset(X))
                    ex = dict(support=X, joint_parity=T.viable(X), exotic_d=T.rank(X, "d", "bd"))
                    if ex["joint_parity"] == "viable" and ex["exotic_d"][0] == ex["exotic_d"][1] > 0:
                        ex["dflat_full"] = T.dflat_full(X)
                        ex["higgs"] = T.rank(X, "bl", "l")
                        mons = T.monomials(X)
                        by = Counter(d for d, _ in mons)
                        ex["monomials_by_degree"] = {str(k): by[k] for k in sorted(by)}
                        ex["lowest_monomials"] = [list(e) for d, e in mons if d == min(by)]
                        ex["first_dependent_degree"] = f_obstruction(mons)
                        Esum["extended_vacua"] += 1
                        Esum["dflat_full"] += ex["dflat_full"]
                        Esum["higgs_massless_all_orders"] += ex["higgs"][0] == 0
                        Esum["higgs_massive"] += ex["higgs"][0] > 0
                        fd = ex["first_dependent_degree"]
                        Esum["independent_through_%d" % DMAX] += fd is None
                        Esum["min_first_dependent_degree"] = min(Esum.get("min_first_dependent_degree", 99), fd or 99)
                        print("E", b, len(X), ex["dflat_full"], ex["higgs"], ex["monomials_by_degree"], fd, flush=True)
                    else:
                        Esum["extended_rejected"] += 1
                    E.setdefault(name, []).append(ex)
            rows.append(rec)
        C[name] = dict(vacua=vac, supports=rows,
                       summary=dict(Counter(r["mssm_parity"] for r in rows)),
                       exotic_d_massless_all=all(r["exotic_d"][0] == 0 for r in rows if r["mssm_parity"] == "viable"),
                       W_S_zero_all=all(r["W_S_lowest_order"] is None for r in rows if r["mssm_parity"] == "viable"),
                       max_linear_order=max((r["linear_terms"][0][0] for r in rows if r.get("linear_terms")), default=None))
        print("C", b, C[name]["summary"], C[name]["exotic_d_massless_all"], C[name]["W_S_zero_all"],
              C[name]["max_linear_order"], flush=True)
    print("E", dict(Esum), flush=True)
    OUT.write_text(json.dumps(dict(pass_id=10979, dmax=DMAX, A=A, B=dict(summary=Bsum, models=B), C=C,
                                   E=dict(summary=dict(Esum), vacua=E)), indent=1, sort_keys=True, default=str))


if __name__ == "__main__":
    main()
