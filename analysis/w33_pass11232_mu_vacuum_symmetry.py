#!/usr/bin/env python3
"""Pass 11232: can a vacuum symmetry forbid the mu term while keeping the top Yukawa?  Exact test over all 215 models.

The paper's scorecard (sec09) lists the mu term as OPEN: "arises at order 3 or 4 in every model; none is mu-free through
order five.  Must come from a vacuum symmetry."  A vacuum symmetry is a symmetry left unbroken by the singlet vacuum
values that cancel the anomalous Fayet-Iliopoulos term.  This pass tests the continuous part exactly.

Setting (the frozen exact-rational ledger of Pass 10960, all 215 Z6-I / Z6-II models; its D-flat machinery reused):
  * vacuum = an FI-cancelling D-flat direction of the SM- and hidden-neutral singlets.  Every such direction is a
    nonnegative combination of extreme rays of the D-flat cone, and its support contains the support of a ray, so the
    unbroken symmetry is largest on a ray: rays suffice for an existence question.
  * the unbroken U(1)s of a vacuum S are the charge directions t with t . q_s = 0 for every s in S.  A coupling with
    total charge Q is forbidden TO ALL ORDERS by an unbroken U(1) iff Q is not in the rational span of {q_s : s in S}.
  * mu-protection of an up-type Higgs H_u = bl_i: every mu_ij = bl_i l_j is forbidden by the vacuum U(1)s, while some
    top Yukawa q_a bl_i bu_b is not (its charge lies in the span: necessary for it to appear at some order).
  * a mu term can then come only from breaking the protecting U(1)' (the U(1)' solution of the mu problem): we record
    whether some SM-singlet s' outside the vacuum has q_s' + q_mu in the span, i.e. an s' bl_i l_j coupling that would
    give mu = lambda <s'> if s' condenses at the TeV scale.

Limits, stated up front: only the continuous gauge U(1)s are tested.  The discrete space-group and R-symmetry rules
can forbid more, never less, so a "protected" verdict here is rigorous, while "not protected" leaves room for a discrete
vacuum symmetry (cf. Kappl et al. arXiv:0812.2120 and Pass 10968, where an unbroken R-symmetry also keeps exotics
massless).  F-flatness is not imposed (Pass 10960/10962 scope).  The unbroken U(1)' that protects mu is a massless Z'
in that vacuum: a phenomenological cost, recorded, not hidden.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass10960_matter_even_dflat_closure as P  # noqa: E402

OUT = ROOT / "data" / "w33_pass11232_mu_vacuum_symmetry.json"


def structural_rank(rows, cols, allowed):
    """generic rank of a mass matrix whose allowed entries are the given pairs (maximum bipartite matching)"""
    match = {}

    def aug(i, seen):
        for j in cols:
            if allowed(i, j) and j not in seen:
                seen.add(j)
                if j not in match or aug(match[j], seen):
                    match[j] = i
                    return True
        return False
    return sum(aug(i, set()) for i in rows)


class Span:
    """exact rational row space with membership test"""

    def __init__(self, vectors):
        self.rows = []          # (pivot, row) in echelon form
        for v in vectors:
            self.add(v)

    def reduce(self, v):
        v = list(v)
        for p, r in self.rows:
            if v[p] != 0:
                c = v[p] / r[p]
                v = [a - c * b for a, b in zip(v, r)]
        return v

    def add(self, v):
        w = self.reduce(v)
        piv = next((i for i, x in enumerate(w) if x != 0), None)
        if piv is not None:
            self.rows.append((piv, w))

    def __contains__(self, v):
        return all(x == 0 for x in self.reduce(v))

    def key(self):
        # canonical RREF
        rows = [list(r) for _, r in sorted(self.rows)]
        for i, r in enumerate(rows):
            p = next(j for j, x in enumerate(r) if x != 0)
            r[:] = [x / r[p] for x in r]
            for k, s in enumerate(rows):
                if k != i and s[p] != 0:
                    c = s[p]
                    s[:] = [a - c * b for a, b in zip(s, r)]
        return tuple(tuple(r) for r in sorted(rows))


def fields(model):
    left = []
    for f in model["left"]:
        dims = f["dim"].split(",")
        import re
        left.append(dict(name=f["name"], base=P.base_of(f["name"]),
                         trivial=all(abs(int(re.match(r"(-?\d+)", d).group(1))) == 1 and "adj" not in d for d in dims),
                         dimprod=abs(eval("*".join(re.match(r"(-?\d+)", d).group(1) for d in dims))),
                         q=[Fraction(x) for x in f["q"]]))
    return left


def sampled_vertices(C, M, samples=3000, seed=11232):
    """vertices of {a >= 0, sum a_i M_i = 0, C.a = -1} from random LP objectives; each re-solved exactly on its support"""
    import numpy as np
    import sympy as sp
    from scipy.optimize import linprog
    rng = np.random.default_rng(seed)
    n = len(C)
    Aeq = np.array([[float(M[i][j]) for i in range(n)] for j in range(len(M[0]))] + [[float(c) for c in C]])
    beq = np.zeros(len(Aeq))
    beq[-1] = -1
    supports, out = set(), []
    for _ in range(samples):
        r = linprog(rng.exponential(size=n), A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs-ds")
        if r.status != 0:
            continue
        sup = tuple(int(i) for i in np.flatnonzero(r.x > 1e-9))
        if sup in supports:
            continue
        supports.add(sup)
        A = sp.Matrix([[sp.Rational(M[i][j].numerator, M[i][j].denominator) for i in sup] for j in range(len(M[0]))]
                      + [[sp.Rational(C[i].numerator, C[i].denominator) for i in sup]])
        b = sp.Matrix([0] * len(M[0]) + [-1])
        try:
            x, params = A.gauss_jordan_solve(b)
        except ValueError:
            continue
        if params.shape[0]:
            continue                      # not a vertex support
        a = [Fraction(0)] * n
        for k, i in enumerate(sup):
            a[i] = Fraction(int(sp.fraction(x[k])[0]), int(sp.fraction(x[k])[1]))
        if all(a[i] > 0 for i in sup):
            assert all(sum((M[i][j] * a[i] for i in range(n)), Fraction(0)) == 0 for j in range(len(M[0])))
            assert P.dot(C, a) == -1
            out.append(a)
    return out


def analyse(name, model, max_types_exact=44, vertex_samples=3000):
    nu1 = model["nu1"]
    left = fields(model)
    tr = [sum((f["dimprod"] * f["q"][i] for f in left), Fraction(0)) for i in range(nu1)]
    rec = dict(anomalous=tr[0] != 0)
    if tr[0] == 0:
        return rec
    lab = [f for f in left if f["base"] in P.SMY]
    y = P.solve_affine([f["q"] for f in lab], [P.SMY[f["base"]] for f in lab])[0]
    sing = [f for f in left if f["trivial"] and P.dot(y, f["q"]) == 0]
    types = {}
    for f in sing:
        types.setdefault(tuple(f["q"]), []).append(f["name"])
    T = list(types)
    s = 1 if tr[0] > 0 else -1
    C = [s * t[0] for t in T]
    M = [list(t[1:]) for t in T]
    if P.farkas(C, M) is not None:
        rec["dflat"] = False
        return rec
    rec["dflat"] = True
    if len(T) <= max_types_exact:
        rays = [r for r in P.cone_rays(M) if P.dot(C, r) < 0]
        rec["method"] = "all extreme rays (exact, cdd)"
    else:
        rays = sampled_vertices(C, M, samples=vertex_samples)
        rec["method"] = f"{vertex_samples} LP-sampled vertices, each verified exactly (sampling: 'none found' is not proof)"
    rec["fi_rays"] = len(rays)
    hu = [f for f in left if f["base"] == "bl"]
    hd = [f for f in left if f["base"] == "l"]
    qs = [f for f in left if f["base"] == "q"]
    us = [f for f in left if f["base"] == "bu"]
    add = lambda *vs: [sum(x) for x in zip(*vs)]
    mu = {(i, j): add(a["q"], b["q"]) for i, a in enumerate(hu) for j, b in enumerate(hd)}
    top = {i: [add(qa["q"], a["q"], ub["q"]) for qa in qs for ub in us] for i, a in enumerate(hu)}
    seen = {}
    for r in rays:
        sup = [T[k] for k, v in enumerate(r) if v != 0]
        V = Span([list(t) for t in sup])
        key = V.key()
        if key in seen:
            seen[key]["rays"] += 1
            continue
        protected = []
        for i in range(len(hu)):
            if all(mu[i, j] not in V for j in range(len(hd))) and any(t in V for t in top[i]):
                gen = sorted({types[t][0] for t in T if t not in sup
                              for j in range(len(hd)) if add(list(t), mu[i, j]) in V})
                rest = [k for k in range(len(hu)) if k != i]
                drank = structural_rank(rest, range(len(hd)), lambda a, b: mu[a, b] in V)
                protected.append(dict(Hu=hu[i]["name"], mu_generating_singlets=gen[:12],
                                      mu_generating_count=len(gen),
                                      other_Hu_all_massive=drank == len(rest)))
        exotic = {}
        for b in ("q", "u", "d", "e", "v", "x"):
            A = [f for f in left if f["base"] == b]
            B = [f for f in left if f["base"] == "b" + b]
            if A and B:
                k = structural_rank(range(len(A)), range(len(B)), lambda a, c: add(A[a]["q"], B[c]["q"]) in V)
                exotic[b] = dict(n=len(A), nbar=len(B), massive_pairs=k, light_vectorlike=min(len(A), len(B)) - k)
        seen[key] = dict(rays=1, rank=len(V.rows), nu1=nu1, support=sorted(types[t][0] for t in sup),
                         exotics=exotic, vectorlike_all_massive=all(e["light_vectorlike"] == 0 for e in exotic.values()),
                         mu_forbidden_pairs=sum(1 for v in mu.values() if v not in V),
                         top_allowed_Hu=sum(1 for i in top if any(t in V for t in top[i])),
                         protected=protected)
    rec["vacuum_spans"] = len(seen)
    rec["mu_protected_vacua"] = sum(1 for v in seen.values() if v["protected"])
    rec["mu_protected_with_generator"] = sum(1 for v in seen.values()
                                             if any(p["mu_generating_count"] for p in v["protected"]))
    rec["mu_protected_clean"] = sum(1 for v in seen.values() if v["vectorlike_all_massive"] and any(
        p["mu_generating_count"] and p["other_Hu_all_massive"] for p in v["protected"]))
    rec["any_mu_forbidden"] = sum(1 for v in seen.values() if v["mu_forbidden_pairs"])
    rec["examples"] = [v for v in seen.values() if v["protected"]][:3]
    rec["Hu_count"], rec["Hd_count"] = len(hu), len(hd)
    return rec


def clean_detail(model, support_names):
    """the unbroken U(1)s of a vacuum: exact null space of the support charges; U(1)' = the part that charges mu"""
    import sympy as sp
    left = fields(model)
    by = {f["name"]: f for f in left}
    Qs = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in by[n]["q"]] for n in support_names])
    null = Qs.nullspace()
    lab = [f for f in left if f["base"] in P.SMY]
    y = P.solve_affine([f["q"] for f in lab], [P.SMY[f["base"]] for f in lab])[0]
    ch = lambda t, f: sum((sp.Rational(a.numerator, a.denominator) * b for a, b in zip(f["q"], t)), sp.Integer(0))
    gens = []
    for t in null:
        t = list(t)
        gens.append(dict(t=[str(v) for v in t],
                         is_hypercharge_multiple=all(ch(t, f) * P.SMY[lab[0]["base"]] == ch(t, lab[0]) * P.SMY[f["base"]]
                                                     for f in lab)))
    # a U(1)' orthogonal-in-charge-space choice: the null vector giving bl_1 + l_j nonzero charge
    out = dict(unbroken_dim=len(null), generators=gens)
    tp = None
    for t in null:
        if any(ch(list(t), by["bl_1"]) + ch(list(t), f) != 0 for f in left if f["base"] == "l"):
            tp = list(t)
            break
    if tp is not None:
        out["U1prime_charges"] = {f["name"]: str(ch(tp, f)) for f in left
                                  if f["base"] in ("bl", "l", "q", "bu", "bd", "be") or f["name"] in support_names}
    return out


def main():
    ledger, sha = P.load_ledger()
    names = sorted(ledger)
    if len(sys.argv) > 1:
        names = [n for n in names if any(a in n for a in sys.argv[1:])]
    res, t0 = {}, time.time()
    for n in names:
        t = time.time()
        res[n] = analyse(n, ledger[n])
        r = res[n]
        print(n, {k: v for k, v in r.items() if k != "examples"}, f"{time.time() - t:.1f}s", flush=True)
    dflat = [k for k, v in res.items() if v.get("dflat")]
    summary = dict(models=len(res), anomalous=sum(1 for v in res.values() if v["anomalous"]), dflat=len(dflat),
                   models_with_mu_protected_vacuum=sorted(k for k in dflat if res[k]["mu_protected_vacua"]),
                   models_with_protected_and_generator=sorted(k for k in dflat if res[k]["mu_protected_with_generator"]),
                   models_with_clean_protected_vacuum=sorted(k for k in dflat if res[k]["mu_protected_clean"]),
                   models_with_any_mu_forbidden=sum(1 for k in dflat if res[k]["any_mu_forbidden"]),
                   ledger_sha256=sha, seconds=round(time.time() - t0, 1))
    details = {}
    if len(sys.argv) == 1:
        for k in summary["models_with_clean_protected_vacuum"]:
            for ex in res[k]["examples"]:
                if ex["vectorlike_all_massive"]:
                    details[k] = dict(support=ex["support"], **clean_detail(ledger[k], ex["support"]))
                    break
    OUT.write_text(json.dumps(dict(pass_id=11232, summary=summary, clean_vacua=details, models=res), indent=1,
                              default=str))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
