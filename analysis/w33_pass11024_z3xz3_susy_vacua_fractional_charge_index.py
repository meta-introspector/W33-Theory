#!/usr/bin/env python3
"""Pass 11024: the Z3xZ3 parity vacua are supersymmetric -- and carry massless fractional charges.

Pass 10979 left one door: 28 Z3xZ3 vacua that keep matter parity through the Fayet-Iliopoulos term, whose
superpotential has independent exponent vectors through degree 11.  This pass settles them.

A. All-orders structure.  On every vacuum S the monomials invariant under all continuous charges form a FREE
   monoid on k = 1 or 2 composites p_i = phi^{r_i} (disjoint field products of degree 3 or 4; lattice index 1),
   and p^a is allowed in W iff a_1 + ... + a_k = 1 (mod 3) -- a condition checked on a full period, hence exact to
   all orders.  So W|_S = F(p) with F(w p) = w F(p), w^3 = 1: a Z3 R-symmetry acting on the composites.  With
   every field nonzero, F_phi = 0 <=> grad_p F = 0.  This explains Pass 10979's 'degree 12' (= 4 x 3) and its
   prediction (degree 16 for the quartic composites) is met.
B. SUSY vacua exist.  For random complex O(1) couplings (truncation |a| <= 4), grad F = 0 has the Bezout number of
   roots (3 for k = 1, ~9 for k = 2), each lifts to an exact F = D = 0 point (FI term xi = g^2 Tr Q0/192 pi^2,
   Holotrade cfdc1f2 convention) with every field nonzero: 28 of 28 vacua.  The VEV is set by the coupling
   ratio, |p|^3 ~ |c_1/c_4|; with O(1) couplings max|phi| ~ 0.9-1.1 M_s.  (Pass 10979's numerical search started
   at the FI scale and missed them; its exact statement -- no balance below degree 12 -- stands.)
C. Massless fractional charges.  A mass term needs a conjugate (SM x hidden) pair, so the excess per class is
   massless in any vacuum preserving the hidden gauge group.  In both parity models the excess is 48 states of
   electric charge +-1/3: six hidden-SINGLET colour-singlet multiplets plus SU(4) 4, 4bar and 6 multiplets.  All
   28 vacua are hidden singlets, so all 28 carry six exactly massless charge-1/3 fermions.  Both hidden groups are
   infrared free (b_SU(4) = -7/-4, b_SU(3)' = -16/-13), so nothing confines them.  Control: the same index
   vanishes in most Z6 models (it is a property of these spectra, not of the method).
D. Hidden condensates.  Adding every hidden meson/baryon (Luty-Taylor) as a composite singlet to the exact parity
   x FI census; auditing which realizable FI supports use a composite; and screening every SU(4) meson added to
   each of the 28 vacua for joint parity, D-flatness and the pairing of the six unconfined states.
"""
from __future__ import annotations

import gzip
import importlib.util
import itertools
import json
import math
from math import lcm
import re
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.optimize import least_squares, linear_sum_assignment, linprog

ROOT = Path(__file__).resolve().parents[1]


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / "analysis" / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


P79 = _load("p10979", "w33_pass10979_z3xz3_parity_vacua_f_obstruction.py")
EXO = _load("exact_orders", "w33_exact_monomial_orders.py")
FastParity = _load("fast_parity", "w33_fast_parity_realizability.py").FastParity
P8, P74, M = P79.P8, P79.P74, P79.M
P = M.P
CERT79 = ROOT / "data" / "w33_pass10979_z3xz3_parity_vacua_f_obstruction.json"
OUT = ROOT / "data" / "w33_pass11024_z3xz3_susy_vacua_fractional_charge_index.json"
CONTROLS = {"Z2xZ6-I": ROOT / "data" / "w33_pass10968_z2xz6_ledger.json.gz",
            "Z3xZ6": ROOT / "data" / "w33_pass10974_z3xz6_ledger.json.gz",
            "Z6xZ6": ROOT / "data" / "w33_pass10974_z6xz6_ledger.json.gz"}
G2 = 0.5
TDYN = {2: Fraction(1, 2), 3: Fraction(1, 2), 4: Fraction(1, 2), 6: Fraction(1), 8: Fraction(3), 15: Fraction(4)}


def model_Z(blob, name):
    b = name.split("|")[1]
    return P8.Z6IIR(name, {name: blob["ledger"][name]}, {b: blob["disc"][b]}, {})


# ---------------------------------------------------------------- A. all-orders structure
def structure(Z, S):
    nq = Z.nq
    Q = sp.Matrix([[sp.Rational(Z.vec[x][c].numerator, Z.vec[x][c].denominator) for x in S] for c in range(nq)])
    ker = Q.nullspace()
    k = len(ker)
    B = sp.Matrix.hstack(*[v * sp.ilcm(*[sp.fraction(e)[1] for e in v]) for v in ker])
    if k == 1:
        v = B[:, 0]
        v = -v if all(e <= 0 for e in v) else v
        assert all(e >= 0 for e in v)
        rays = [v / sp.igcd(*[int(e) for e in v])]
    else:
        assert k == 2
        cands = set()
        for i in range(len(S)):
            t = sp.Matrix([B[i, 1], -B[i, 0]])
            for s in (1, -1):
                r = B * (s * t)
                if any(r) and all(e >= 0 for e in r):
                    g = sp.igcd(*[int(e) for e in r])
                    cands.add(tuple(int(e) // g for e in r))
        rays = [sp.Matrix(r) for r in sorted(cands)]
    assert len(rays) == k
    R = sp.Matrix.hstack(*rays)
    minors = [R.extract(list(c), list(range(k))).det() for c in itertools.combinations(range(len(S)), k)]
    index = abs(sp.igcd(*[int(m) for m in minors if m != 0]))
    dims = len(Z.wvec)

    def allowed(a):
        e = [sum(a[i] * int(R[j, i]) for i in range(k)) for j in range(len(S))]
        tot = [sum(Fraction(e[j]) * Z.vec[S[j]][c] for j in range(len(S))) - Z.wvec[c] for c in range(dims)]
        return all(t == 0 for t in tot[:nq]) and all(t.denominator == 1 for t in tot[nq:])
    period = list(itertools.product(range(3), repeat=k))
    allowed_mod3 = [a for a in period if allowed(a)]
    periodic = all(allowed(tuple(x + 3 * y for x, y in zip(a, sh))) == allowed(a)
                   for a in period for sh in itertools.product(range(3), repeat=k))
    z3_graded = sorted(allowed_mod3) == sorted(a for a in period if sum(a) % 3 == 1)
    disjoint = all(sum(1 for i in range(k) if R[j, i] != 0) <= 1 for j in range(len(S)))
    return dict(k=k, rays=[[int(x) for x in R[:, i]] for i in range(k)], ray_degrees=[int(sum(R[:, i])) for i in range(k)],
                lattice_index=int(index), allowed_mod3=[list(a) for a in allowed_mod3], periodic=periodic,
                z3_graded=z3_graded, disjoint_supports=disjoint)


# ---------------------------------------------------------------- B. explicit SUSY roots
def susy_roots(Z, S, st, seeds=4):
    k, R = st["k"], np.array(st["rays"], dtype=float)
    n = len(S)
    Q = np.array([[float(Z.vec[x][c]) for x in S] for c in range(Z.nq)])
    xi = G2 * abs(float(Z.tr0)) / (192 * math.pi ** 2)
    sgn = 1 if Z.tr0 > 0 else -1
    terms = [(1,), (4,)] if k == 1 else [(1, 0), (0, 1)] + [(i, 4 - i) for i in range(5)]

    def dF(p, coef):
        g = np.zeros(k, complex)
        for a, c in zip(terms, coef):
            for i in range(k):
                if a[i]:
                    e = list(a); e[i] -= 1
                    g[i] += c * a[i] * np.prod(p ** np.array(e))
        return g
    out = []
    for seed in range(seeds):
        rng = np.random.default_rng(seed)
        coef = rng.normal(size=len(terms)) + 1j * rng.normal(size=len(terms))
        roots = []
        for _ in range(60):
            z0 = rng.normal(size=2 * k)
            f = lambda z: (lambda g: np.concatenate([g.real, g.imag]))(dF(z[:k] + 1j * z[k:], coef))
            sol = least_squares(f, z0, xtol=1e-15, ftol=1e-15, gtol=1e-15)
            p = sol.x[:k] + 1j * sol.x[k:]
            if np.linalg.norm(f(sol.x)) < 1e-11 and np.all(np.abs(p) > 1e-6) and not any(np.allclose(p, r, atol=1e-8) for r in roots):
                roots.append(p)
        for p in roots:
            def eqs(s):
                D = Q @ np.exp(2 * s)
                D[0] += sgn * xi
                return np.concatenate([R @ s - np.log(np.abs(p)), D])
            best = None
            for _ in range(30):
                sol = least_squares(eqs, rng.normal(size=n) + np.log(0.5), xtol=1e-15, ftol=1e-15, gtol=1e-15)
                r = np.linalg.norm(eqs(sol.x))
                if best is None or r < best[0]:
                    best = (r, sol.x)
                if r < 1e-12:
                    break
            s = best[1]
            th = np.linalg.lstsq(R, np.angle(p), rcond=None)[0]
            phi = np.exp(s) * np.exp(1j * th)
            pv = np.array([np.prod(phi ** R[i]) for i in range(k)])
            g = dF(pv, coef)
            Fphi = np.array([sum(g[i] * R[i, j] * pv[i] / phi[j] for i in range(k)) for j in range(n)])
            D = Q @ np.abs(phi) ** 2
            D[0] += sgn * xi
            out.append(dict(F=float(np.linalg.norm(Fphi)), D=float(np.linalg.norm(D)), min_abs=float(np.abs(phi).min()),
                            max_abs=float(np.abs(phi).max())))
    ok = [o for o in out if o["F"] < 1e-9 and o["D"] < 1e-9 and o["min_abs"] > 1e-6]
    mx = sorted(o["max_abs"] for o in ok)
    return dict(xi=xi, roots=len(out), verified=len(ok),
                max_abs_median=float(np.median(mx)) if mx else None, max_abs_range=[mx[0], mx[-1]] if mx else None)


# ---------------------------------------------------------------- C. hidden-refined index
def dimof(tok):
    """dimension of a representation token ('-3', '8v', '16', ...)."""
    return int(re.match(r"-?(\d+)", tok).group(1))


def conj_token(tok, present, pseudo_real):
    """conjugate representation: '-r' <-> 'r' when both occur in the slot; otherwise self-conjugate."""
    if pseudo_real or dimof(tok) == 1:
        return tok
    if tok.startswith("-"):
        return tok[1:]
    return "-" + tok if "-" + tok in present else tok


def slots(model):
    left = model["left"]
    dims = {f["name"]: f["dim"].split(",") for f in left}
    q0 = next(f for f in left if P.base_of(f["name"]) == "q")
    l0 = next(f for f in left if P.base_of(f["name"]) == "l")
    qd, ld = dims[q0["name"]], dims[l0["name"]]
    c3 = [i for i, x in enumerate(qd) if dimof(x) == 3][0]
    w2 = [i for i, x in enumerate(ld) if dimof(x) == 2 and dimof(qd[i]) == 2][0]
    return dims, qd, c3, w2


def chiral_index(model):
    left = model["left"]
    q = {f["name"]: [Fraction(x) for x in f["q"]] for f in left}
    lab = [f for f in left if P.base_of(f["name"]) in P.SMY]
    y = P.solve_affine([q[f["name"]] for f in lab], [P.SMY[P.base_of(f["name"])] for f in lab])[0]
    dims, qd, c3, w2 = slots(model)
    nslot = len(qd)
    present = defaultdict(set)
    for d in dims.values():
        for i, x in enumerate(d):
            present[i].add(x)
    conj = lambda i, x: conj_token(x, present[i], i == w2)
    cls = Counter((tuple(dims[f["name"]]), P.dot(y, q[f["name"]])) for f in left)
    cq = qd[c3]
    cqb = conj_token(cq, present[c3], False)
    hs = lambda col, wk: tuple(col if i == c3 else (wk if i == w2 else "1") for i in range(nslot))
    mssm = {(hs(cq, "2"), Fraction(1, 6)): 3, (hs(cqb, "1"), Fraction(-2, 3)): 3, (hs(cqb, "1"), Fraction(1, 3)): 3,
            (hs("1", "2"), Fraction(-1, 2)): 3, (hs("1", "1"), Fraction(1)): 3}
    ck = lambda K: (tuple(conj(i, x) for i, x in enumerate(K[0])), -K[1])
    out, frac_states, all_states = [], 0, 0
    for K in sorted(set(cls) | {ck(K) for K in cls}, key=str):
        d, Y = K
        col, wk = dimof(d[c3]), dimof(d[w2])
        if (col == 1 and wk == 1 and Y == 0) or ck(K) == K:
            continue
        ex = cls.get(K, 0) - cls.get(ck(K), 0) - mssm.get(K, 0) + mssm.get(ck(K), 0)
        if ex <= 0:
            continue
        hid = math.prod(dimof(x) for i, x in enumerate(d) if i not in (c3, w2))
        charges = [Y + Fraction(t_, 2) for t_ in range(-(wk - 1), wk, 2)]
        if col == 1:
            frac = any(c.denominator != 1 for c in charges)
        elif col == 3:
            tgt = Fraction(2, 3) if d[c3] == cq else Fraction(1, 3)
            frac = any((c - tgt).denominator != 1 for c in charges)
        else:
            frac = False
        n = ex * hid * col * wk
        out.append(dict(rep=list(d), Y=str(Y), excess=ex, hidden_dim=hid, fractional=frac, states=n,
                        hidden_singlet=hid == 1))
        all_states += n
        frac_states += n if frac else 0
    return dict(classes=out, fractional_states=frac_states, exotic_states=all_states,
                unconfinable_fractional=sum(c["states"] for c in out if c["fractional"] and c["hidden_singlet"]))


def hidden_beta(model):
    """b = 3 C2(G) - sum T(R) for each hidden SU(N) slot (N = smallest complex rep present; None if not SU(N))."""
    dims, qd, c3, w2 = slots(model)
    out = {}
    for i in range(len(qd)):
        if i in (c3, w2):
            continue
        vals = {d[i] for d in dims.values()}
        cplx = sorted(dimof(v) for v in vals if not v.startswith("-") and dimof(v) > 1 and "-" + v in vals)
        if not cplx or any(dimof(v) not in TDYN and dimof(v) != 1 for v in vals):
            out[str(i)] = dict(N=None, reps=sorted(vals))
            continue
        N = cplx[0]
        s = Fraction(0)
        for d in dims.values():
            a = dimof(d[i])
            if a == 1:
                continue
            s += TDYN[a] * math.prod(dimof(x) for j, x in enumerate(d) if j != i)
        out[str(i)] = dict(N=N, sum_T=str(s), b=str(3 * N - s), infrared_free=3 * N - s < 0)
    return out


def break_slots(model, slots_):
    """the spectrum with the given hidden factors broken: their multiplets split into singlet components."""
    new = []
    for f in model["left"]:
        d = f["dim"].split(",")
        k = 1
        for s_ in slots_:
            k *= dimof(d[s_])
            d[s_] = "1"
        new += [dict(f, dim=",".join(d))] * k
    return dict(model, left=new)


def curing_slots(model):
    """hidden factors whose breaking alone makes the fractional-charge index vanish."""
    dims, qd, c3, w2 = slots(model)
    hid = [i for i in range(len(qd)) if i not in (c3, w2)]
    table = {}
    for i in hid:
        ix = chiral_index(break_slots(model, [i]))
        table[str(i)] = dict(fractional=ix["fractional_states"], unconfinable=ix["unconfinable_fractional"])
    ix = chiral_index(break_slots(model, hid))
    table["all"] = dict(fractional=ix["fractional_states"], unconfinable=ix["unconfinable_fractional"])
    return [i for i in hid if table[str(i)]["fractional"] == 0], table


# ---------------------------------------------------------------- D. hidden condensates (Luty-Taylor composites)
ORIG_FIELDS = M.fields


def hidden_composites(model):
    """holomorphic invariants of one hidden factor built from SM-singlet fields charged under that factor only."""
    left = model["left"]
    q = {f["name"]: [Fraction(x) for x in f["q"]] for f in left}
    lab = [f for f in left if P.base_of(f["name"]) in P.SMY]
    y = P.solve_affine([q[f["name"]] for f in lab], [P.SMY[P.base_of(f["name"])] for f in lab])[0]
    dims, qd, c3, w2 = slots(model)
    hid = [i for i in range(len(qd)) if i not in (c3, w2)]
    present = defaultdict(set)
    for d in dims.values():
        for i, x in enumerate(d):
            present[i].add(x)
    group = {}
    for i in hid:
        cplx = sorted(dimof(v) for v in present[i] if not v.startswith("-") and dimof(v) > 1 and "-" + v in present[i])
        group[i] = cplx[0] if cplx else (2 if "2" in present[i] else None)
    by = defaultdict(list)
    for f in left:
        n, d = f["name"], dims[f["name"]]
        if dimof(d[c3]) != 1 or dimof(d[w2]) != 1 or P.dot(y, q[n]) != 0:
            continue
        nz = [i for i in hid if dimof(d[i]) != 1]
        if len(nz) == 1:
            tok = d[nz[0]]
            r = -dimof(tok) if tok.startswith("-") else (dimof(tok) if tok.isdigit() else None)
            if r is None:
                continue
            by[(nz[0], r)].append(n)
    comps, skipped = [], Counter()
    for (i, r), fs in sorted(by.items()):
        N = group[i]
        if N is not None and N >= 3 and r == N:
            comps += [("meson", (a, b), i) for a in fs for b in by.get((i, -N), [])]
            comps += [("baryon", c, i) for c in itertools.combinations(fs, N)]
        elif N is not None and N >= 3 and r == -N:
            comps += [("antibaryon", c, i) for c in itertools.combinations(fs, N)]
        elif N == 2 and r == 2:
            comps += [("su2pair", c, i) for c in itertools.combinations(fs, 2)]
        elif N == 4 and r == 6:
            comps += [("six2", c, i) for c in itertools.combinations_with_replacement(fs, 2)]
        else:
            skipped[f"{i}:{r}"] += len(fs)
    return comps, dict(skipped)


QUADRATIC = ("meson", "su2pair", "six2")


def augment(model, disc, kinds=QUADRATIC):
    """disc + composite singlets (summed U(1), space-group and R charges) and a fields() patch adding them with
    zero multiplicity (so Tr Q is unchanged).  Default scope: the QUADRATIC hidden invariants (as Pass 10962 Tier B);
    baryons enlarge the discrete-choice search beyond reach and are left as named residue."""
    comps, skipped = hidden_composites(model)
    all_kinds = dict(Counter(c[0] for c in comps))
    comps = [c for c in comps if c[0] in kinds]
    skipped = dict(skipped, **{"not_in_scope:" + k: v for k, v in all_kinds.items() if k not in kinds})
    d = json.loads(json.dumps(disc))
    orders, nR = disc["nonR_orders"], len(disc["R"])
    qmap = {f["name"]: [Fraction(x) for x in f["q"]] for f in model["left"]}
    extra = []
    for j, (kind, cs, slot) in enumerate(comps):
        cname = f"n_{90000 + j}"
        d["fields"][cname] = dict(susy=2, nonR=[str(sum(Fraction(disc["fields"][c]["nonR"][t]) for c in cs)) for t in range(len(orders))],
                                  R=[str(sum(Fraction(disc["fields"][c]["R"][t]) for c in cs)) for t in range(nR)])
        extra.append(dict(name=cname, kind=kind, of=list(cs), slot=slot,
                          q=[sum(qmap[c][t] for c in cs) for t in range(len(qmap[cs[0]]))]))

    def fields_aug(mdl, dsc):
        ords, out = ORIG_FIELDS(mdl, dsc)
        for e in extra:
            ch = dsc["fields"][e["name"]]["nonR"]
            out.append(dict(name=e["name"], base="n", k=int(Fraction(ch[0])) % ords[0], trivial=True, dimprod=0,
                            q=list(e["q"]), ext=[Fraction(c) / o for c, o in zip(ch, ords)]))
        return ords, out
    return d, fields_aug, extra, skipped


_POOL = None


def _realizable_job(args):
    vecs, odd = args
    return M.realizable(vecs, odd)


def pmap_realizable(jobs):
    """M.realizable over many independent (vecs, odd) jobs, in parallel when a pool is available."""
    if _POOL is None:
        return [M.realizable(v, o) for v, o in jobs]
    return _POOL.map(_realizable_job, jobs, chunksize=4)


def _ray_job(args):
    """depth-first realizability search for one FI ray (exact; FastParity engine, validated against Pass 10967's
    Smith-normal-form test): returns the chosen vectors or None."""
    choices, odd_vecs, even_vecs, nq, orders = args
    F = FastParity(nq, orders, odd_vecs, even_vecs)

    def dfs(i, chosen):
        # realizability is monotone (more even constraints only remove parity elements): prune partial choices
        if chosen and not F.realizable(chosen):
            return None
        if i == len(choices):
            return chosen
        for ch in choices[i]:
            r_ = dfs(i + 1, chosen + [ch])
            if r_ is not None:
                return r_
        return None
    return dfs(0, [])


def condensate_census(name, model, disc):
    """exact parity x FI census with the hidden composites added (the Pass 10967 census, mirrored so that each
    realizable FI-cancelling extreme ray also reports whether its support needs a composite-only U(1) type)."""
    d, fields_aug, extra, skipped = augment(model, disc)
    RN = [(int(x.split(":")[0]), Fraction(x.split(":")[1])) for x in disc["R"]]
    d2 = json.loads(json.dumps(d))
    d2["nonR_orders"] = d["nonR_orders"] + [n for n, _ in RN]
    for v in d2["fields"].values():
        v["nonR"] = v["nonR"] + v["R"]
    orders, left = fields_aug(model, d2)
    nq, nd = len(left[0]["q"]), len(orders)
    wvec = [Fraction(0)] * (nq + len(disc["nonR_orders"])) + [w / n for n, w in RN]
    real = lambda vecs, odd: M.realizable(vecs + [wvec], odd)
    vec = lambda f: f["q"] + f["ext"]
    virtual = [[Fraction(0)] * nq + [Fraction(int(i == j)) for j in range(nd)] for i in range(nd)]
    tr0 = sum((f["dimprod"] * f["q"][0] for f in left), Fraction(0))
    lab = [f for f in left if f["base"] in P.SMY]
    y = P.solve_affine([f["q"] for f in lab], [P.SMY[f["base"]] for f in lab])[0]
    sing = [f for f in left if f["trivial"] and P.dot(y, f["q"]) == 0]
    mv = [list(k) for k in dict.fromkeys(tuple(vec(f)) for f in left if f["base"] in P.MATTER)]
    base, odd = mv + virtual, list(range(len(mv)))
    out = dict(composites=len(extra), kinds=dict(Counter(e["kind"] for e in extra)), skipped=skipped)
    if not real(base, odd):
        out["census"] = dict(parity_exists=False)
        return out
    comp_names = {e["name"] for e in extra}
    ext_of = {e["name"]: e for e in extra}
    sv = list(dict.fromkeys(tuple(vec(f)) for f in sing))
    Fp = FastParity(nq, orders, mv, virtual + [wvec])   # validated against realizable() (Pass 11024 E)
    flags = [Fp.realizable([list(k)]) for k in sv]
    ever = [k for k, ok in zip(sv, flags) if ok]
    real_types = {tuple(vec(f))[:nq] for f in sing if f["name"] not in comp_names}
    utypes = sorted({k[:nq] for k in ever})
    s = 1 if tr0 > 0 else -1
    lam = P.farkas([s * u[0] for u in utypes], [list(u[1:]) for u in utypes]) if utypes else []
    if lam is not None:
        out["census"] = dict(parity_exists=True, verdict="closed: Farkas on the even-able singlets and composites")
        return out
    rays = P.cone_rays([list(u[1:]) for u in utypes])
    neg = [r for r in rays if sum((s * utypes[i][0] * r[i] for i in range(len(utypes))), Fraction(0)) < 0]
    hits = with_comp = 0
    found = []
    byvec = {}
    for f in sing:
        byvec.setdefault(tuple(vec(f)), f["name"])
    ext_of = {e["name"]: e for e in extra}
    cure, cure_table = curing_slots(model)
    curing_vecs = {tuple(vec(f)) for f in sing if f["name"] in comp_names and ext_of[f["name"]]["slot"] in cure}
    curing_types = {k[:nq] for k in curing_vecs if k in set(ever)}
    sups = [[utypes[i] for i in range(len(utypes)) if r[i] != 0] for r in neg]
    jobs, owners = [], []
    for ri, sup in enumerate(sups):
        for pos, u in enumerate(sup):
            if u not in curing_types:
                continue
            ch = [[list(k) for k in ever if k[:nq] == w] for w in sup]
            ch[pos] = [list(k) for k in ever if k[:nq] == u and k in curing_vecs]
            jobs.append((ch, mv, virtual + [wvec], nq, orders))
            owners.append(ri)
    print("   census", name.split("|")[1], "utypes", len(utypes), "FI rays", len(neg), "curing slots", cure,
          "rays with a curing condensate", len(set(owners)), "jobs", len(jobs), flush=True)
    results = _POOL.map(_ray_job, jobs, chunksize=1) if _POOL is not None else [_ray_job(j) for j in jobs]
    done = set()
    for ri, got in zip(owners, results):
        if got is not None and ri not in done:
            done.add(ri)
            hits += 1
            with_comp += 1
            found.append(dict(needs_composite=True, fields=[byvec[tuple(v)] for v in got]))
    out["curing_slots"] = cure
    out["index_with_one_factor_broken"] = cure_table
    out["rays_with_a_curing_condensate"] = len(set(owners))
    out["census"] = dict(parity_exists=True, fi_cancelling_rays=len(neg), realizable_curing_rays=hits,
                         realizable_rays_needing_a_composite=with_comp,
                         verdict="COUNTEREXAMPLE" if hits else "closed: every FI-cancelling ray obstructed")
    out["realizable_supports"] = [dict(needs_composite=s_["needs_composite"],
                                       fields=[dict(name=n_, kind=ext_of[n_]["kind"], of=ext_of[n_]["of"]) if n_ in ext_of
                                               else dict(name=n_) for n_ in s_["fields"]]) for s_ in found]
    out["_vecs"] = {n_: [str(x) for x in v] for v, n_ in byvec.items()}
    return out


def composite_vacua_parity(blob, name, census):
    """MSSM-viable parity (exact matcher) on every realizable FI support that needs a hidden composite."""
    Z = model_Z(blob, name)
    T = P79.Tools(Z)
    # the superpotential's own charge must be parity-even (P.w in Z): a non-R matter parity
    realz = lambda ev, od: M.realizable([list(t) for t in od] + [list(t) for t in ev] + [list(t) for t in Z.virtual]
                                        + [Z.wvec], list(range(len(od))))
    res = Counter()
    rows = []
    for s_ in census.get("realizable_supports", []):
        if not s_["needs_composite"]:
            continue
        ev = [[Fraction(x) for x in census["_vecs"][f["name"]]] for f in s_["fields"]]
        v = P8.match_search(realz, ev, T.species, T.leptons)[0]
        res[v] += 1
        rows.append(dict(fields=s_["fields"], mssm_parity=v))
    return dict(summary=dict(res), supports=rows)


def exact_for(Z, S):
    return EXO.ExactOrders([Z.vec[s] for s in S], Z.nq,
                           [x.denominator for v in list(Z.vec.values()) + [Z.wvec] + Z.virtual for x in v])


def assignment_rank(B):
    if B.size == 0:
        return 0
    r, c = linear_sum_assignment(-B)
    return int(B[r, c].sum())


def su4_rescue(blob, name, vacua):
    """every SU(4) meson a.b (a in 4, b in 4bar; SM x SU(3)' singlets) added to each vacuum: joint parity (necessary
    Smith-normal-form test), full-support FI D-flatness, and the pairing of the unconfined fractional class with its
    same-class partners plus the SU(3)-singlet components of the SU(4)-charged partners (exact orders, all orders)."""
    Z = model_Z(blob, name)
    T = P79.Tools(Z)
    lab = [f for f in Z.left if f["base"] in P.SMY]
    y = P.solve_affine([f["q"] for f in lab], [P.SMY[f["base"]] for f in lab])[0]
    dim = lambda n: [int(x) for x in Z.raw[n]["dim"].split(",")]
    Y = {f["name"]: P.dot(y, f["q"]) for f in Z.left}
    n4 = [f["name"] for f in Z.left if f["base"] == "n" and dim(f["name"]) == [4, 1, 1, 1]]
    n4b = [f["name"] for f in Z.left if f["base"] == "n" and dim(f["name"]) == [-4, 1, 1, 1]]
    plus = [f["name"] for f in Z.left if dim(f["name"]) == [1, 1, 1, 1] and Y[f["name"]] == Fraction(1, 3)]
    minus = [f["name"] for f in Z.left if dim(f["name"]) == [1, 1, 1, 1] and Y[f["name"]] == Fraction(-1, 3)]
    X, sgn_ = (minus, -1) if len(minus) > len(plus) else (plus, 1)
    same = plus if sgn_ == -1 else minus
    p4 = [f["name"] for f in Z.left if dim(f["name"]) == [4, 1, 1, 1] and Y[f["name"]] == -sgn_ * Fraction(1, 3)]
    p4b = [f["name"] for f in Z.left if dim(f["name"]) == [-4, 1, 1, 1] and Y[f["name"]] == -sgn_ * Fraction(1, 3)]
    need_su4 = len(X) - len(same)
    sgn = 1 if Z.tr0 > 0 else -1
    # the superpotential's own charge must be parity-even (P.w in Z): a non-R matter parity
    realz = lambda ev, od: M.realizable([list(t) for t in od] + [list(t) for t in ev] + [list(t) for t in Z.virtual]
                                        + [Z.wvec], list(range(len(od))))
    tally, rows = Counter(), []
    for v in vacua:
        S = v["support"]
        vt = Counter()
        E0 = exact_for(Z, S)
        base = assignment_rank(np.array([[int(E0.coupling([Z.vec[xn], Z.vec[cn]], Z.wvec) is not None) for cn in same]
                                         for xn in X]))
        for a in n4:
            for b in n4b:
                cvec = [x + z for x, z in zip(Z.vec[a], Z.vec[b])]
                if not Z.R([Z.vec[s] for s in S] + [cvec]):
                    vt["parity"] += 1
                    continue
                vecs = [Z.vec[s] for s in S] + [cvec]
                A = np.array([[float(q[k]) for q in vecs] for k in range(1, Z.nq)])
                c0 = np.array([float(sgn * q[0]) for q in vecs])
                full = True
                for i in range(len(vecs)):
                    o = np.zeros(len(vecs))
                    o[i] = -1
                    rr = linprog(o, A_eq=np.vstack([A, c0]), b_eq=np.concatenate([np.zeros(A.shape[0]), [-1.0]]),
                                 bounds=[(0, 10)] * len(vecs), method="highs")
                    if not (rr.status == 0 and -rr.fun > 1e-9):
                        full = False
                        break
                if not full:
                    vt["dflat"] += 1
                    continue
                E = exact_for(Z, list(S) + [a, b])
                Bsu4 = np.array([[int(E.coupling([Z.vec[xn], Z.vec[cn], Z.vec[b if cn in p4 else a]], Z.wvec) is not None)
                                  for cn in p4 + p4b] for xn in X])
                r4 = assignment_rank(Bsu4)
                if r4 < need_su4:
                    vt["su4_block_rank_%d_of_%d" % (r4, need_su4)] += 1
                    continue
                Bs = np.array([[int(E.coupling([Z.vec[xn], Z.vec[cn]], Z.wvec) is not None) for cn in same] for xn in X])
                rk = assignment_rank(np.hstack([Bs, Bsu4]))
                if rk < len(X):
                    vt["full_rank_%d_of_%d" % (rk, len(X))] += 1
                    continue
                par = P8.match_search(realz, vecs, T.species, T.leptons)[0]
                vt["all_paired_mssm_parity_" + par] += 1
        rows.append(dict(support=S, tally=dict(vt), baseline_rank=base, massless_unconfined_fractional=len(X) - base))
        tally.update(vt)
        tally["baseline_massless_min"] = min(tally.get("baseline_massless_min", 999), len(X) - base)
    return dict(n4=len(n4), n4bar=len(n4b), excess=len(X), same_class=len(same), su4_partners=[len(p4), len(p4b)],
                need_su4=need_su4, tally=dict(tally), vacua=rows)


# ---------------------------------------------------------------- E. exact certificate of the 10978-10979 ILPs
def _rank_orders(O):
    B = (O < 999).astype(int)
    rk = assignment_rank(B)
    reach = None
    for k in sorted(set(O.flatten())):
        if k == 999:
            break
        if assignment_rank((O <= k).astype(int)) == rk:
            reach = int(k)
            break
    return rk, reach


def certify_10979(blob, cert):
    diffs = checked = 0

    def rank(Z, E, Xb, Xs):
        R_, C_ = Z.lab(Xb), Z.lab(Xs)
        if not R_:
            return (0, 0, None)
        O = np.full((len(R_), len(C_)), 999)
        for i, a in enumerate(R_):
            for j, c in enumerate(C_):
                o = E.coupling([Z.vec[a], Z.vec[c]], Z.wvec)
                if o is not None:
                    O[i, j] = o + 2
        rk, reach = _rank_orders(O)
        return (rk, len(R_), reach)
    for name, r in cert["C"].items():
        Z = model_Z(blob, name)
        for s in r["supports"]:
            if s["mssm_parity"] != "viable":
                continue
            S = s["support"]
            E = exact_for(Z, S)
            lin = sorted((E.coupling([Z.vec[f["name"]]], Z.wvec) + 1, f["name"]) for f in Z.sing
                         if f["name"] not in S and E.coupling([Z.vec[f["name"]]], Z.wvec) is not None)
            new = (rank(Z, E, "d", "bd"), rank(Z, E, "bl", "l"), E.coupling([], Z.wvec), [tuple(x) for x in lin[:6]])
            old = (tuple(s["exotic_d"]), tuple(s["higgs"]), s["W_S_lowest_order"], [tuple(x) for x in s["linear_terms"]])
            checked += 1
            diffs += new != old
    for name, vs in cert["E"]["vacua"].items():
        Z = model_Z(blob, name)
        for v in vs:
            E = exact_for(Z, v["support"])
            new = (rank(Z, E, "d", "bd"), rank(Z, E, "bl", "l") if "higgs" in v else None)
            old = (tuple(v["exotic_d"]), tuple(v["higgs"]) if "higgs" in v else None)
            checked += 1
            diffs += new != old
    return dict(checked=checked, differences=diffs)


def certify_10978():
    P78 = _load("p10978", "w33_pass10978_z12_e6_admissible_r_rule.py")
    with gzip.open(P78.Z12, "rt") as fh:
        z12 = json.load(fh)
    prev = json.loads(P78.P10968.read_text())["z12_models"]
    old = json.loads((ROOT / "data" / "w33_pass10978_z12_e6_admissible_r_rule.json").read_text())
    diffs = checked = 0

    def rank(mdl, E, Xb, X):
        R_, C_ = mdl.lab(Xb), mdl.lab(X)
        if not R_:
            return [0, 0]
        B = np.array([[int(E.coupling([mdl.V[a], mdl.V[b]], mdl.W) is not None) for b in C_] for a in R_])
        return [assignment_rank(B), len(R_)]
    for name, r in sorted(prev.items()):
        if "support" not in r:
            continue
        S = r["support"]
        o = old["B"]["models"][name]
        for tag, flag in (("weak", False), ("z2", True)):
            mdl = P78.Z12Model(z12[name], flag)
            E = EXO.ExactOrders([mdl.V[s] for s in S], mdl.nq, [mdl.den])
            checked += 2
            diffs += rank(mdl, E, "d", "bd") != list(o["exotic_d_ilp"][tag])
            diffs += rank(mdl, E, "bl", "l") != list(o["higgs_ilp"][tag])
    name = "Z12-I|Z12I_1063__SM_20260926_303"
    z2 = P78.Z12Model(z12[name], True)
    S = prev[name]["support"]
    E = EXO.ExactOrders([z2.V[s] for s in S], z2.nq, [z2.den])
    lin = sorted((E.coupling([z2.V[x]], z2.W) + 1, x) for x in z2.sing()
                 if x not in S and E.coupling([z2.V[x]], z2.W) is not None)
    checked += 1
    diffs += [list(x) for x in lin[:5]] != [list(x) for x in old["C"]["linear_terms"][:5]]
    return dict(checked=checked, differences=diffs, c1063_first_linear=list(lin[0]) if lin else None)


def certify_10980():
    """rerun Pass 10980 with its ILP `order` replaced by the exact solver (same seeds) and compare the frozen summary."""
    import tempfile
    src = (ROOT / "analysis" / "w33_pass10980_rpv_flavour_upper_bound.py").read_text(encoding="utf-8")
    s, e = src.index("    def order(fields):"), src.index("    rng = np.random.default_rng(1)")
    exact_order = (
        "    EX = __EXO__.ExactOrders([vec[s] for s in S], nq, [den])\n\n"
        "    def order(fields):\n"
        "        key = tuple(sorted(fields))\n"
        "        if key not in cache:\n"
        "            o = EX.coupling([vec[f] for f in fields], [0] * dim)\n"
        "            cache[key] = None if o is None else o + len(fields)\n"
        "        return cache[key]\n\n")
    src = src[:s] + exact_order + src[e:]
    tmp = Path(tempfile.mkdtemp()) / "p10980_exact.json"
    src = src.replace('OUT = ROOT / "data" / "w33_pass10980_rpv_flavour_upper_bound.json"', "OUT = __OUT__")
    ns = {"__EXO__": EXO, "__OUT__": tmp, "__name__": "p10980_exact"}
    import contextlib, io
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(src, "p10980_exact", "exec"), ns)
    new = json.loads(tmp.read_text())["summary"]
    old = json.loads((ROOT / "data" / "w33_pass10980_rpv_flavour_upper_bound.json").read_text())["summary"]
    return dict(identical_summary=new == old, summary=new)



def _centre(tok, N):
    """centre charge (mod 1) of a representation token of SU(N): fundamental 1/N, antifundamental -1/N,
    SU(4) 6 -> 1/2, SU(2) doublet -> 1/2, singlet 0 (other reps: 0, flagged by the caller)."""
    d = dimof(tok)
    if d == 1:
        return Fraction(0)
    if N == 2 and d == 2:
        return Fraction(1, 2)
    if d == N:
        return Fraction(-1, N) if tok.startswith("-") else Fraction(1, N)
    if N == 4 and d == 6:
        return Fraction(1, 2)
    return None


def _components(tok, N, broken):
    """SU(N-1) (or, for SU(2), trivial) component tokens of a representation of a broken SU(N) factor."""
    d = dimof(tok)
    if not broken or d == 1:
        return [tok]
    if N == 2:
        return ["1", "1"]
    neg = tok.startswith("-")
    if d == N:
        return [("-" if neg else "") + str(N - 1), "1"]
    if N == 4 and d == 6:
        return ["3", "-3"]
    return None


def unbroken_group(rows_all, rows_vac):
    """L_all / L_vac for integer charge lattices (rows): (number of free directions = unbroken U(1)s,
    torsion invariant factors = the unbroken discrete symmetry of the vacuum)."""
    from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
    H = hermite_normal_form(sp.Matrix(rows_all).T).T
    basis = [list(H.row(i)) for i in range(H.rows) if any(H.row(i))]
    Bm = sp.Matrix(basis)
    coords = []
    for v in rows_vac:
        sol = Bm.T.gauss_jordan_solve(sp.Matrix(v))[0]
        sol = sol.subs({s_: 0 for s_ in sol.free_symbols})
        assert all(x.is_integer for x in sol)
        coords.append([int(x) for x in sol])
    S_ = smith_normal_form(sp.Matrix(coords), domain=sp.ZZ)
    inv = [abs(int(S_[i, i])) for i in range(min(S_.shape)) if S_[i, i] != 0]
    return len(basis) - len(inv), [d for d in inv if d > 1]


def condensate_candidate(blob, name, support):
    """A viable composite vacuum analysed with its hidden factor(s) broken (minimal aligned breaking:
    SU(N) -> SU(N-1) by 4.4bar-type mesons, SU(2) fully by doublet pairs):
      (1) parity with the broken Cartan aligned across all condensates of one factor;
      (2) every SM-charged class after breaking (components expanded), mass ranks from exact orders on the vacuum
          fields with centre charges of every hidden SU(N) as extra discrete charges (N-ality; this can only
          over-count couplings, so 'light' is conservative);
      (3) W on the vacuum and linear terms."""
    Z = model_Z(blob, name)
    T = P79.Tools(Z)
    dims = {n: Z.raw[n]["dim"].split(",") for n in Z.raw}
    _, qd, c3, w2 = slots(blob["ledger"][name])
    hid = [i for i in range(len(qd)) if i not in (c3, w2)]
    present = defaultdict(set)
    for d in dims.values():
        for i, x in enumerate(d):
            present[i].add(x)
    Nof = {}
    for i in hid:
        cplx = sorted(dimof(v) for v in present[i] if not v.startswith("-") and dimof(v) > 1 and "-" + v in present[i])
        Nof[i] = cplx[0] if cplx else (2 if "2" in present[i] else None)
    comps = [f for f in support if "of" in f]
    if any(f["kind"] not in ("su2pair", "meson") for f in comps):
        return dict(analysed=False, reason="only SU(2)-pair and SU(N)-meson condensates are analysed")
    singles = [f["name"] for f in support if "of" not in f]
    slot_of = lambda n: [i for i in hid if dimof(dims[n][i]) != 1][0]
    broken = sorted({slot_of(f["of"][0]) for f in comps})
    realw = lambda ev, od: M.realizable([list(t) for t in od] + [list(t) for t in ev] + [list(t) for t in Z.virtual]
                                        + [Z.wvec], list(range(len(od))))
    base = [Z.vec[s] for s in singles] + [[a + b for a, b in zip(Z.vec[f["of"][0]], Z.vec[f["of"][1]])] for f in comps]
    par_composite = P8.match_search(realw, base, T.species, T.leptons)[0]
    aligned = par_composite
    by_slot = defaultdict(list)
    for f in comps:
        by_slot[slot_of(f["of"][0])].append(f)
    if par_composite == "viable" and any(len(v) > 1 for v in by_slot.values()):
        aligned = "none"
        opts = []
        for s_, fs in by_slot.items():
            if len(fs) < 2:
                continue
            first = fs[0]["of"]
            # meson: the fundamental constituent is the one whose token is positive; SU(2) pairs may swap
            fund = lambda pr: pr[0] if not dims[pr[0]][s_].startswith("-") else pr[1]
            for f in fs[1:]:
                if f["kind"] == "meson":
                    opts.append([[x - z for x, z in zip(Z.vec[fund(first)], Z.vec[fund(f["of"])])]])
                else:
                    opts.append([[x - z for x, z in zip(Z.vec[first[0]], Z.vec[f["of"][k]])] for k in (0, 1)])
        for choice in itertools.product(*[range(len(o)) for o in opts]):
            extra = [o[c] for o, c in zip(opts, choice)]
            if P8.match_search(realw, base + extra, T.species, T.leptons)[0] == "viable":
                aligned = "viable"
                break
    S = singles + [x for f in comps for x in f["of"]]
    cen_slots = [i for i in hid if Nof[i]]

    def vx(n):
        return list(Z.vec[n]) + [(_centre(dims[n][i], Nof[i]) or Fraction(0)) for i in cen_slots]
    wv = list(Z.wvec) + [Fraction(0)] * len(cen_slots)
    E = EXO.ExactOrders([vx(s) for s in S], Z.nq,
                        [x.denominator for v in list(Z.vec.values()) + [Z.wvec] + Z.virtual for x in v] + [2, 3, 4])
    lab = [f for f in Z.left if f["base"] in P.SMY]
    y = P.solve_affine([f["q"] for f in lab], [P.SMY[f["base"]] for f in lab])[0]
    Yv = {f["name"]: P.dot(y, f["q"]) for f in Z.left}
    cls = defaultdict(list)
    undecomposed = 0
    for f in Z.left:
        n, d = f["name"], dims[f["name"]]
        if dimof(d[c3]) == 1 and dimof(d[w2]) == 1 and Yv[n] == 0:
            continue
        parts = [[]]
        ok = True
        for i in hid:
            cs = _components(d[i], Nof[i], i in broken)
            if cs is None:
                ok = False
                break
            parts = [p + [c] for p in parts for c in cs]
        if not ok:
            undecomposed += 1
            continue
        for p in parts:
            cls[(d[c3], d[w2], Yv[n], tuple(p))].append(n)

    def conj(k):
        hp = []
        for i, x in zip(hid, k[3]):
            if i in broken:
                hp.append(x if (Nof[i] == 2 or dimof(x) == 1) else (x[1:] if x.startswith("-") else "-" + x))
            else:
                hp.append(conj_token(x, present[i], Nof[i] == 2))
        return (conj_token(k[0], present[c3], False), k[1], -k[2], tuple(hp))
    rows, seen = [], set()
    light_frac = light_frac_sym = 0
    cq = qd[c3]
    for K in sorted(cls, key=str):
        Kb = conj(K)
        if K in seen or Kb in seen:
            continue
        seen.update([K, Kb])
        A_, B_ = cls[K], cls.get(Kb, [])
        rk = rl = 0
        if B_:
            rk = assignment_rank(np.array([[int(E.coupling([vx(a), vx(b)], wv) is not None) for b in B_] for a in A_]))
            rl = assignment_rank(np.array([[int(E.coupling_lattice([vx(a), vx(b)], wv)) for b in B_] for a in A_]))
        la, lb = len(A_) - rk, len(B_) - rk
        wk = dimof(K[1])
        charges = [K[2] + Fraction(t_, 2) for t_ in range(-(wk - 1), wk, 2)]
        col = dimof(K[0])
        frac = (col == 1 and any(c.denominator != 1 for c in charges)) or \
               (col == 3 and any((c - (Fraction(2, 3) if K[0] == cq else Fraction(1, 3))).denominator != 1 for c in charges))
        rows.append(dict(cls=str(K), n=len(A_), nbar=len(B_), rank=rk, rank_ignoring_holomorphy=rl, light=[la, lb],
                         fractional=frac))
        light_frac += (la + lb) if frac else 0
        light_frac_sym += (len(A_) + len(B_) - 2 * rl) if frac else 0
    den_ = lcm(*[x.denominator for v in list(Z.vec.values()) + [Z.wvec] + Z.virtual for x in v], 12)
    ncen = len(cen_slots)
    row_ = lambda n: [int(x * den_) for x in vx(n)]
    per = [[int(x * den_) for x in v] + [0] * ncen for v in Z.virtual] + \
          [[0] * len(Z.wvec) + [den_ if j == k else 0 for j in range(ncen)] for k in range(ncen)]
    free_u1, torsion = unbroken_group([row_(f["name"]) for f in Z.left] + per, [row_(s) for s in S] + per)
    lin = sorted((E.coupling([vx(f["name"])], wv) + 1, f["name"]) for f in Z.sing
                 if f["name"] not in S and E.coupling([vx(f["name"])], wv) is not None)
    return dict(analysed=True, broken_slots=broken, parity_composites_even=par_composite, parity_aligned=aligned,
                classes=rows, light_fractional_multiplets=light_frac,
                light_fractional_protected_by_unbroken_symmetry=light_frac_sym,
                unbroken_u1s=free_u1, unbroken_discrete_invariant_factors=torsion, undecomposed_fields=undecomposed,
                W_on_vacuum=E.coupling([], wv), linear_terms=lin[:8])


def certify_fast_parity(blob, cert):
    """FastParity vs Pass 10967's realizable() (via Z.R): every certified FI support (known positives) plus 60 seeded
    random singlet sets per parity model."""
    import random
    tot = agree = 0
    for name in ("z3z3|Z3Z3_0001_c3__SM_20260934_739", "z3z3|Z3Z3_0001_c4__SM_20260935_3729",
                 "z3z3|Z3Z3_0001_c1__SM_20260932_2822"):
        Z = model_Z(blob, name)
        F = FastParity(Z.nq, Z.d2["nonR_orders"], Z.matter, Z.virtual + [Z.wvec])
        sing = [f["name"] for f in Z.sing]
        rng = random.Random(7)
        tests = [s["support"] for s in cert["C"].get(name, {}).get("supports", [])]
        tests += [rng.sample(sing, rng.randint(1, 6)) for _ in range(60)]
        for S in tests:
            ev = [Z.vec[s] for s in S]
            tot += 1
            agree += Z.R(ev) == F.realizable(ev)
    return dict(tests=tot, agree=agree)


def main():
    blob = P79.load()
    cert = json.loads(CERT79.read_text())
    A, B = {}, {}
    Asum, Bsum = Counter(), Counter()
    for name, vs in cert["E"]["vacua"].items():
        Z = model_Z(blob, name)
        for v in vs:
            if "dflat_full" not in v:
                continue
            st = structure(Z, v["support"])
            A.setdefault(name, []).append(dict(support=v["support"], **st))
            Asum["vacua"] += 1
            Asum["z3_graded_all_orders"] += st["z3_graded"] and st["periodic"] and st["lattice_index"] == 1
            Asum["k=%d" % st["k"]] += 1
            Asum["disjoint_composites"] += st["disjoint_supports"]
            # all exponent vectors live in the k-dim span of the rays, so the first dependency is the first
            # monomial beyond the k linear ones: the lowest |a| = 4 term, degree 4 * min(ray degree)
            st["predicted_first_dependent_degree"] = 4 * min(st["ray_degrees"])
            Asum["prediction_matches_10979"] += (v["first_dependent_degree"] == st["predicted_first_dependent_degree"]
                                                 or (v["first_dependent_degree"] is None and st["predicted_first_dependent_degree"] > P79.DMAX))
            r = susy_roots(Z, v["support"], st)
            B.setdefault(name, []).append(r)
            Bsum["vacua"] += 1
            Bsum["with_verified_susy_roots"] += r["verified"] > 0 and r["verified"] == r["roots"]
            print("AB", name.split("|")[1][-9:], st["k"], st["ray_degrees"], st["z3_graded"], r["verified"], "/", r["roots"],
                  r["max_abs_median"], flush=True)
    maxes = [r["max_abs_median"] for rs in B.values() for r in rs if r["max_abs_median"]]
    Bsum["max_abs_median_range"] = [min(maxes), max(maxes)]
    print("A", dict(Asum), "\nB", dict(Bsum), flush=True)
    # C: index for all 13 Z3xZ3 models, hidden beta functions, controls
    C = {}
    for name in sorted(blob["ledger"]):
        ix = chiral_index(blob["ledger"][name])
        C[name] = dict(index=ix, beta=hidden_beta(blob["ledger"][name]))
        print("C", name.split("|")[1], ix["fractional_states"], ix["unconfinable_fractional"], C[name]["beta"], flush=True)
    ctrl = {}
    L6, _ = P.load_ledger()
    fams = {"Z6-I/II": L6}
    for fam, path in CONTROLS.items():
        with gzip.open(path, "rt") as fh:
            fams[fam] = json.load(fh)["ledger"]
    for fam, LL in fams.items():
        cs = Counter()
        for n, m in LL.items():
            try:
                ix = chiral_index(m)
            except (IndexError, StopIteration, TypeError, ValueError, AssertionError):
                cs["slot_identification_failed"] += 1
                continue
            cs["models"] += 1
            cs["zero_fractional_excess"] += ix["fractional_states"] == 0
            cs["zero_unconfinable"] += ix["unconfinable_fractional"] == 0
        ctrl[fam] = dict(cs)
        print("C control", fam, dict(cs), flush=True)
    # D: hidden condensates
    D = {}
    for name in sorted(blob["ledger"]):
        b = name.split("|")[1]
        D[name] = condensate_census(name, blob["ledger"][name], blob["disc"][b])
        if D[name]["census"].get("realizable_rays_needing_a_composite"):
            D[name]["composite_vacua"] = composite_vacua_parity(blob, name, D[name])
            D[name]["candidates"] = [condensate_candidate(blob, name, s_["fields"])
                                     for s_ in D[name]["composite_vacua"]["supports"] if s_["mssm_parity"] == "viable"]
        vecs_ = D[name].pop("_vecs", None)
        print("D", b, D[name]["composites"], D[name]["census"], D[name].get("composite_vacua", {}).get("summary"), flush=True)
    rescue = {}
    for name, vs in cert["E"]["vacua"].items():
        rescue[name] = su4_rescue(blob, name, [v for v in vs if "dflat_full" in v])
        print("D rescue", name.split("|")[1], rescue[name]["tally"], flush=True)
    E = dict(p10979=certify_10979(blob, cert), p10978=certify_10978(), p10980=certify_10980(),
             fast_parity=certify_fast_parity(blob, cert))
    print("E", E, flush=True)
    OUT.write_text(json.dumps(dict(pass_id=11024, A=dict(summary=dict(Asum), vacua=A), B=dict(summary=dict(Bsum), vacua=B),
                                   C=dict(models=C, controls=ctrl), D=dict(census=D, su4_rescue=rescue), E=E),
                              indent=1, sort_keys=True, default=str))


if __name__ == "__main__":
    from multiprocessing import Pool
    with Pool(8) as _p:
        _POOL = _p
        main()
