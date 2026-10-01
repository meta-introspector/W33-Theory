#!/usr/bin/env python3
"""Pass 11241: mu from the FULL vacuum symmetry -- continuous U(1)s, their discrete remnants, space group and R.

Pass 11232 tested only the continuous U(1)s left unbroken by an FI-cancelling singlet vacuum: in the Z6-I class (the
W33 vacuum's) none can forbid mu (0/23 models).  The scorecard sentence "mu must come from a vacuum symmetry" can still
hold through a DISCRETE symmetry.  The full symmetry of a vacuum S is decided exactly by a lattice:
  a superpotential term with fields F is allowed at some order iff   sum_{f in F} v_f - w  lies in
  Lambda_S = Z-span{ v_s : s in S }  +  Z-span{ integrality vectors of the discrete charges },
where v_f = (U(1) charges, discrete charges / order, R charges / order) is the frozen field vector of Pass 10967
(data/w33_pass10967_space_group_discrete_charges.json.gz; for Z6-II the geometric R-symmetries Z6^R x Z3^R x Z2^R with
W charge -1 each, as used in Pass 10968; for Z6-I the point-group Z6 only -- no R charges are recorded) and w is the
superpotential's R vector.  Integer (not rational) spans capture the discrete remnants of broken U(1)s (a U(1) broken by
a charge-2 field leaves Z2).  Positivity of exponents is ignored, so "forbidden" verdicts are exact (no order allows the
term), "allowed" verdicts are necessary conditions.

Per vacuum (an FI-cancelling extreme ray of the D-flat cone -- Pass 11232 -- with a choice, for each U(1) type in its
support, of which field of that type condenses):
  * mu-protected H_u = bl_i: every bl_i l_j forbidden, some q_a bl_i bu_b allowed;
  * clean: other bl_k all paired (generic matching rank), all vector-like pairs paired, and a mu-generating s' exists;
  * F-flat by lattice: no singlet o outside S with v_o - w in Lambda (F_o = 0 at all orders) and -w not in Lambda
    (W restricted to the vacuum vanishes identically) -- sufficient for F-flatness to all orders.
Z6-I is enumerated exhaustively (every FI ray, every field choice); Z6-II is capped (rays and choices sampled; stated).
"""
from __future__ import annotations

import gzip
import importlib.util
import itertools
import json
import random
import sys
import time
from fractions import Fraction
from math import gcd, lcm
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass10960_matter_even_dflat_closure as P  # noqa: E402
import w33_pass11232_mu_vacuum_symmetry as MU  # noqa: E402

_spec = importlib.util.spec_from_file_location("p10967", ROOT / "analysis" / "w33_pass10967_fi_vacuum_regenerates_rpv.py")
M67 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(M67)

DISC = ROOT / "data" / "w33_pass10967_space_group_discrete_charges.json.gz"
OUT = ROOT / "data" / "w33_pass11241_discrete_mu_symmetry.json"


# ----------------------------------------------------------------------------- exact integer lattices
class Lattice:
    """Z-span of integer vectors in row-echelon (Hermite-like) form; exact membership by back-substitution"""

    def __init__(self, rows):
        rows = [list(r) for r in rows if any(r)]
        self.ech = []
        dim = len(rows[0]) if rows else 0
        col = 0
        while rows and col < dim:
            nz = [r for r in rows if r[col] != 0]
            if not nz:
                col += 1
                continue
            while len(nz) > 1:
                nz.sort(key=lambda r: abs(r[col]))
                p = nz[0]
                new = [p]
                for r in nz[1:]:
                    k = r[col] // p[col]
                    rr = [a - k * b for a, b in zip(r, p)]
                    if rr[col] != 0:
                        new.append(rr)
                    elif any(rr):
                        rows.append(rr)
                rows = [r for r in rows if r[col] == 0 and any(r)]
                nz = new
            p = nz[0]
            if p[col] < 0:
                p = [-a for a in p]
            self.ech.append((col, p))
            rows = [r for r in rows if r[col] == 0 and any(r)]
            col += 1

    def __contains__(self, v):
        v = list(v)
        for c, p in self.ech:
            if v[c] % p[c]:
                return False
            k = v[c] // p[c]
            if k:
                v = [a - k * b for a, b in zip(v, p)]
        return not any(v)


# ----------------------------------------------------------------------------- one model
class Model:
    def __init__(self, name, model, disc):
        self.name = name
        d = disc[name.split("|")[1]]
        RN = [(int(x.split(":")[0]), Fraction(x.split(":")[1])) for x in d["R"]]
        d2 = json.loads(json.dumps(d))
        d2["nonR_orders"] = d["nonR_orders"] + [n for n, _ in RN]
        for v in d2["fields"].values():
            v["nonR"] = v["nonR"] + v["R"]
        self.orders, self.left = M67.fields(model, d2)
        self.nq = len(self.left[0]["q"])
        nd = len(self.orders)
        self.R_orders = [n for n, _ in RN]
        w = [Fraction(0)] * (self.nq + len(d["nonR_orders"])) + [x / n for n, x in RN]
        allv = [f["q"] + f["ext"] for f in self.left] + [w]
        self.den = lcm(*[x.denominator for v in allv for x in v])
        I = lambda v: tuple(int(x * self.den) for x in v)
        self.vec = {f["name"]: I(f["q"] + f["ext"]) for f in self.left}
        self.w = I(w)
        self.virtual = [tuple(self.den * int(i == j + self.nq) for i in range(self.nq + nd)) for j in range(nd)]
        lab = [f for f in self.left if f["base"] in P.SMY]
        y = P.solve_affine([f["q"] for f in lab], [P.SMY[f["base"]] for f in lab])[0]
        self.sing = [f for f in self.left if f["trivial"] and P.dot(y, f["q"]) == 0]
        self.byb = {}
        for f in self.left:
            self.byb.setdefault(f["base"], []).append(f["name"])
        self.tr0 = sum((f["dimprod"] * f["q"][0] for f in self.left), Fraction(0))

    def u1_types(self):
        types = {}
        for f in self.sing:
            types.setdefault(tuple(f["q"]), []).append(f["name"])
        # distinct full vectors only
        return {t: list(dict.fromkeys(sorted(ns, key=lambda n: self.vec[n]))) for t, ns in types.items()}

    def target(self, names):
        return tuple(sum(self.vec[n][k] for n in names) - self.w[k] for k in range(len(self.w)))

    def analyse_vacuum(self, S):
        L = Lattice([self.vec[s] for s in S] + self.virtual)
        hu, hd = self.byb.get("bl", []), self.byb.get("l", [])
        qs, us = self.byb.get("q", []), self.byb.get("bu", [])
        prot = []
        for h in hu:
            if any(self.target([h, l]) in L for l in hd):
                continue
            if not any(self.target([q, h, u]) in L for q in qs for u in us):
                continue
            prot.append(h)
        if not prot:
            return None
        rec = dict(support=list(S), protected=prot)
        clean = []
        for h in prot:
            rest = [k for k in hu if k != h]
            r = MU.structural_rank(rest, hd, lambda a, b: self.target([a, b]) in L)
            gen = [o["name"] for o in self.sing if o["name"] not in S
                   and any(self.target([o["name"], h, l]) in L for l in hd)]
            ex = {}
            for b in ("q", "u", "d", "e", "v", "x"):
                A, B = self.byb.get(b, []), self.byb.get("b" + b, [])
                if A and B:
                    k = MU.structural_rank(A, B, lambda a, c: self.target([a, c]) in L)
                    ex[b] = min(len(A), len(B)) - k
            ok = r == len(rest) and gen and all(v == 0 for v in ex.values())
            clean.append(dict(Hu=h, other_Hu_all_massive=r == len(rest), generators=len(gen),
                              light_vectorlike=ex, clean=bool(ok)))
        rec["detail"] = clean
        rec["clean"] = any(c["clean"] for c in clean)
        # lattice F-flatness (sufficient, all orders)
        minus_w = tuple(-x for x in self.w)
        lin = [o["name"] for o in self.sing if o["name"] not in S and self.target([o["name"]]) in L]
        rec["W_on_vacuum_possible"] = minus_w in L
        rec["F_linear_outside"] = len(lin)
        rec["F_flat_lattice"] = (minus_w not in L) and not lin
        return rec


def model_rays(Mo, max_types_exact=44):
    types = Mo.u1_types()
    T = list(types)
    s = 1 if Mo.tr0 > 0 else -1
    C = [s * t[0] for t in T]
    Mx = [list(t[1:]) for t in T]
    if P.farkas(C, Mx) is not None:
        return T, types, None
    if len(T) <= max_types_exact:
        rays = [r for r in P.cone_rays(Mx) if P.dot(C, r) < 0]
    else:
        rays = MU.sampled_vertices(C, Mx)
    return T, types, rays


def analyse_model(name, model, disc, ray_cap=None, choice_cap=None, seed=11241):
    rng = random.Random(seed)
    Mo = Model(name, model, disc)
    if Mo.tr0 == 0:
        return dict(anomalous=False)
    T, types, rays = model_rays(Mo)
    if rays is None:
        return dict(anomalous=True, dflat=False)
    rec = dict(anomalous=True, dflat=True, fi_rays=len(rays), discrete_orders=Mo.orders, R_orders=Mo.R_orders)
    supports = {tuple(k for k, v in enumerate(r) if v != 0) for r in rays}
    supports = sorted(supports)
    exhaustive = True
    if ray_cap and len(supports) > ray_cap:
        supports = rng.sample(supports, ray_cap)
        exhaustive = False
    vac = prot = clean = cleanF = 0
    examples, clean_examples = [], []
    for sup in supports:
        choices = [types[T[k]] for k in sup]
        combos = itertools.product(*choices)
        n = 1
        for c in choices:
            n *= len(c)
        if choice_cap and n > choice_cap:
            combos = (tuple(rng.choice(c) for c in choices) for _ in range(choice_cap))
            exhaustive = False
        for S in combos:
            vac += 1
            r = Mo.analyse_vacuum(S)
            if r is None:
                continue
            prot += 1
            if len(examples) < 3:
                examples.append(r)
            if r["clean"]:
                clean += 1
                cleanF += r["F_flat_lattice"]
                if len(clean_examples) < 5:
                    clean_examples.append(r)
    rec.update(vacua=vac, exhaustive=exhaustive, mu_protected=prot, clean=clean, clean_and_F_flat=cleanF,
               examples=examples, clean_examples=clean_examples)
    return rec


def main():
    ledger, sha = P.load_ledger()
    disc = M67.load(DISC)
    names = sorted(ledger)
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if args:
        names = [n for n in names if any(a in n for a in args)]
    res, t0 = {}, time.time()
    for n in names:
        t = time.time()
        z6i = n.startswith("Z6-I|")
        res[n] = analyse_model(n, ledger[n], disc, ray_cap=None if z6i else 150, choice_cap=None if z6i else 4)
        r = res[n]
        print(n, {k: v for k, v in r.items() if k not in ("examples", "clean_examples")}, f"{time.time() - t:.1f}s",
              flush=True)
    dflat = [k for k, v in res.items() if v.get("dflat")]
    z6i = [k for k in dflat if k.startswith("Z6-I|")]
    summary = dict(models=len(res), dflat=len(dflat),
                   z6i_dflat=len(z6i), z6i_exhaustive=all(res[k]["exhaustive"] for k in z6i),
                   z6i_models_mu_protected=sorted(k for k in z6i if res[k]["mu_protected"]),
                   z6i_models_clean=sorted(k for k in z6i if res[k]["clean"]),
                   models_mu_protected=sorted(k for k in dflat if res[k]["mu_protected"]),
                   models_clean=sorted(k for k in dflat if res[k]["clean"]),
                   models_clean_and_F_flat=sorted(k for k in dflat if res[k]["clean_and_F_flat"]),
                   ledger_sha256=sha, seconds=round(time.time() - t0, 1))
    if not args:
        OUT.write_text(json.dumps(dict(pass_id=11241, summary=summary, models=res), indent=1, default=str))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
