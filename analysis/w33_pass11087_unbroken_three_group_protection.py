#!/usr/bin/env python3
"""Pass 11087: what symmetry protects the massless fractional charges of the Z3xZ3 vacua?

Pass 11024 found that in every Z3xZ3 parity vacuum with the hidden SU(4) condensed, at least 120 fractionally
charged multiplets stay massless at all orders, most of them because an exact symmetry of the vacuum forbids their
mass terms (no integer exponents of any sign exist) -- an unbroken discrete 3-group of order 3^6 or 3^4.  Here
the protecting symmetry is taken apart.  The charge coordinates fall into blocks:
    U1     the continuous U(1) charges (the VEVs break them to discrete remnants),
    PG     the two point-group (sector) charges theta^k omega^l,
    SG     the three space-group translation charges (the fixed-point labels, one Z3 per torus),
    R      the three Z3^R plane-rotation charges,
    centre the centre (N-ality) charges of the hidden SU(4), SU(3)'.
For every vacuum:
  (a) the unbroken group G = torsion(L_all / L_vac) with all blocks, and the order it has when one block is dropped
      (a block whose removal shrinks |G| carries part of the unbroken symmetry);
  (b) for every lattice-forbidden mass term in the hidden-singlet charge-1/3 class, which blocks are NECESSARY
      (dropping the block makes the term allowed) -- or whether the protection is redundant (no single block is
      necessary: independent symmetries forbid it separately).
Vacua: the 55 hidden-condensate parity vacua of Pass 11024 (c1 1, c3 39, c4 15) and its 28 singlet vacua.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from math import lcm
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11024_z3xz3_susy_vacua_fractional_charge_index as P  # noqa: E402

CERT = ROOT / "data" / "w33_pass11024_z3xz3_susy_vacua_fractional_charge_index.json"
CERT79 = ROOT / "data" / "w33_pass10979_z3xz3_parity_vacua_f_obstruction.json"
OUT = ROOT / "data" / "w33_pass11087_unbroken_three_group_protection.json"
BLOB = None


def _blob():
    global BLOB
    if BLOB is None:
        BLOB = P.P79.load()
    return BLOB


def setup(name, S):
    blob = _blob()
    Z = P.model_Z(blob, name)
    dims = {n: Z.raw[n]["dim"].split(",") for n in Z.raw}
    _, qd, c3, w2 = P.slots(blob["ledger"][name])
    hid = [i for i in range(len(qd)) if i not in (c3, w2)]
    present = defaultdict(set)
    for dd in dims.values():
        for i, x in enumerate(dd):
            present[i].add(x)
    Nof = {}
    for i in hid:
        cp = sorted(P.dimof(v) for v in present[i] if not v.startswith("-") and P.dimof(v) > 1 and "-" + v in present[i])
        Nof[i] = cp[0] if cp else (2 if "2" in present[i] else None)
    cen = [i for i in hid if Nof[i]]
    vx = lambda n: list(Z.vec[n]) + [(P._centre(dims[n][i], Nof[i]) or Fraction(0)) for i in cen]
    wv = list(Z.wvec) + [Fraction(0)] * len(cen)
    nq = Z.nq
    nnr = len(Z.d2["nonR_orders"]) - len(Z.RN)
    blocks = {"U1": list(range(nq)), "PG": [nq, nq + 1], "SG": list(range(nq + 2, nq + nnr)),
              "R": list(range(nq + nnr, len(Z.wvec))), "centre": list(range(len(Z.wvec), len(Z.wvec) + len(cen)))}
    # periodicities: every discrete coordinate is defined mod 1 (virtual unit vectors), centre coords mod 1
    per = [[Fraction(int(i == j)) for i in range(len(wv))] for j in range(nq, len(wv))]
    return Z, dims, vx, wv, blocks, per


def group_order(Z, vx, S, per, keep):
    den = lcm(*[x.denominator for v in list(Z.vec.values()) + [Z.wvec] for x in v], 36)
    row = lambda v: [int(v[i] * den) for i in keep]
    pk = [row(p) for p in per if any(p[i] for i in keep)]
    free, tors = P.unbroken_group([row(vx(f["name"])) for f in Z.left] + pk, [row(vx(s)) for s in S] + pk)
    o = 1
    for d in tors:
        o *= d
    return free, tors, o


def analyse(job):
    name, S, tag = job
    Z, dims, vx, wv, blocks, per = setup(name, S)
    allc = sum(blocks.values(), [])
    free, tors, order = group_order(Z, vx, S, per, allc)
    drop = {}
    for k, cs in blocks.items():
        if cs:
            f2, t2, o2 = group_order(Z, vx, S, per, [i for i in allc if i not in cs])
            drop[k] = dict(order=o2, invariant_factors=t2)
    cache = {}

    def forbidden(pair, keep, nq_keep):
        key = tuple(keep)
        if key not in cache:
            cache[key] = P.EXO.ExactOrders([[vx(s)[i] for i in keep] for s in S], nq_keep, [2, 3, 4, 9])
        E = cache[key]
        return not E.coupling_lattice([[vx(f)[i] for i in keep] for f in pair], [wv[i] for i in keep])
    lab = [f for f in Z.left if f["base"] in P.P.SMY]
    y = P.P.solve_affine([f["q"] for f in lab], [P.P.SMY[f["base"]] for f in lab])[0]
    Yv = {f["name"]: P.P.dot(y, f["q"]) for f in Z.left}
    single = lambda n: all(P.dimof(t) == 1 for t in dims[n])
    X = [f["name"] for f in Z.left if single(f["name"]) and Yv[f["name"]] == Fraction(-1, 3)]
    Xb = [f["name"] for f in Z.left if single(f["name"]) and Yv[f["name"]] == Fraction(1, 3)]
    nec = Counter()
    for a in X:
        for b in Xb:
            if not forbidden([a, b], allc, Z.nq):
                nec["allowed"] += 1
                continue
            need = [k for k, cs in blocks.items() if cs and
                    not forbidden([a, b], [i for i in allc if i not in cs], 0 if k == "U1" else Z.nq)]
            nec["forbidden"] += 1
            nec["necessary: " + "+".join(need) if need else "redundant (no single block necessary)"] += 1
    return dict(model=name, tag=tag, support=S, unbroken_u1s=free, invariant_factors=tors, order=order,
                order_with_block_dropped=drop, charge_third_pairs=dict(nec))


def jobs():
    cert = json.loads(CERT.read_text())
    cert79 = json.loads(CERT79.read_text())
    out = []
    for name, D in cert["D"]["census"].items():
        if "composite_vacua" not in D:
            continue
        for s in D["composite_vacua"]["supports"]:
            if s["mssm_parity"] != "viable":
                continue
            f = s["fields"]
            S = [x["name"] for x in f if "of" not in x] + [c for x in f if "of" in x for c in x["of"]]
            out.append((name, S, "condensate"))
    for name, vs in cert79["E"]["vacua"].items():
        for v in vs:
            if "dflat_full" in v:
                out.append((name, v["support"], "singlet"))
    return out


def extra_u1_couplings():
    """In every singlet vacuum: an unbroken U(1)' beyond hypercharge (t.q(s) = 0 on the vacuum, t independent of Y);
    does Standard-Model matter carry charge under it (a massless gauge boson coupled to quarks and leptons)?"""
    import sympy as sp
    cert79 = json.loads(CERT79.read_text())
    out = Counter()
    for name, vs in cert79["E"]["vacua"].items():
        Z = P.model_Z(_blob(), name)
        lab = [f for f in Z.left if f["base"] in P.P.SMY]
        y = P.P.solve_affine([f["q"] for f in lab], [P.P.SMY[f["base"]] for f in lab])[0]
        yv = sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in y])
        for v in vs:
            if "dflat_full" not in v:
                continue
            Qm = sp.Matrix([[sp.Rational(Z.vec[s][k].numerator, Z.vec[s][k].denominator) for k in range(Z.nq)]
                            for s in v["support"]])
            extra = [k for k in Qm.nullspace() if sp.Matrix.hstack(yv, k).rank() == 2]
            if not extra:
                out["no extra U(1)"] += 1
                continue
            t = extra[0]
            ch = lambda n: sum(sp.Rational(Z.vec[n][k].numerator, Z.vec[n][k].denominator) * t[k] for k in range(Z.nq))
            sm = [f["name"] for f in Z.left if f["base"] in ("q", "bu", "bd", "l", "be", "bl")]
            out["extra U(1)' coupling to SM matter" if any(ch(n) != 0 for n in sm) else "extra U(1)', SM neutral"] += 1
    return dict(out)


def main():
    from multiprocessing import Pool
    if "--zprime-only" in sys.argv and OUT.exists():
        d = json.loads(OUT.read_text())
        d["singlet_vacua_extra_u1"] = extra_u1_couplings()
        print(d["singlet_vacua_extra_u1"])
        OUT.write_text(json.dumps(d, indent=1, sort_keys=True, default=str))
        return
    js = jobs()
    with Pool(8) as pool:
        res = pool.map(analyse, js, chunksize=1)
    summ = {}
    for tag in ("condensate", "singlet"):
        rs = [r for r in res if r["tag"] == tag]
        tot = Counter()
        for r in rs:
            tot.update(r["charge_third_pairs"])
        summ[tag] = dict(vacua=len(rs), orders=dict(Counter(str(r["order"]) for r in rs)),
                         groups=dict(Counter(str(r["invariant_factors"]) for r in rs)),
                         unbroken_u1s=dict(Counter(r["unbroken_u1s"] for r in rs)),
                         order_drop=dict(Counter(tuple(sorted((k, v["order"]) for k, v in r["order_with_block_dropped"].items()))
                                                  .__str__() for r in rs)),
                         charge_third_pairs=dict(tot))
        print(tag, json.dumps(summ[tag], indent=1)[:2500], flush=True)
    OUT.write_text(json.dumps(dict(pass_id=11087, summary=summ, vacua=res, singlet_vacua_extra_u1=extra_u1_couplings()),
                              indent=1, sort_keys=True, default=str))


if __name__ == "__main__":
    main()
