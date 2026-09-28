#!/usr/bin/env python3
"""Pass 11088: baryonic and mixed hidden condensates in Z3xZ3 -- the named residue of Pass 11024.

Pass 11024 searched every FI-cancelling direction containing a QUADRATIC condensate of the curing factor (the hidden
factor whose breaking alone makes the fractional-charge index vanish: SU(4) in c3/c4, SU(2)_1 in c1).  Baryonic
invariants (epsilon contractions of N distinct fundamentals; for c3 430 antibaryons + 155 baryons) were left out: they
made the discrete-choice search explode.  With the FastParity engine (Pass 11024 E, validated 217/217) the full census
is now run:
  * composites: every single-factor hidden invariant (mesons, baryons, antibaryons, SU(2) pairs, 6.6);
  * every FI-cancelling extreme ray containing a condensate of the curing factor is searched for a realizable
    parity -- exactly, by a subset-sum dynamic programme over the finite group of residue vectors
    (`realizable_choice`, validated against the depth-first search 99/99);
  * each realizable ray with an MSSM-viable parity (composite-even: a NECESSARY condition when a factor is broken
    by several condensates) is analysed with the curing factor broken:
      - aligned breaking SU(N) -> SU(N-1) if only mesons of that factor condense,
      - FULL breaking (every multiplet splits into singlets) if any baryon/antibaryon/6.6 of that factor condenses;
    centre (N-ality) charges enter as extra discrete charges, and mass ranks come from exact orders.
Both simplifications (composite-even parity, N-ality-allowed couplings with full breaking) can only OVER-state
viability and masses, so a vacuum found to keep massless fractional charges is excluded without further assumptions.
"""
from __future__ import annotations

import itertools
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from math import lcm
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11024_z3xz3_susy_vacua_fractional_charge_index as P  # noqa: E402
from w33_fast_parity_realizability import realizable_choice  # noqa: E402

OUT = ROOT / "data" / "w33_pass11088_baryonic_hidden_condensates.json"
ALL_KINDS = ("meson", "baryon", "antibaryon", "su2pair", "six2")
MODELS = ("z3z3|Z3Z3_0001_c3__SM_20260934_739", "z3z3|Z3Z3_0001_c4__SM_20260935_3729",
          "z3z3|Z3Z3_0001_c1__SM_20260932_2822")


def _ray_job_dp(args):
    """exact realizability of one FI ray by the subset-sum programme (validated against the depth-first search:
    60/60 on the c4 quadratic census, 39/39 on c3)."""
    choices, odd_vecs, even_vecs, nq, orders = args
    k = realizable_choice(nq, orders, odd_vecs, even_vecs, choices)
    return None if k is None else [choices[u][ku] for u, ku in enumerate(k)]


def census_all_kinds(name, model, disc):
    """Pass 11024's condensate census with every single-factor invariant in scope."""
    orig, orig_job = P.augment.__defaults__, P._ray_job
    P.augment.__defaults__ = (ALL_KINDS,)
    P.condensate_census.__globals__["_ray_job"] = _ray_job_dp
    try:
        return P.condensate_census(name, model, disc)
    finally:
        P.augment.__defaults__ = orig
        P.condensate_census.__globals__["_ray_job"] = orig_job


def candidate_full(blob, name, support):
    Z = P.model_Z(blob, name)
    T = P.P79.Tools(Z)
    dims = {n: Z.raw[n]["dim"].split(",") for n in Z.raw}
    _, qd, c3, w2 = P.slots(blob["ledger"][name])
    hid = [i for i in range(len(qd)) if i not in (c3, w2)]
    present = defaultdict(set)
    for d in dims.values():
        for i, x in enumerate(d):
            present[i].add(x)
    Nof = {}
    for i in hid:
        cp = sorted(P.dimof(v) for v in present[i] if not v.startswith("-") and P.dimof(v) > 1 and "-" + v in present[i])
        Nof[i] = cp[0] if cp else (2 if "2" in present[i] else None)
    comps = [f for f in support if "of" in f]
    singles = [f["name"] for f in support if "of" not in f]
    slot_of = lambda n: [i for i in hid if P.dimof(dims[n][i]) != 1][0]
    kinds_by_slot = defaultdict(set)
    for f in comps:
        kinds_by_slot[slot_of(f["of"][0])].add(f["kind"])
    full = {s for s, ks in kinds_by_slot.items() if ks - {"meson"} or Nof[s] == 2}
    aligned = {s for s in kinds_by_slot if s not in full}
    realw = lambda ev, od: P.M.realizable([list(t) for t in od] + [list(t) for t in ev] + [list(t) for t in Z.virtual]
                                          + [Z.wvec], list(range(len(od))))
    base = [Z.vec[s] for s in singles] + [[sum(Z.vec[c][k] for c in f["of"]) for k in range(len(Z.wvec))] for f in comps]
    parity = P.P8.match_search(realw, base, T.species, T.leptons)[0]
    S = singles + [x for f in comps for x in f["of"]]
    cen = [i for i in hid if Nof[i]]
    vx = lambda n: list(Z.vec[n]) + [(P._centre(dims[n][i], Nof[i]) or Fraction(0)) for i in cen]
    wv = list(Z.wvec) + [Fraction(0)] * len(cen)
    E = P.EXO.ExactOrders([vx(s) for s in S], Z.nq,
                          [x.denominator for v in list(Z.vec.values()) + [Z.wvec] + Z.virtual for x in v] + [2, 3, 4])
    lab = [f for f in Z.left if f["base"] in P.P.SMY]
    y = P.P.solve_affine([f["q"] for f in lab], [P.P.SMY[f["base"]] for f in lab])[0]
    Yv = {f["name"]: P.P.dot(y, f["q"]) for f in Z.left}

    def components(tok, i):
        if i in full:
            return ["1"] * P.dimof(tok)
        return P._components(tok, Nof[i], i in aligned)
    cls = defaultdict(list)
    for f in Z.left:
        n, d = f["name"], dims[f["name"]]
        if P.dimof(d[c3]) == 1 and P.dimof(d[w2]) == 1 and Yv[n] == 0:
            continue
        parts = [[]]
        for i in hid:
            parts = [p + [c] for p in parts for c in components(d[i], i)]
        for p in parts:
            cls[(d[c3], d[w2], Yv[n], tuple(p))].append(n)

    def conj(k):
        hp = []
        for i, x in zip(hid, k[3]):
            if i in full or P.dimof(x) == 1:
                hp.append(x)
            elif i in aligned:
                hp.append(x if Nof[i] == 2 else (x[1:] if x.startswith("-") else "-" + x))
            else:
                hp.append(P.conj_token(x, present[i], Nof[i] == 2))
        return (P.conj_token(k[0], present[c3], False), k[1], -k[2], tuple(hp))
    cq = qd[c3]
    light = light_sym = 0
    seen = set()
    for K in sorted(cls, key=str):
        Kb = conj(K)
        if K in seen or Kb in seen:
            continue
        seen.update([K, Kb])
        wk = P.dimof(K[1])
        charges = [K[2] + Fraction(t_, 2) for t_ in range(-(wk - 1), wk, 2)]
        col = P.dimof(K[0])
        frac = (col == 1 and any(c.denominator != 1 for c in charges)) or \
            (col == 3 and any((c - (Fraction(2, 3) if K[0] == cq else Fraction(1, 3))).denominator != 1 for c in charges))
        if not frac:
            continue
        A_, B_ = cls[K], cls.get(Kb, [])
        rk = rl = 0
        if B_:
            rk = P.assignment_rank(np.array([[int(E.coupling([vx(a), vx(b)], wv) is not None) for b in B_] for a in A_]))
            rl = P.assignment_rank(np.array([[int(E.coupling_lattice([vx(a), vx(b)], wv)) for b in B_] for a in A_]))
        light += len(A_) + len(B_) - 2 * rk
        light_sym += len(A_) + len(B_) - 2 * rl
    return dict(parity_composites_even=parity, full_breaking_slots=sorted(full), aligned_slots=sorted(aligned),
                kinds=sorted({f["kind"] for f in comps}), light_fractional=light,
                light_fractional_protected_by_symmetry=light_sym, W_on_vacuum=E.coupling([], wv))


def main():
    from multiprocessing import Pool
    blob = P.P79.load()
    res = {}
    with Pool(8) as pool:
        P._POOL = pool
        for name in MODELS:
            b = name.split("|")[1]
            r = census_all_kinds(name, blob["ledger"][name], blob["disc"][b])
            r.pop("_vecs", None)
            # MSSM parity + analysis for every realizable support that contains a baryonic (or any new) condensate
            sups = r.get("realizable_supports", [])
            new = [s for s in sups if any(f.get("kind") in ("baryon", "antibaryon", "six2") for f in s["fields"])]
            cands = []
            for s in new:
                c = candidate_full(blob, name, s["fields"])
                c["fields"] = s["fields"]
                cands.append(c)
            summ = Counter((c["parity_composites_even"], c["light_fractional"] > 0) for c in cands)
            res[name] = dict(census=r["census"], curing_slots=r.get("curing_slots"), composites=r["composites"],
                             kinds=r["kinds"], realizable_supports=len(sups), with_baryonic=len(new),
                             candidates=cands, summary={str(k): v for k, v in summ.items()},
                             min_light_fractional_viable=min((c["light_fractional"] for c in cands
                                                              if c["parity_composites_even"] == "viable"), default=None))
            print(name.split("|")[1], r["composites"], r["census"], "baryonic supports", len(new), dict(summ),
                  "min light (viable)", res[name]["min_light_fractional_viable"], flush=True)
    OUT.write_text(json.dumps(dict(pass_id=11088, models=res), indent=1, sort_keys=True, default=str))


if __name__ == "__main__":
    main()
