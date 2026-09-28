#!/usr/bin/env python3
"""Pass 11090: blowing up the fixed points cannot free the massless fractional charges.

A vacuum in which twisted singlets take values is the orbifold limit of a resolved (blown-up) compactification: the
singlet VEVs are the blow-up moduli of the fixed points they live at.  At a generic point of the resolution the
fixed-point labels lose their meaning, so the space-group selection rules (and with them the point-group labels and
the discrete R-rotations) need not survive; the gauge charges -- the U(1)s and the hidden centres -- always do.
Dropping selection rules can only ADD couplings, so the number of massless states computed with a weaker rule set is
a lower bound for every resolution that preserves the rules kept.

For every Z3xZ3 parity vacuum (the 55 hidden-condensate vacua and the 28 singlet vacua of Passes 11024/10979) the exact
holomorphic mass ranks of every SM-charged vector-like class are recomputed with four rule sets:
    all            U1 + PG + SG + R + centres
    no SG          drop the space-group translations (fixed-point labels)
    no SG, PG      drop the point-group sector labels as well
    gauge only     U1 + centres
and the fractionally charged states left massless are counted.

Why no resolution preserving the Standard Model can reach them: on every fixed point (sector + translation labels) of
all 13 Z3xZ3 Standard Models the charge class c = 3Q + t (mod 3), t the colour triality, is CONSTANT, and it is a
linear Z3 character of the space group (unique per model).  Every state at a fixed point with c != 0 is fractionally
charged, so those fixed points host no neutral field: they can be blown up only by VEVs that break colour or
electromagnetism.
"""
from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11024_z3xz3_susy_vacua_fractional_charge_index as P  # noqa: E402
import w33_pass11087_unbroken_three_group_protection as Q  # noqa: E402

OUT = ROOT / "data" / "w33_pass11090_resolution_cannot_free_fractional_charges.json"
RULES = {"all": ("U1", "PG", "SG", "R", "centre"), "no SG": ("U1", "PG", "R", "centre"),
         "no SG, PG": ("U1", "R", "centre"), "gauge only": ("U1", "centre")}


def analyse(job):
    name, S, tag, broken_fields = job
    blob = Q._blob()
    Z, dims, vx, wv, blocks, per = Q.setup(name, S)
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
    slot_of = lambda n: [i for i in hid if P.dimof(dims[n][i]) != 1][0]
    broken = sorted({slot_of(c) for c in broken_fields})
    lab = [f for f in Z.left if f["base"] in P.P.SMY]
    y = P.P.solve_affine([f["q"] for f in lab], [P.P.SMY[f["base"]] for f in lab])[0]
    Yv = {f["name"]: P.P.dot(y, f["q"]) for f in Z.left}
    cls = defaultdict(list)
    for f in Z.left:
        n, d = f["name"], dims[f["name"]]
        if P.dimof(d[c3]) == 1 and P.dimof(d[w2]) == 1 and Yv[n] == 0:
            continue
        parts = [[]]
        for i in hid:
            cs = P._components(d[i], Nof[i], i in broken)
            parts = [p + [c] for p in parts for c in cs]
        for p in parts:
            cls[(d[c3], d[w2], Yv[n], tuple(p))].append(n)

    def conj(k):
        hp = []
        for i, x in zip(hid, k[3]):
            if i in broken:
                hp.append(x if (Nof[i] == 2 or P.dimof(x) == 1) else (x[1:] if x.startswith("-") else "-" + x))
            else:
                hp.append(P.conj_token(x, present[i], Nof[i] == 2))
        return (P.conj_token(k[0], present[c3], False), k[1], -k[2], tuple(hp))
    cq = qd[c3]
    frac_of = {}
    for K in cls:
        wk = P.dimof(K[1])
        charges = [K[2] + Fraction(t_, 2) for t_ in range(-(wk - 1), wk, 2)]
        col = P.dimof(K[0])
        frac_of[K] = (col == 1 and any(c.denominator != 1 for c in charges)) or \
            (col == 3 and any((c - (Fraction(2, 3) if K[0] == cq else Fraction(1, 3))).denominator != 1 for c in charges))
    res = {}
    isolated_loc = Counter()
    vac_loc = Counter((tuple(vx(s)[Z.nq:Z.nq + 2]), tuple(vx(s)[Z.nq + 2:Z.nq + 5])) for s in S)
    for rname, keep_blocks in RULES.items():
        keep = [i for b in keep_blocks for i in blocks[b]]
        nqk = Z.nq if "U1" in keep_blocks else 0
        E = P.EXO.ExactOrders([[vx(s)[i] for i in keep] for s in S], nqk, [2, 3, 4, 9])
        wk_ = [wv[i] for i in keep]
        light, seen = 0, set()
        for K in sorted(cls, key=str):
            Kb = conj(K)
            if K in seen or Kb in seen:
                continue
            seen.update([K, Kb])
            if not frac_of[K]:
                continue
            A_, B_ = cls[K], cls.get(Kb, [])
            M = np.array([[int(E.coupling([[vx(a)[i] for i in keep], [vx(b)[i] for i in keep]], wk_) is not None)
                           for b in B_] for a in A_]) if B_ else np.zeros((len(A_), 0), dtype=int)
            rk = P.assignment_rank(M) if B_ else 0
            light += len(A_) + len(B_) - 2 * rk
            if rname == "all" and B_:
                for i_, a in enumerate(A_):
                    if not M[i_].any():
                        loc = (tuple(vx(a)[Z.nq:Z.nq + 2]), tuple(vx(a)[Z.nq + 2:Z.nq + 5]))
                        isolated_loc["at a fixed point the vacuum resolves" if loc in vac_loc else "elsewhere"] += 1
        res[rname] = light
    return dict(model=name, tag=tag, light_fractional=res, isolated_states=dict(isolated_loc))


def charge_character(name):
    """the electric-charge class c = 3Q + t (mod 3), t = colour triality, on every fixed point (sector + translation
    labels): constant per fixed point?  a linear Z3 character of the space group?  which fixed points are fractional,
    and do any SM-singlet fields (of any hidden representation) live there?"""
    import itertools
    blob = Q._blob()
    Z = P.model_Z(blob, name)
    lab = [f for f in Z.left if f["base"] in P.P.SMY]
    y = P.P.solve_affine([f["q"] for f in lab], [P.P.SMY[f["base"]] for f in lab])[0]
    dims = {n: Z.raw[n]["dim"].split(",") for n in Z.raw}
    _, qd, c3, w2 = P.slots(blob["ledger"][name])
    cq = qd[c3]
    byfp = defaultdict(set)
    singlet_at = defaultdict(int)
    for f in Z.left:
        n, d = f["name"], dims[f["name"]]
        wk = P.dimof(d[w2])
        Qn = P.P.dot(y, Z.vec[n][:Z.nq]) + Fraction(wk - 1, 2)
        tri = 0 if P.dimof(d[c3]) == 1 else (1 if d[c3] == cq else -1)
        l5 = tuple(int(x * 3) % 3 for x in Z.vec[n][Z.nq:Z.nq + 5])
        byfp[l5].add(int((3 * Qn + tri) % 3))
        if P.dimof(d[c3]) == 1 and wk == 1 and P.P.dot(y, Z.vec[n][:Z.nq]) == 0:
            singlet_at[l5] += 1
    constant = all(len(v) == 1 for v in byfp.values())
    cval = {l: next(iter(v)) for l, v in byfp.items() if len(v) == 1}
    chars = [a for a in itertools.product(range(3), repeat=5)
             if all(sum(ai * li for ai, li in zip(a, l)) % 3 == c for l, c in cval.items())]
    frac_fp = [l for l, c in cval.items() if c != 0]
    return dict(fixed_points=len(byfp), constant_on_every_fixed_point=constant, fractional_fixed_points=len(frac_fp),
                space_group_characters=[list(a) for a in chars],
                sm_singlets_at_fractional_fixed_points=sum(singlet_at[l] for l in frac_fp),
                untwisted_class=cval.get((0, 0, 0, 0, 0)))


def jobs():
    cert = json.loads(Q.CERT.read_text())
    cert79 = json.loads(Q.CERT79.read_text())
    out = []
    for name, D in cert["D"]["census"].items():
        if "composite_vacua" not in D:
            continue
        for s in D["composite_vacua"]["supports"]:
            if s["mssm_parity"] != "viable":
                continue
            f = s["fields"]
            S = [x["name"] for x in f if "of" not in x] + [c for x in f if "of" in x for c in x["of"]]
            out.append((name, S, "condensate", [x["of"][0] for x in f if "of" in x]))
    for name, vs in cert79["E"]["vacua"].items():
        for v in vs:
            if "dflat_full" in v:
                out.append((name, v["support"], "singlet", []))
    return out


def main():
    from multiprocessing import Pool
    chars = {name: charge_character(name) for name in sorted(Q._blob()["ledger"])}
    for k, v in chars.items():
        print(k.split("|")[1], v, flush=True)
    if "--characters-only" in sys.argv and OUT.exists():
        d = json.loads(OUT.read_text())
        d["charge_character"] = chars
        OUT.write_text(json.dumps(d, indent=1, sort_keys=True, default=str))
        return
    with Pool(8) as pool:
        res = pool.map(analyse, jobs(), chunksize=1)
    summ = {}
    for tag in ("condensate", "singlet"):
        rs = [r for r in res if r["tag"] == tag]
        summ[tag] = dict(vacua=len(rs),
                         min_light_fractional={k: min(r["light_fractional"][k] for r in rs) for k in RULES},
                         max_light_fractional={k: max(r["light_fractional"][k] for r in rs) for k in RULES},
                         isolated=dict(sum((Counter(r["isolated_states"]) for r in rs), Counter())))
        print(tag, json.dumps(summ[tag]), flush=True)
    OUT.write_text(json.dumps(dict(pass_id=11090, rules=RULES, summary=summ, vacua=res, charge_character=chars),
                              indent=1, sort_keys=True, default=str))


if __name__ == "__main__":
    main()
