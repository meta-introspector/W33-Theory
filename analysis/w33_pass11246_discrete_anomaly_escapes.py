#!/usr/bin/env python3
"""Pass 11246: can a DISCRETE vacuum symmetry forbid mu while every colored state is massive?

Pass 11242 proved that an anomaly-free continuous U(1)' forbidding mu leaves a colored state massless (812 cases, 0
violations).  For a discrete symmetry the same pair-counting gives only a congruence -- the anomaly coefficient vanishes
mod N, or equals a universal Green-Schwarz value -- so the argument can fail.  Here the full vacuum symmetry (Pass
11241's exact integer lattice: U(1)s with their discrete remnants, space group, and for Z6-II the geometric R-symmetries)
is used on sampled FI-cancelling vacua, and for every mu-protected H_u the colored perfect matching of Pass 11242 is
tested with lattice-allowed mass terms (vector-like, and Yukawas with H_u or some H_d).  An "escape" is a vacuum where
mu is forbidden to all orders and yet every colored component can be massive.
"""
from __future__ import annotations

import json
import random
import re
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass10960_matter_even_dflat_closure as P  # noqa: E402
import w33_pass11232_mu_vacuum_symmetry as MU  # noqa: E402
import w33_pass11241_discrete_mu_symmetry as D  # noqa: E402

OUT = ROOT / "data" / "w33_pass11246_discrete_anomaly_escapes.json"
ab = lambda s: abs(int(re.match(r"(-?\d+)", s).group(1)))


def colored_components(model, Mo):
    raw = {f["name"]: f["dim"].split(",") for f in model["left"]}
    q0, l0 = raw[Mo.byb["q"][0]], raw[Mo.byb["l"][0]]
    c3 = [i for i, x in enumerate(q0) if ab(x) == 3][0]
    c2 = [i for i, x in enumerate(l0) if ab(x) == 2 and ab(q0[i]) == 2][0]
    hkey = lambda n: tuple(ab(x) for i, x in enumerate(raw[n]) if i not in (c3, c2))

    def comps(b):
        out = []
        for n in Mo.byb.get(b, []):
            m = 1
            for x in hkey(n):
                m *= x
            out += [(b, n)] * m
        return out
    return comps, hkey


def escapes_in_vacuum(Mo, comps, hkey, S):
    L = D.Lattice([Mo.vec[s] for s in S] + Mo.virtual)
    ok = lambda *names: Mo.target(list(names)) in L
    hu, hd = Mo.byb.get("bl", []), Mo.byb.get("l", [])
    up_r, up_c = comps("q") + comps("u"), comps("bq") + comps("bu")
    dn_r, dn_c = comps("q") + comps("d"), comps("bq") + comps("bd")
    out = []
    for H in hu:
        if any(ok(H, l) for l in hd):
            continue
        if not any(ok(q, H, u) for q in Mo.byb.get("q", []) for u in Mo.byb.get("bu", [])):
            continue
        esc = False
        for Hd in hd:
            def up(i, j):
                (ta, a), (tb, b) = up_r[i], up_c[j]
                if hkey(a) != hkey(b):
                    return False
                return ok(a, b, H) if (ta, tb) == ("q", "bu") else (ok(a, b, Hd) if (ta, tb) == ("u", "bq") else ok(a, b))

            def dn(i, j):
                (ta, a), (tb, b) = dn_r[i], dn_c[j]
                if hkey(a) != hkey(b):
                    return False
                return ok(a, b, Hd) if (ta, tb) == ("q", "bd") else (ok(a, b, H) if (ta, tb) == ("d", "bq") else ok(a, b))
            if (len(up_r) == len(up_c) and MU.structural_rank(range(len(up_r)), range(len(up_c)), up) == len(up_r)
                    and len(dn_r) == len(dn_c)
                    and MU.structural_rank(range(len(dn_r)), range(len(dn_c)), dn) == len(dn_r)):
                esc = True
                break
        out.append((H, esc))
    return out


def run(ray_cap=60, choice_cap=2, seed=11246):
    rng = random.Random(seed)
    ledger, sha = P.load_ledger()
    disc = D.M67.load(D.DISC)
    res, t0 = {}, time.time()
    for name in sorted(ledger):
        model = ledger[name]
        Mo = D.Model(name, model, disc)
        if Mo.tr0 == 0:
            continue
        T, types, rays = D.model_rays(Mo)
        if rays is None:
            continue
        comps, hkey = colored_components(model, Mo)
        sups = sorted({tuple(k for k, v in enumerate(r) if v != 0) for r in rays})
        if len(sups) > ray_cap:
            sups = rng.sample(sups, ray_cap)
        prot = esc = vac = 0
        examples = []
        for sup in sups:
            for _ in range(choice_cap):
                S = tuple(rng.choice(types[T[k]]) for k in sup)
                vac += 1
                for H, e in escapes_in_vacuum(Mo, comps, hkey, S):
                    prot += 1
                    esc += e
                    if e and len(examples) < 3:
                        examples.append(dict(support=list(S), Hu=H))
        res[name] = dict(vacua=vac, protected=prot, escapes=esc, examples=examples,
                         rules="U(1)+remnants+space group" + ("+R" if Mo.R_orders else ""))
        print(name, res[name]["vacua"], res[name]["protected"], res[name]["escapes"], flush=True)
    summary = dict(models=len(res), protected=sum(v["protected"] for v in res.values()),
                   escapes=sum(v["escapes"] for v in res.values()),
                   models_with_escape=sorted(k for k, v in res.items() if v["escapes"]),
                   z6i_protected=sum(v["protected"] for k, v in res.items() if k.startswith("Z6-I|")),
                   sampling=dict(ray_cap=ray_cap, choice_cap=choice_cap, seed=seed), ledger_sha256=sha,
                   seconds=round(time.time() - t0, 1))
    return dict(pass_id=11246, summary=summary, models=res)


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res["summary"], indent=1))


if __name__ == "__main__":
    main()
