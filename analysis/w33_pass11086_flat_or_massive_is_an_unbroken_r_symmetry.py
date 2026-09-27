#!/usr/bin/env python3
"""Pass 11086: the 'flat-or-massive' pattern is an unbroken R-symmetry, not a law.

Since Pass 10968 the programme has reported that parity-preserving vacua with W|_S = 0 keep an exotic vector-like
colour triplet massless, and the manuscript says the triplet is massless *because* the superpotential vanishes.
This pass tests the mechanism exactly, with the lattice solver of Pass 11024, on every parity vacuum with a
decided exotic spectrum:
  * Z12-I, 14 vacua (Pass 10978), under the weakest rules and under the admissible Z2^R;
  * Z3xZ3, the 36 MSSM-viable minimal supports and the 28 extended vacua (Pass 10979).
For each vacuum:
  W status      -- W|_S has a monomial (nonzero), or is forbidden even with integer exponents of any sign
                   ('lattice-forbidden': an exact unbroken symmetry of the vacuum under which W is charged, i.e.
                   an unbroken discrete R-symmetry), or only by non-negativity ('holomorphy');
  exotic d      -- the holomorphic mass rank of the d.bd matrix and its rank ignoring holomorphy;
  R-neutral     -- d.bd pairs whose total charge lies in the vacuum lattice: their mass term carries exactly W's
                   charge, so it is forbidden iff W is.
Result: where both W and the exotic masses vanish, both are (almost always) forbidden by the same unbroken
symmetry; but neither implies the other -- exotics are massless with W != 0 and massive with W lattice-forbidden.
"""
from __future__ import annotations

import gzip
import importlib.util
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / "analysis" / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


EXO = _load("exact_orders", "w33_exact_monomial_orders.py")
P78 = _load("p10978", "w33_pass10978_z12_e6_admissible_r_rule.py")
P79 = _load("p10979", "w33_pass10979_z3xz3_parity_vacua_f_obstruction.py")
OUT = ROOT / "data" / "w33_pass11086_flat_or_massive_is_an_unbroken_r_symmetry.json"


def rank(B):
    if B.size == 0:
        return 0
    r, c = linear_sum_assignment(-B)
    return int(B[r, c].sum())


def classify(vec, W, nq, dens, S, d_list, bd_list):
    E = EXO.ExactOrders([vec[s] for s in S], nq, dens)
    if E.coupling([], W) is not None:
        w = "nonzero"
    elif E.coupling_lattice([], W):
        w = "zero (holomorphy)"
    else:
        w = "zero (unbroken symmetry)"
    hol = rank(np.array([[int(E.coupling([vec[a], vec[b]], W) is not None) for b in bd_list] for a in d_list]))
    lat = rank(np.array([[int(E.coupling_lattice([vec[a], vec[b]], W)) for b in bd_list] for a in d_list]))
    neutral = sum(E.lattice_solvable([-(x + y) for x, y in zip(vec[a], vec[b])]) for a in d_list for b in bd_list)
    n = len(d_list)
    if hol == n:
        ex = "massive"
    elif lat < n:
        ex = "massless (unbroken symmetry)"
    else:
        ex = "massless (holomorphy)"
    return dict(W=w, exotic=ex, rank_holomorphic=hol, rank_lattice=lat, n_exotic=n, r_neutral_pairs=int(neutral),
                pairs=n * len(bd_list))


def main():
    out, table = {}, {}
    # Z12-I
    z12 = json.load(gzip.open(P78.Z12, "rt"))
    prev = json.loads(P78.P10968.read_text())["z12_models"]
    for tag, flag in (("Z12-I weakest", False), ("Z12-I Z2^R", True)):
        cnt = Counter()
        for name, r in sorted(prev.items()):
            if "support" not in r:
                continue
            m = P78.Z12Model(z12[name], flag)
            c = classify(m.V, m.W, m.nq, [m.den], r["support"], m.lab("d"), m.lab("bd"))
            out.setdefault(tag, {})[name] = c
            cnt[(c["W"], c["exotic"])] += 1
        table[tag] = {" | ".join(k): v for k, v in sorted(cnt.items())}
        print(tag, table[tag], flush=True)
    # Z3xZ3
    blob = P79.load()
    cert = json.loads((ROOT / "data" / "w33_pass10979_z3xz3_parity_vacua_f_obstruction.json").read_text())
    for tag, source in (("Z3xZ3 minimal", "C"), ("Z3xZ3 extended", "E")):
        cnt = Counter()
        items = cert["C"].items() if source == "C" else cert["E"]["vacua"].items()
        for name, r in items:
            b = name.split("|")[1]
            Z = P79.P8.Z6IIR(name, {name: blob["ledger"][name]}, {b: blob["disc"][b]}, {})
            dens = [x.denominator for v in list(Z.vec.values()) + [Z.wvec] + Z.virtual for x in v]
            rows = r["supports"] if source == "C" else r
            for s in rows:
                if source == "C" and s["mssm_parity"] != "viable":
                    continue
                if source == "E" and "dflat_full" not in s:
                    continue
                c = classify(Z.vec, Z.wvec, Z.nq, dens, s["support"], Z.lab("d"), Z.lab("bd"))
                out.setdefault(tag, {}).setdefault(name, []).append(c)
                cnt[(c["W"], c["exotic"])] += 1
        table[tag] = {" | ".join(k): v for k, v in sorted(cnt.items())}
        print(tag, table[tag], flush=True)
    OUT.write_text(json.dumps(dict(pass_id=11086, table=table, vacua=out), indent=1, sort_keys=True, default=str))


if __name__ == "__main__":
    main()
