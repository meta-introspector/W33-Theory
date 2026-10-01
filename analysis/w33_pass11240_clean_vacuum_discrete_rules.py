#!/usr/bin/env python3
"""Pass 11240: does Pass 11232's clean U(1)' mu-vacuum (model Z6II_06_1204) survive the discrete selection rules?

Pass 11232 found one vacuum, support {n_2, n_13, n_41, n_55, n_64} of model Z6-II|Z6II_06__SM_20260917_1204, in which
an unbroken U(1)' forbids every mu = bl_1 l_j to all orders while the top Yukawa q_1 bl_1 bu is allowed and every other
Higgs pair and vector-like exotic can be massive -- at the level of the gauge U(1)s.  The string has more rules: the
space-group charges (Z6 x Z3 x Z2 x Z2 for this Z6-II lattice) and the geometric R-symmetries Z6^R x Z3^R x Z2^R with
superpotential charge -1 each (frozen in Pass 10967; used in Pass 10968).  The singlets named in the support stand for
U(1) TYPES; several fields share each type and differ in discrete charges, so every choice of condensing field is
enumerated exhaustively.

For each choice the full vacuum lattice of Pass 11241 decides every coupling exactly (forbidden at all orders) or as a
necessary condition (allowed).  Three rule sets are compared: gauge U(1)s with their discrete remnants only; + space
group; + space group + R.  Also per choice: lattice F-flatness (sufficient).  Matter parity is not re-tested here
(Pass 10967/10968: no surviving parity in this model).
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass10960_matter_even_dflat_closure as P  # noqa: E402
import w33_pass11241_discrete_mu_symmetry as D  # noqa: E402

OUT = ROOT / "data" / "w33_pass11240_clean_vacuum_discrete_rules.json"
NAME = "Z6-II|Z6II_06__SM_20260917_1204"
SUPPORT_TYPES = ["n_13", "n_2", "n_41", "n_55", "n_64"]


class Restricted(D.Model):
    """the same model with some coordinates (space group and/or R) switched off"""

    def __init__(self, name, model, disc, keep_space_group=True, keep_R=True):
        super().__init__(name, model, disc)
        nq = self.nq
        nsg = len(self.orders) - len(self.R_orders)
        keep = list(range(nq)) + [nq + i for i in range(nsg) if keep_space_group] + \
            [nq + nsg + i for i in range(len(self.R_orders)) if keep_R]
        self.vec = {k: tuple(v[i] for i in keep) for k, v in self.vec.items()}
        self.w = tuple(self.w[i] for i in keep)
        self.virtual = [tuple(v[i] for i in keep) for v in self.virtual if any(v[i] for i in keep)]


def run():
    ledger, _ = P.load_ledger()
    disc = D.M67.load(D.DISC)
    model = ledger[NAME]
    full = D.Model(NAME, model, disc)
    types = full.u1_types()
    byname = {f["name"]: tuple(f["q"]) for f in full.left}
    choices = [types[byname[n]] for n in SUPPORT_TYPES]
    rulesets = {"U(1) lattice only": dict(keep_space_group=False, keep_R=False),
                "+ space group": dict(keep_space_group=True, keep_R=False),
                "+ space group + R": dict(keep_space_group=True, keep_R=True)}
    models = {k: Restricted(NAME, model, disc, **v) for k, v in rulesets.items()}
    res = dict(pass_id=11240, model=NAME, support_types=SUPPORT_TYPES,
               fields_per_type=[len(c) for c in choices], choices=1)
    for c in choices:
        res["choices"] *= len(c)
    table = {k: dict(protected=0, clean=0, clean_F_flat=0, light=dict()) for k in rulesets}
    examples = {}
    for S in itertools.product(*choices):
        for k, Mo in models.items():
            r = Mo.analyse_vacuum(S)
            if r is None:
                continue
            t = table[k]
            t["protected"] += 1
            t["clean"] += r["clean"]
            t["clean_F_flat"] += r["clean"] and r["F_flat_lattice"]
            for dd in r["detail"]:
                if dd["Hu"] == "bl_1":
                    for b, n in dd["light_vectorlike"].items():
                        if n:
                            t["light"][b] = t["light"].get(b, 0) + 1
                    if not dd["other_Hu_all_massive"]:
                        t["light"]["extra Higgs"] = t["light"].get("extra Higgs", 0) + 1
            if r["clean"] and k not in examples:
                examples[k] = r
    res["by_ruleset"] = table
    res["clean_examples"] = examples
    return res


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({k: v for k, v in res.items() if k != "clean_examples"}, indent=1))
    for k, v in res["clean_examples"].items():
        print(k, json.dumps(v)[:800])


if __name__ == "__main__":
    main()
