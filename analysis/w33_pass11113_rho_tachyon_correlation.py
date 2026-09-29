#!/usr/bin/env python3
"""Pass 11113: the 'tachyon-free at rho => untwisted quarks' correlation at scale -- suggestive, not significant.

Pass 11108 found every model that stays tachyon-free with its Wilson-line tori at the SU(3) point rho to have untwisted
quarks.  Here the quark sector of ALL 387 new models of the Pass 11108 rescan is classified (tree-level up-type Yukawa
triples, Pass 11098 engine, on full spectra from both engines), and combined with the 104 models of Pass 11095:

                       twisted up-quarks   untwisted up-quarks
    rho-safe                 0                    5
    tachyonic at rho       152                  334

Under independence the expected number of rho-safe twisted models is 152 x 5/491 = 1.55, so P(0) = exp(-1.55) = 0.21
(Poisson; Fisher exact one-sided p = 0.16).  The correlation is SUGGESTIVE BUT NOT SIGNIFICANT: the twisted-quark sector
is not proven to be tachyonic at rho; the rho-safe rate is simply low (1.0%) everywhere.  (Pass 11108's '0/31' was a
convenience sample.)  The rho-level spectrum is universal: Delta in {-1/18, -1/6, none} for all 491 models.
"""
from __future__ import annotations

import json
from math import comb, exp
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SECT = ROOT / "data" / "w33_pass11113_quark_sectors_387.json"
LEV104 = ROOT / "data" / "w33_pass11108_tachyon_levels_104.json"
YUK104 = ROOT / "data" / "w33_pass11098_yukawa_104.json"
OUT = ROOT / "data" / "w33_pass11113_rho_tachyon_correlation.json"


def table():
    t = {('twisted', True): 0, ('twisted', False): 0, ('untwisted', True): 0, ('untwisted', False): 0}
    s = json.loads(SECT.read_text())
    for v in s.values():
        cls = 'twisted' if v['up_triples'] and v['up_untwisted'] == 0 else 'untwisted'
        t[(cls, v['rho'] is None)] += 1
    lev, yk = json.loads(LEV104.read_text()), json.loads(YUK104.read_text())
    for m, v in yk.items():
        cls = 'twisted' if v['up']['untwisted'] == 0 else 'untwisted'
        t[(cls, lev[m]['rho'] is None)] += 1
    return t


def fisher_one_sided(a, b, c, d):
    """P(X <= a) for the (twisted, safe) cell, margins fixed"""
    n1, k, n = a + b, a + c, a + b + c + d
    return sum(comb(k, x) * comb(n - k, n1 - x) for x in range(0, a + 1)) / comb(n, n1)


def summarize():
    t = table()
    a, b = t[('twisted', True)], t[('twisted', False)]
    c, d = t[('untwisted', True)], t[('untwisted', False)]
    n = a + b + c + d
    expct = (a + b) * (a + c) / n
    res = dict(pass_id=11113, table={f"{k[0]}|rho_safe={k[1]}": v for k, v in t.items()}, models=n,
               expected_twisted_safe=expct, poisson_p0=exp(-expct), fisher_one_sided=fisher_one_sided(a, b, c, d),
               significant_at_5pct=fisher_one_sided(a, b, c, d) < 0.05)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    summarize()
