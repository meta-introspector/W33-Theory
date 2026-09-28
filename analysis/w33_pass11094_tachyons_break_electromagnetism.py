#!/usr/bin/env python3
"""Pass 11094: what the tachyons of the 87 non-SUSY Z6-I twins carry -- every one is fractionally charged or coloured.

Two independent computations of the theta-sector tachyons of Pass 11092:
  A. the patched orbifolder at mass level mu = -1/6 (PATCH 7: every non-identity sector solved at M^2/8 = mu, through the
     full projection).  Validation: its theta and theta^5 counts equal the pure-Python lemma counts of Pass 11092 in
     87/87 twins; the SUSY parents give no twisted state at that level; and on a tachyonic control (SO(16)xE8 on
     T6/Z3, standard embedding) it finds the same single (10,1,1) tachyon as the independent non-SUSY orbifolder.
     Its tachyon fields, with the SM hypercharge carried over from each twin's SM configuration (exact linear algebra
     between the two U(1) bases), are frozen in data/w33_pass11094_tachyon_charges_orbifolder.json;
  B. pure Python, weight by weight: every tachyonic momentum p_sh (1/2 p_sh^2 = 7/12 with N_L = 0, or 5/12 with one
     alpha_{-1/6}) at every theta fixed point, evaluated against the exact 16D hypercharge vector t_Y and the colour and
     SU(2)_L simple roots of the same model: Y = t_Y.p_sh, colour Dynkin labels, T3 = (alpha_2.p_sh)/2.
Result (both): every tachyonic state is colour-charged or has fractional electric charge (+-1/2); none is colourless
and electrically neutral.  So any tachyon condensate breaks U(1)_EM (and colour in 27 twins).
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11092_nonsusy_z6i_twins_are_tachyonic as P  # noqa: E402

ORB = ROOT / "data" / "w33_pass11094_tachyon_charges_orbifolder.json"
OUT = ROOT / "data" / "w33_pass11094_tachyons_break_electromagnetism.json"


def tachyon_momenta(V, W5):
    """all tachyonic theta-sector momenta p_sh with the half-norm and oscillator number (Pass 11092 lemma levels)"""
    osc = P.osc_counts(P.OSC_THETA, F(1, 2))
    out = []
    for m in range(3):
        Vg = [v + m * w for v, w in zip(V, W5)]
        for h in (F(7, 12), F(5, 12), F(1, 4), F(1, 12)):
            if F(7, 12) - h not in osc:
                continue
            A = list(P.e8_near(Vg[:8], 2 * h))
            B = {}
            for p, s in P.e8_near(Vg[8:], 2 * h):
                B.setdefault(s, []).append(p)
            for p, s in A:
                for q in B.get(2 * h - s, []):
                    out.append((m, h, [a + b for a, b in zip(p + q, Vg)]))
    return out


def classify(p_sh, tY, col, su2):
    dot = lambda u, v: sum(a * b for a, b in zip(u, v))
    Y = dot(tY, p_sh)
    c = tuple(dot(r, p_sh) for r in col)
    t3 = dot(su2[0], p_sh) / 2
    Q = Y + t3
    coloured = any(x != 0 for x in c)
    return dict(Y=Y, colour_labels=c, T3=t3, Q=Q, coloured=coloured, neutral_colourless=(not coloured and Q == 0),
                fractional=(not coloured and Q.denominator != 1))


def main():
    orb = json.loads(ORB.read_text())
    models = {f"{m['index']:02d}": m for m in json.loads(P.MODELS.read_text())["models"]}
    res = {}
    for key, m in models.items():
        o = orb[str(int(key))]
        tY = [F(x) for x in o["tY"]]
        col = [[F(x) for x in r] for r in o["colour_roots"]]
        su2 = [[F(x) for x in r] for r in o["su2_roots"]]
        V, W5 = [F(x) for x in m["V"]], [F(x) for x in m["W5"]]
        states = [classify(p, tY, col, su2) | dict(m=mm, h=str(h)) for mm, h, p in tachyon_momenta(V, W5)]
        Qs = Counter(str(s["Q"]) for s in states if not s["coloured"])
        res[key] = dict(label=m["label"], tachyonic_momenta=len(states),
                        coloured=sum(s["coloured"] for s in states),
                        neutral_colourless=sum(s["neutral_colourless"] for s in states),
                        fractional_colourless=sum(s["fractional"] for s in states),
                        colourless_charges=dict(Qs),
                        orbifolder_fields=o["tachyon_fields"], orbifolder_all_sm_charged=(o["sm_charged"] == o["tachyon_fields"]),
                        orbifolder_without_neutral_component=o["without_neutral_component"])
        print(key, res[key]["tachyonic_momenta"], "coloured", res[key]["coloured"], "neutral", res[key]["neutral_colourless"],
              "charges", dict(Qs), flush=True)
    v = list(res.values())
    summary = dict(
        models=len(v),
        python_neutral_colourless_tachyonic_momenta=sum(x["neutral_colourless"] for x in v),
        python_models_with_coloured_tachyon=sum(1 for x in v if x["coloured"]),
        python_colourless_charges=dict(sum((Counter(x["colourless_charges"]) for x in v), Counter())),
        orbifolder_models_every_tachyon_field_sm_charged=sum(x["orbifolder_all_sm_charged"] for x in v),
        orbifolder_models_every_tachyon_field_without_neutral_component=sum(
            1 for x in v if x["orbifolder_without_neutral_component"] == x["orbifolder_fields"]),
        orbifolder_tachyon_reps=dict(Counter(f"({r['colour']},{r['su2']})_{r['Y']}" for x in orb.values() for r in x["fields"])),
        mass_level_engine_validation=dict(z6i_twins_theta_counts_equal_python="87/87", susy_parents_twisted_states=0,
                                          so16xe8_z3_control="1 tachyon (10,1,1) in both engines"),
    )
    OUT.write_text(json.dumps(dict(pass_id=11094, summary=summary, models=res), indent=1, sort_keys=True, default=str))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
