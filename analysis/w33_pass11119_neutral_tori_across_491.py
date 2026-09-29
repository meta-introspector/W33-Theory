#!/usr/bin/env python3
"""Pass 11119: models whose winding tachyons separate onto their own Wilson-line torus -- across all 491 SM-like models,
and which of them survive the full gauntlet.

For every model the two Wilson-line tori are shrunk one at a time to Im T = 1 (the other at 3i, family torus at rho) and
the tachyons that appear are classified as SM-neutral (Y = 0, colour and SU(2) singlets) or charged (explicit enumerator,
Pass 11107; frozen: data/w33_pass11119_torus_resolved_tachyons_491.json).  A torus is a 'neutral torus' if every tachyon
it produces is SM-neutral.  'neutral-only' models: one neutral torus, the other tachyon-free or neutral (any tachyon met
by shrinking a single torus preserves the SM); 'separated': one neutral and one charged torus (like 53, 57).
The new-model candidates are put through the gauntlet of the survivors (field dumps of both engines; frozen:
data/w33_pass11119_gauntlet_new_candidates.json): (a) all fractionally charged fermions massive (hidden group unbroken,
all rules; Pass 11097), (b) a single heavy top from a localized Higgs, (c) a single heavy bottom (Passes 11101-11102).
Results: 60/491 models place their neutral tachyons on a torus of their own (38 neutral-only, 22 separated), ALL with
twisted quarks; 19 new models pass the gauntlet, giving 21 with 53 and 57.  Along the B = 0 diagonal (both Wilson-line
tori equal; data/w33_pass11119_diagonal_onset_21.json) 8 of the 21 meet ONLY neutral tachyons, just below their onset
(Im T = 1.513 or 1.729) and at Im T = 1; 6 of them stay neutral-only down to Im T = 0.6, charged tachyons appearing only
with a B-field at small radius (0.2 + 0.9i) and at the SU(3) point rho (data/w33_pass11119_extended_moduli_8.json).
In 2 of the 21 (3818, 5273) the first diagonal tachyon is charged.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TC = ROOT / "data" / "w33_pass11119_torus_resolved_tachyons_491.json"
GA = ROOT / "data" / "w33_pass11119_gauntlet_new_candidates.json"
YK = ROOT / "data" / "w33_pass11098_yukawa_104.json"
S104 = ROOT / "data" / "w33_pass11108_model_shifts_104.json"
SECT = ROOT / "data" / "w33_pass11113_quark_sectors_387.json"
FRAC = ROOT / "data" / "w33_pass11097_symmetry_test_104.json"
OUT = ROOT / "data" / "w33_pass11119_neutral_tori_across_491.json"
DIAG = ROOT / "data" / "w33_pass11119_diagonal_onset_21.json"
EXT = ROOT / "data" / "w33_pass11119_extended_moduli_8.json"
KEY = 'A_hidden_unbroken|gauge+Z2W+PG+SG+R'


def classify(v):
    k1, k2 = v['torus1']['kinds'], v['torus2']['kinds']
    if (k1 == ['neutral'] and k2 in ([], ['neutral'])) or (k2 == ['neutral'] and k1 == []):
        return 'neutral-only'
    if (k1 == ['neutral'] and k2 == ['charged']) or (k1 == ['charged'] and k2 == ['neutral']):
        return 'separated'
    return 'other'


def summarize():
    tc = json.loads(TC.read_text())
    lab104 = {v['label']: k for k, v in json.loads(S104.read_text())['models'].items()}
    yk, sect, frac = json.loads(YK.read_text()), json.loads(SECT.read_text()), json.loads(FRAC.read_text())
    pattern = Counter((tuple(v['torus1']['kinds']), tuple(v['torus2']['kinds'])) for v in tc.values())
    rows = {}
    for lab, v in tc.items():
        cls = classify(v)
        if cls == 'other':
            continue
        if lab in lab104:
            m = lab104[lab]
            rows[lab] = dict(cls=cls, model=m, quarks='twisted' if yk[m]['up']['untwisted'] == 0 else 'untwisted',
                             light_fractional=frac[m]['results'][KEY]['light_fractional_states'])
        else:
            s = sect[lab]
            rows[lab] = dict(cls=cls, model=lab, quarks='twisted' if s['up_triples'] and s['up_untwisted'] == 0 else 'untwisted')
    ga = json.loads(GA.read_text()) if GA.exists() else {}
    for lab, g in ga.items():
        rows[lab].update(light_fractional=g['light_fractional'], single_heavy_top=g['single_heavy_top'],
                         single_heavy_bottom=g['single_heavy_bottom'])
    passing = sorted(l for l, r in rows.items() if r.get('light_fractional') == 0 and r['quarks'] == 'twisted'
                     and (r.get('single_heavy_top', r['model'] in ('53', '57')) and r.get('single_heavy_bottom', r['model'] in ('53', '57'))))
    res = dict(pass_id=11119, models=len(tc), pattern={f"{k[0]}|{k[1]}": v for k, v in pattern.items()},
               counts=dict(Counter((r['cls'], r['quarks']) for r in rows.values())).__repr__(),
               neutral_only=sum(1 for r in rows.values() if r['cls'] == 'neutral-only'),
               separated=sum(1 for r in rows.values() if r['cls'] == 'separated'),
               all_twisted=all(r['quarks'] == 'twisted' for r in rows.values()),
               passing_gauntlet=passing,
               passing_neutral_only=[l for l in passing if rows[l]['cls'] == 'neutral-only'],
               rows=rows)
    if DIAG.exists():
        dg = json.loads(DIAG.read_text())
        res['diagonal_neutral_only'] = sorted(l for l, v in dg.items() if v['crit'] and set(v['just_below']) == {'neutral'} and set(v['at_1']) == {'neutral'})
        res['diagonal_charged_first'] = sorted(l for l, v in dg.items() if v['crit'] and set(v['just_below']) == {'charged'})
        res['diagonal_crit'] = {l: v['crit'] for l, v in dg.items()}
    if EXT.exists():
        ex = json.loads(EXT.read_text())
        res['extended_charged_points'] = {l: sorted(p for p, c in d.items() if 'charged' in c) for l, d in ex.items()}
        res['neutral_down_to_0p6_on_axis'] = sorted(l for l, d in ex.items() if 'charged' not in d['0.8'] and 'charged' not in d['0.6'])
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True, default=str))
    print(json.dumps({k: v for k, v in res.items() if k != 'rows'}, indent=1))
    return res


if __name__ == "__main__":
    summarize()
