#!/usr/bin/env python3
"""Pass 11108: which SO(16)xSO(16) A8 models stay tachyon-free where the one-loop potential sends their moduli?

The one-loop potential shrinks the Wilson-line tori (Passes 11106, 11110), and below a critical radius a fractionally
charged winding state becomes tachyonic (Pass 11107).  Here the explicit tachyon enumerator (Pass 11107) is run on every
model, with both Wilson-line tori at T_WL = rho (the SU(3) point, B-field 1/2), i, 1.5i, 2i and the family torus at rho:

  A. the 104 SM-like models of Pass 11095 (frozen: data/w33_pass11108_tachyon_levels_104.json);
  B. a fresh Wilson-line scan with the non-SUSY orbifolder (8 x 50,000 tries, seeds 20260941-20260952, the 4 productive
     Witten-shift classes; WSL ~/orb/p1109x/a8/rescan), SM-like and tachyon-free at large radius, converted with
     nso2orb (the unused shift row dropped) and deduplicated by shift vectors (frozen:
     data/w33_pass11108_rescan_levels.json).

A result: tachyon-free at 2i: 104/104; at 1.5i and at i: 55/104; at rho: 1/104 (model 18).  The lowest level at rho
takes only two values, Delta = -1/18 (76 models) and -1/6 (27).  Model 18 fails the earlier tests: 48 fractionally
charged states stay massless (Pass 11097, all rule sets), its up Yukawas are untwisted (top = charm, Pass 11098) and it
has no down-type Yukawa.  None of the 12 survivors is tachyon-free at rho.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
LEV104 = ROOT / "data" / "w33_pass11108_tachyon_levels_104.json"
RESCAN = ROOT / "data" / "w33_pass11108_rescan_levels.json"
GAUNTLET = ROOT / "data" / "w33_pass11108_rescan_candidates_gauntlet.json"
YUK = ROOT / "data" / "w33_pass11098_yukawa_104.json"
OUT = ROOT / "data" / "w33_pass11108_tachyon_free_moduli_rescan.json"
SURVIVORS = ['2', '10', '13', '14', '15', '35', '53', '57', '69', '77', '78', '102']
RADII = (('rho', None), ('1.0', 1.0j), ('1.5', 1.5j), ('2.0', 2.0j))


def parse_scan(paths):
    """scan SM files -> {key: rows} (V0, V, W1..W6; the unused third row dropped), deduplicated"""
    import glob
    out = {}
    for p in sorted(x for pat in paths for x in glob.glob(pat)):
        t = open(p).read()
        for blk in t.split('begin model')[1:]:
            label = blk.split('Label:')[1].split()[0]
            rows = blk.split('Shifts and Wilsonlines:')[1].split('end model')[0].strip().splitlines()
            assert len(rows) == 9 and all(x.strip(',') in ('0/1', '0') for x in rows[2].split())
            rows = [r.strip() for r in rows[:2] + rows[3:]]
            key = '|'.join(' '.join(r.split()) for r in rows)
            out.setdefault(key, dict(label=label, rows=rows, source=Path(p).name))
    return list(out.values())


def levels(rows, rho_first=False):
    import w33_pass11107_winding_tachyons as T7
    res = {}
    for name, T in RADII:
        if rho_first and name != 'rho' and res.get('rho') is not None:
            res[name] = 'skipped'
            continue
        T = T7.RHO if T is None else T
        r = T7.tachyons('x', [T, T], T7.RHO, rows=rows)
        res[name] = round(r[0]['Delta'], 4) if r else None
    return res


def _job(m):
    return m['label'], dict(source=m['source'], rows=m['rows'], levels=levels(m['rows'], rho_first=True))


def rescan(paths):
    from multiprocessing import Pool
    models = parse_scan(paths)
    with Pool(8) as p:
        out = dict(p.map(_job, models))
    RESCAN.write_text(json.dumps(out, indent=1))
    return out


def summarize():
    lev = json.loads(LEV104.read_text())
    A = {k: sum(1 for v in lev.values() if v[k] is None) for k, _ in RADII}
    res = dict(pass_id=11108, A_tachyon_free_counts_104=A,
               A_free_at_rho=[m for m, v in lev.items() if v['rho'] is None],
               A_levels_at_rho=dict(Counter(str(v['rho']) for v in lev.values())),
               A_survivors_free_at_rho=[m for m in SURVIVORS if lev[m]['rho'] is None],
               A_survivors_free_on_axis_to_i=[m for m in SURVIVORS if lev[m]['1.0'] is None])
    if RESCAN.exists():
        rs = json.loads(RESCAN.read_text())
        res['B_new_models'] = len(rs)
        res['B_tachyon_free_at_rho'] = sum(1 for v in rs.values() if v['levels']['rho'] is None)
        res['B_free_at_rho'] = sorted(l for l, v in rs.items() if v['levels']['rho'] is None)
        res['B_levels_at_rho'] = dict(Counter(str(v['levels']['rho']) for v in rs.values()))
    yk = json.loads(YUK.read_text())
    twisted_up = [m for m, v in yk.items() if v['up']['untwisted'] == 0]
    res['A_twisted_up_models'] = len(twisted_up)
    res['A_twisted_up_free_at_rho'] = [m for m in twisted_up if lev[m]['rho'] is None]
    if GAUNTLET.exists():
        g = json.loads(GAUNTLET.read_text())
        res['B_candidates'] = {k: dict(Q_sectors=v['Q_l'],
                                       light_fractional_hidden_unbroken=v['fractional_light']['results']['A_hidden_unbroken|gauge+Z2W+PG+SG+R']['light_fractional_states'],
                                       light_fractional_hidden_broken=v['fractional_light']['results']['B_hidden_broken|gauge+Z2W+PG+SG+R']['light_fractional_states'],
                                       up_exponents=[h['exps'] for h in v['up']], down_tree_level=len(v['down']))
                               for k, v in g.items()}
        res['rho_safe_models_all_untwisted_quarks'] = all(v['Q_sectors'] == [0, 0, 0] for v in res['B_candidates'].values())
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    if len(sys.argv) > 1:
        rescan(sys.argv[1:])
    summarize()
