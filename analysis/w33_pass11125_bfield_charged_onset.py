#!/usr/bin/env python3
"""Pass 11125: where charged tachyons begin for the six neutral-exit survivors once the B-field is turned on.

Both Wilson-line tori at T = x + i y (x = B-field, 0 <= x <= 1/2 in steps of 0.1; y = 0.7 ... 1.4), family torus at
rho.  Each point is classified by the explicit Witten-twisted winding enumerator of Pass 11107 (level-matched states,
beta-projected) with the torus-resolved SM classification of Pass 11119: N = only SM-neutral tachyons, X = charged and
neutral, '.' = none.  Frozen map: data/w33_pass11125_bfield_map_six.json.

Results:
  * the map is the SAME in all six models (up to the two '.' cells): charged tachyons need BOTH a B-field and a small
    radius -- the charged region is {x >= 0.3, y <= 1.0} U {x >= 0.2, y <= 0.9} U {x >= 0.1, y <= 0.8}; at x = 0 the
    whole range down to y = 0.7 is neutral-only;
  * the charged onset lies at y <= 1.0-1.1 for every B in [0, 1/2], at least 0.4 below the neutral onset on the B = 0
    diagonal (y_c = 1.51 or 1.73, Pass 11123): the neutral instability comes first along every constant-B line;
  * at y = 1.4 a B-field x >= 0.4 removes the neutral tachyon altogether in four of the six (46043, 24165, 40521, 5904):
    these four are exactly the models with the golden critical radius phi^2/sqrt3 (Pass 11123), the two sqrt3 models keep
    it -- B shifts the onset, and the lower (golden) onset is the one B can push below 1.4.
Scope: equal T on the two Wilson-line tori; B in [0, 1/2] (B -> -B is a reflection of the lattice sum); grid spacing 0.1
in both directions, so the boundary is located to +-0.1.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11125_bfield_map_six.json"
OUT = ROOT / "data" / "w33_pass11125_bfield_charged_onset.json"
XS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
YS = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4]
CODE = {(): '.', ('neutral',): 'N', ('charged',): 'C', ('charged', 'neutral'): 'X'}


def grid(d):
    return {f"{x}|{y}": CODE[tuple(d[f"{x}|{y}"])] for x in XS for y in YS}


def summarize():
    D = json.loads(DATA.read_text())
    maps = {m: grid(d) for m, d in D.items()}
    onset = {}
    for m, g in maps.items():
        onset[m] = {str(x): max((y for y in YS if g[f"{x}|{y}"] in 'XC'), default=None) for x in XS}
    charged_cells = {m: sorted(k for k, v in g.items() if v in 'XC') for m, g in maps.items()}
    empty_cells = {m: sorted(k for k, v in g.items() if v == '.') for m, g in maps.items()}
    ref = next(iter(charged_cells.values()))
    res = dict(pass_id=11125, maps=maps, highest_charged_y=onset, charged_cells=charged_cells, empty_cells=empty_cells,
               same_charged_region_all_six=all(v == ref for v in charged_cells.values()),
               neutral_only_at_zero_B=all(all(g[f"0.0|{y}"] == 'N' for y in YS) for g in maps.values()),
               charged_onset_at_most_1p0=all(all(v is None or v <= 1.0 for v in o.values()) for o in onset.values()),
               lifted_by_B_at_1p4=sorted(m for m, e in empty_cells.items() if e))
    crit = json.loads((ROOT / "data" / "w33_pass11123_critical_point_amplitude.json").read_text())['critical_radii']
    res['lifted_are_golden'] = sorted(m for m, y in crit.items() if abs(y - 1.5115226281523415) < 1e-6) == res['lifted_by_B_at_1p4']
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    print(json.dumps({k: v for k, v in r.items() if k not in ('maps', 'charged_cells')}, indent=1))
