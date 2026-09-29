#!/usr/bin/env python3
"""Pass 11155: the temporal CGLMP optimum climbs toward the algebraic maximum 4 with dimension, the spatial one barely
moves -- d = 2..7.

Data: data/w33_pass11147_cglmp_optima_d2_5.json (d = 2..5) and data/w33_pass11155_cglmp_optima_d6_7.json (d = 6, 7; 8
restarts each, so the d = 6, 7 values are lower bounds on the optima).
      d     spatial     temporal    4 - temporal
      2     2.82843     2.82843     1.172
      3     2.91485     3.16281     0.837
      4     2.97270     3.28306     0.717
      5     3.01571     3.37730     0.623
      6     3.04970     3.45309     0.547
      7     3.07754     3.51682     0.483
  * the spatial optima reproduce the literature values (d = 6: 3.0497, d = 7: 3.0776);
  * the temporal deficit 4 - I_t(d) falls like d^(-0.65) over d = 3..7 (log-log fit) -- consistent with the temporal optimum
    approaching the algebraic maximum 4 as d -> infinity (cf. Budroni-Emary for Leggett-Garg); the spatial deficit falls
    far more slowly;
  * the gap temporal - spatial grows monotonically: 0, 0.248, 0.310, 0.362, 0.403, 0.439.
Reading: the more levels a single system has, the more it can 'tell its future self' through a projective measurement,
while two separated systems gain almost nothing from extra levels.  Scope: projective Lueders measurements, two settings.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / "data" / "w33_pass11147_cglmp_optima_d2_5.json"
B = ROOT / "data" / "w33_pass11155_cglmp_optima_d6_7.json"
OUT = ROOT / "data" / "w33_pass11155_temporal_cglmp_trend.json"


def summarize():
    d = {**json.loads(A.read_text()), **json.loads(B.read_text())}
    rows = {k: dict(spatial=d[f"{k}|spatial"]['value'], temporal=d[f"{k}|temporal"]['value']) for k in range(2, 8)}
    ds = np.arange(3, 8)
    tdef = np.array([4 - rows[k]['temporal'] for k in ds])
    sdef = np.array([4 - rows[k]['spatial'] for k in ds])
    slope_t = float(np.polyfit(np.log(ds), np.log(tdef), 1)[0])
    slope_s = float(np.polyfit(np.log(ds), np.log(sdef), 1)[0])
    gaps = [rows[k]['temporal'] - rows[k]['spatial'] for k in range(2, 8)]
    res = dict(pass_id=11155, table={str(k): v for k, v in rows.items()}, temporal_deficit_loglog_slope=slope_t,
               spatial_deficit_loglog_slope=slope_s, gaps=gaps, gaps_monotone=all(a < b for a, b in zip(gaps, gaps[1:])))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
