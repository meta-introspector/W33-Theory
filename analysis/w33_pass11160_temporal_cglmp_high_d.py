#!/usr/bin/env python3
"""Pass 11160: temporal vs spatial CGLMP up to d = 10 -- the temporal optimum keeps climbing toward 4.

Scan: analysis/w33_pass11160_scan_cglmp_d8_10.py (batched evaluator: unitaries from batched Hermitian
eigendecomposition, central-difference gradients in one batch per step; method check: reproduces the committed d = 3
optima 3.16280651 / 2.91485422 and the d = 4 spatial 2.97269827 from their frozen parameters).  5 restarts per case, so
the d = 8..10 temporal values are lower bounds.  Frozen: data/w33_pass11160_cglmp_optima_d8_10.json.
      d     spatial     temporal    4 - temporal
      8     3.10128     3.57620     0.424
      9     3.12168     3.62078     0.379
     10     3.13959     3.65763     0.342
(d = 2..7 in Passes 11147 and 11155.)  Over d = 3..10 the temporal deficit falls like d^-0.75 (log-log fit), the spatial
like d^-0.2; the gap temporal - spatial grows to 0.518 at d = 10.  Consistent with the temporal optimum approaching the
algebraic maximum 4 (Budroni-Emary for Leggett-Garg); not a proof.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FILES = [ROOT / "data" / "w33_pass11147_cglmp_optima_d2_5.json", ROOT / "data" / "w33_pass11155_cglmp_optima_d6_7.json"]
NEW = ROOT / "data" / "w33_pass11160_cglmp_optima_d8_10.json"
OUT = ROOT / "data" / "w33_pass11160_temporal_cglmp_high_d.json"


def summarize():
    rows = {}
    for f in FILES:
        for k, v in json.loads(f.read_text()).items():
            d, kind = k.split('|'); rows.setdefault(int(d), {})[kind] = v['value']
    for k, v in json.loads(NEW.read_text()).items():
        d, kind = k.split('|'); rows.setdefault(int(d), {})[kind] = max(v, rows.get(int(d), {}).get(kind, -9))
    ds = np.array([d for d in sorted(rows) if d >= 3])
    tdef = np.array([4 - rows[d]['temporal'] for d in ds]); sdef = np.array([4 - rows[d]['spatial'] for d in ds])
    st = float(np.polyfit(np.log(ds), np.log(tdef), 1)[0]); ss = float(np.polyfit(np.log(ds), np.log(sdef), 1)[0])
    gaps = {str(d): rows[d]['temporal'] - rows[d]['spatial'] for d in sorted(rows)}
    res = dict(pass_id=11160, table={str(d): rows[d] for d in sorted(rows)}, temporal_deficit_slope=st, spatial_deficit_slope=ss,
               gaps=gaps, gaps_monotone=all(gaps[str(a)] < gaps[str(b)] for a, b in zip(sorted(rows), sorted(rows)[1:])),
               method_check_d3=json.loads(NEW.read_text())['3|temporal'])
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
