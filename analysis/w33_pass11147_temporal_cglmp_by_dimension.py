#!/usr/bin/env python3
"""Pass 11147: the temporal advantage grows with dimension -- one qudit measured at two times beats two qudits in space
for every d >= 3, and the gap widens (d = 2 is the only coincidence).

Scan: analysis/w33_pass11147_scan_cglmp_by_dimension.py (CGLMP I_d, two settings, d outcomes; temporal: pure start,
projective Lueders measurements, P(a,b|x,y) = <a_x|rho|a_x> |<b_y|a_x>|^2; spatial: pure bipartite state, projective
measurements; L-BFGS-B, 24 restarts for d <= 4, 16 for d = 5).  Frozen (with the optimal parameters):
data/w33_pass11147_cglmp_optima_d2_5.json.
Results:
      d   spatial optimum        temporal optimum    gap
      2   2.8284271 (2 sqrt2)    2.8284271          0          (Fritz: coincide)
      3   2.9148542 (1+sqrt(11/3)) 3.1628065        0.2480
      4   2.9726983              3.2830617          0.3104
      5   3.0157105              3.3773031          0.3616
  * validation: the spatial optima reproduce the known values (d = 3: exactly 1 + sqrt(11/3), recovered by an integer-
    relation search; d = 4, 5: Acin et al. / Chen et al.), and d = 2 reproduces Fritz's temporal = spatial;
  * no low-height algebraic form was found for the temporal optima (degree <= 6, coefficients <= 2000: only spurious
    high-coefficient fits at 13-digit precision) -- they are reported numerically;
  * the gap grows monotonically with d (0, 0.248, 0.310, 0.362), in line with Budroni-Emary's finding that temporal
    (Leggett-Garg-type) quantum values grow toward the algebraic maximum with dimension; the algebraic maximum of I_d is 4.
Reading: for qubits space and time are interchangeable (Fritz); for every larger dimension one system measured twice
is 'more correlated' than any pair -- because the measurement can hand information forward in time, which space forbids.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11147_cglmp_optima_d2_5.json"
OUT = ROOT / "data" / "w33_pass11147_temporal_cglmp_by_dimension.json"


def summarize():
    d = json.loads(DATA.read_text())
    rows = {}
    for dim in (2, 3, 4, 5):
        s, t = d[f"{dim}|spatial"]['value'], d[f"{dim}|temporal"]['value']
        rows[str(dim)] = dict(spatial=s, temporal=t, gap=t - s)
    gaps = [rows[str(k)]['gap'] for k in (2, 3, 4, 5)]
    res = dict(pass_id=11147, table=rows, d3_spatial_closed_form=1 + (11 / 3) ** 0.5,
               gap_monotone=all(a < b for a, b in zip(gaps, gaps[1:])), qubit_coincidence=abs(gaps[0]) < 1e-9)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
