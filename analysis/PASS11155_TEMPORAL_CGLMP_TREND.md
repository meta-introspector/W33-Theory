# Pass 11155 — time climbs toward the algebraic maximum with dimension; space barely moves

Producer: `analysis/w33_pass11155_temporal_cglmp_trend.py`
Scan (d = 6, 7): `analysis/w33_pass11155_scan_cglmp_d6_7.py`
Frozen: `data/w33_pass11147_cglmp_optima_d2_5.json`, `data/w33_pass11155_cglmp_optima_d6_7.json`
Regression: `tests/test_w33_pass11155_temporal_cglmp_trend.py`

| d | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| space | 2.8284 | 2.9149 | 2.9727 | 3.0157 | 3.0497 | 3.0775 |
| time (one qudit, two times) | 2.8284 | 3.1628 | 3.2831 | 3.3773 | 3.4531 | 3.5168 |
| gap | 0 | 0.248 | 0.310 | 0.362 | 0.403 | 0.439 |

* The spatial optima reproduce the literature values at every d.
* The temporal deficit from the algebraic maximum, 4 − I_t(d), falls like **d^−0.65** (log–log fit over d = 3…7). The
  spatial deficit falls like d^−0.19. This is consistent with the temporal optimum approaching 4 as d → ∞, as
  Budroni–Emary found for Leggett–Garg.
* The d = 6 and 7 temporal values come from 8 restarts each, so they are lower bounds.

Reading: extra levels let one system "tell its future self" more through a projective measurement. Two separated
systems gain almost nothing from them. That is the dimensional face of Passes 11144 and 11147.
