# Pass 11160 — temporal vs spatial CGLMP up to d = 10: time keeps climbing toward 4

Producer: `analysis/w33_pass11160_temporal_cglmp_high_d.py`
Scan: `analysis/w33_pass11160_scan_cglmp_d8_10.py` (batched evaluator; the method check reproduces the frozen d = 3, 4 optima)
Frozen: `data/w33_pass11160_cglmp_optima_d8_10.json`
Regression: `tests/test_w33_pass11160_temporal_cglmp_high_d.py`

| d | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|
| space | 2.828 | 2.915 | 2.973 | 3.016 | 3.050 | 3.078 | 3.101 | 3.122 | 3.140 |
| time | 2.828 | 3.163 | 3.283 | 3.377 | 3.453 | 3.517 | 3.576 | 3.621 | 3.658 |

* Over d = 3…10 the temporal deficit 4 − I_t falls like **d^−0.75**; the spatial deficit falls like d^−0.19.
* The gap grows monotonically, to 0.518 at d = 10.
* The d = 8–10 values used 5 restarts each, so the temporal values are lower bounds.
* This is consistent with the temporal optimum approaching the algebraic maximum 4, as Budroni–Emary found for
  Leggett–Garg. It is not a proof.

**Correction (Pass 11171).** The d^−0.75 deficit fit (d ≤ 10) is superseded. Data to d = 24 give 4 − I_d ≈ 3.4/d
(exponent 0.97 for d ≥ 6), and I_d < 4 is proved for every finite d. The d = 8, 9, 10 values here are reproduced.
