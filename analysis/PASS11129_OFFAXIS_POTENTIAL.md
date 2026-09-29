# Pass 11129 — the B-field is a near-flat direction: B = 0 is a shallow minimum just above the onset

Producer: `analysis/w33_pass11129_offaxis_potential.py`
Scan: `analysis/w33_pass11129_scan_offaxis.py`
Frozen: `data/w33_pass11129_offaxis_beta_integrands.json`
Regression: `tests/test_w33_pass11129_offaxis_potential.py`

The Pass 11122 engine was rerun with both Wilson-line tori at x + iy (x = the B-field) for one √3 model (36621) and one
golden model (40521). The x = 0 values come from 11122 on the same grid and fit, so common quadrature errors cancel in the
differences.

| model | y | Λ(0.25) − Λ(0) | Λ(0.5) − Λ(0) |
|---|---|---|---|
| 36621 | 1.799 | +0.006 | +0.011 |
| 36621 | 2.0 | −0.002 | −0.004 |
| 40521 | 1.583 | +0.015 | +0.027 |
| 40521 | 2.0 | −0.003 | −0.006 |

* **Shape.** Λ(x) − Λ(0) = A(1 − cos 2πx)/2: the x = 0.5 difference is about twice the x = 0.25 one in all four cases.
* **Just above the onset, B = 0 is a minimum** (A > 0). The flow reaches the neutral exit, which is a B = 0 feature
  (Pass 11125).
* **At y = 2, B = 0 is a slight maximum** (A < 0).
* **B is effectively frozen.** |A| is at most 5·10⁻⁵ of the radial force dΛ/dy. The B-dependence comes only from
  winding states and is exponentially suppressed in Im T, so the potential does not select B during the fall.

Scope:
* Two models, two radii, three values of x; the reflection x → −x covers the rest.
* A is about 2·10⁻⁵ of Λ. Its sign is trustworthy only because the grid and fit are shared.
* No kinetic analysis.
