# Pass 11122 — the six neutral-exit survivors fall steeply to the neutral onset under a volume law

Producer: `analysis/w33_pass11122_six_potentials.py`
Frozen integrands: `data/w33_pass11122_six_beta_integrands.json`
Certificate: `data/w33_pass11122_six_potentials.json`
Regression: `tests/test_w33_pass11122_six_potentials.py`

The Pass 11106 one-loop engine was run for the six models whose only nearby instability is Standard-Model-neutral
(36621, 46043, 24165, 40521, 5904, 17224). Setup:
* both Wilson-line tori at i·y with B = 0, and the family torus at ρ;
* Witten-twisted β sectors with K = 3 cosets, on the 864-point fundamental-domain grid, with the tail fit;
* three radii per model: y_c + 0.07 (just above the Pass 11119 onset), 2.0 and 2.5.

| model | y | Λ_β / M⁴ | Λ/y² |
|---|---|---|---|
| 36621 | 1.799 / 2.0 / 2.5 | 510.0 / 633.1 / 993.7 | 157.6 / 158.3 / 159.0 |
| 46043 | 1.583 / 2.0 / 2.5 | 396.0 / 635.8 / 995.3 | 158.0 / 158.9 / 159.2 |
| 24165 | 1.583 / 2.0 / 2.5 | 392.8 / 633.7 / 994.0 | 156.7 / 158.4 / 159.0 |
| 40521 | 1.583 / 2.0 / 2.5 | 392.8 / 633.7 / 994.0 | 156.7 / 158.4 / 159.0 |
| 5904 | 1.583 / 2.0 / 2.5 | 392.8 / 633.7 / 994.0 | 156.7 / 158.4 / 159.0 |
| 17224 | 1.799 / 2.0 / 2.5 | 507.5 / 631.1 / 992.4 | 156.8 / 157.8 / 158.8 |

Findings:
* **All six decrease monotonically toward the onset.** There is no barrier and no interior minimum along the diagonal,
  so the flow reaches the neutral boundary.
* **Volume law.** Λ_β = a y² + b with a = 160.0–160.8 and b = −4.8 to −12.6. The drive is the large-volume
  (n_B − n_F) term, the same in all six. The model-dependent part is a small *negative* correction, largest near the
  onset: the potential falls slightly faster than the volume law there, never flatter.
* **Steepness.** dΛ/dy ≈ 575–615 M⁴ per unit Im T near the onset. This is a string-scale force, not slow roll.

Scope:
* This is the B = 0 diagonal only, at three radii per model.
* The moduli-independent θ sectors are omitted: they shift Λ by a constant, and Pass 11106 has the total positive.
* 24165, 40521 and 5904 agree to 10⁻⁴. Their β-sector data coincide on this diagonal; this is not a claim that the
  models are equivalent.
