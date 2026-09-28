# Pass 11110 — the one-loop potential over the reachable moduli: no minimum inside the tachyon-free region

Producer: `analysis/w33_pass11110_moduli_map.py`
Data: `data/w33_pass11110_beta_integrands_model2.json` (new points); `data/w33_pass11106_beta_integrands_model2.json`
Certificate: `data/w33_pass11110_moduli_map.json`
Regression: `tests/test_w33_pass11110_moduli_map.py`

This pass computes model 2's potential with the full one-loop engine (Pass 11106). Every point is first certified
tachyon-free with the explicit enumerator of Pass 11107, because the modular integral diverges inside the tachyonic
region; model 2's T_WL = 1.5i point is dropped for that reason. The potential is Λ₄ = −½M⁴I in string frame at fixed
dilaton.

## Results (20 points)

| moduli (T_WL1, T_WL2; T\*) | Λ₄/M⁴ |
|---|---|
| (1.75i, 1.75i; ρ) — just above the critical radius √3 | **501.7 (lowest)** |
| (1.8i, 1.8i; ρ) | 531.7 |
| (0.5+1.8i, 0.5+1.8i; ρ) | 531.9 |
| (2i, 2i; ρ) | 657.6 |
| (0.3+2i, 0.3+2i; ρ) / (0.5+2i, 0.5+2i; ρ) | 657.7 / 657.7 |
| (2.6i, 1.8i; ρ) / (1.8i, 2.6i; ρ) | 759.7 / 769.1 |
| (2.2i, 2.2i; ρ) | 794.8 |
| (1.8i, 3i; ρ) | 884.4 |
| (2.5i, 2.5i; ρ) / (3i, 3i; ρ) | 1022.9 / 1464.9 |
| (2i, 2i; T\* = 0.25+i, i, 0.5+1.2i, 1.2i, 1.5i, 2i) | 676.5, 687.4, 703.1, 710.7, 800.0, 1006.1 |
| (2i, 2i; **T\* = 3.2i**, the hierarchy point of Pass 11109) | **1563.6** |
| (2i, 2i; T\* = 4i) | 1943.2 |

## What the map shows

* **Λ > 0 at every point.**
* **The family modulus has its minimum at T\* = ρ.** The point Im T\* ≈ 3.2 that the charm/top ratio needs is higher in
  energy by a factor of 2.4. So the vacuum does not choose the hierarchy (Pass 11109).
* **The Wilson-line moduli have no minimum.** Λ falls as either radius shrinks, separately or together. The B-field
  directions are almost flat, with a slight rise away from Re T = 0.
* **The lowest tachyon-free point sits on the boundary**, at Im T_WL = 1.75, just above the critical radius √3. The
  gradient points into the region where the fractionally charged winding state is tachyonic (Pass 11107).
* **The dilaton runs away.** In the 4D Einstein frame the one-loop potential is V_E = e^{4φ₄}Λ₄(T). With Λ₄ > 0
  everywhere on the tachyon-free region it has no stationary point in φ₄: the dilaton runs to weak coupling.

## Reading

Inside the region where they are tachyon-free, the survivors have no one-loop vacuum. The family modulus alone is
stabilised, and at the wrong place for the quark hierarchy. The Wilson-line moduli and the dilaton run: the former into
a charged tachyon, the latter to zero coupling.

Scope: one model (2) and 20 points of the moduli space, sampled on the directions the potential selects. This is not a
global minimisation. The WL-torus complex-structure moduli are fixed by the Z3 (U = ρ).
