# Pass 11120 — toward a closed-form ρ-tachyon criterion: a sharp necessary condition, and why no norm invariant suffices

Producer: `analysis/w33_pass11120_rho_cell_criterion.py`
Data: `data/w33_pass11120_rho_cell_invariants_491.json`
Certificate: `data/w33_pass11120_rho_cell_criterion.json`
Regression: `tests/test_w33_pass11120_rho_cell_criterion.py`

## What is computed

The data cover all 491 models (the 104 of Pass 11095 plus the 387 of the rescan). A Wilson-line cell is a pair (N, κ):
* N ∈ Z3² is the winding class of the two Wilson-line tori;
* κ = (π + V0)·W mod 1.

For every cell, the minimal norm ℓ²_min of the β-projected gauge coset E8² + V0 + N·W is computed. Only ℓ² < 2 can carry
a tachyon, since Δ = −½ + p_R²/2 and ℓ² + p_L² − p_R² = 1.

## Results

Patterns of (overall minimum, minimum over single-torus N, minimum over mixed N), over N ≠ 0:

| pattern | tachyonic at ρ | ρ-safe |
|---|---|---|
| (1, 1, 1) | 185 | 0 |
| (1, 1, 5/3) | 95 | 0 |
| (1, 5/3, 1) | 115 | **5** |
| (5/3, 5/3, 5/3) | 91 | 0 |

* **Sharper necessary condition.** Every ρ-safe model has exactly **2 cells of norm 1** and **10 of norm 5/3**, the
  sparsest structure present. Both norm-1 cells have mixed N, and κ is nonzero in one torus only.
* **Not sufficient.** 15 tachyonic models share that sparse structure. Model 10968 even has the same norm-1 cells as the
  safe model 8666.
* **The universal ρ levels.** Δ = −1/6 corresponds to p_R² = 2/3 and Δ = −1/18 to p_R² = 8/9. These are Narain vectors of
  the Wilson-line torus at the SU(3) point with the order-3 offsets.

## Reading

The gauge lattice gives a sharp necessary condition but not a criterion. The decision also needs:
* the torus offsets, o_t(N, κ) = κ_t/3 + V0·W_t + W_t·(N·W)/2;
* the torus norms at ρ.

That is the finite cell × torus lookup the explicit enumerator of Pass 11107 performs. No single gauge-lattice norm
invariant separates the ρ-safe models, so the "closed form" is that finite lookup and not a norm inequality.
