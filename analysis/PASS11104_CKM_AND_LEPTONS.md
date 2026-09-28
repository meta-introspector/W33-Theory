# Pass 11104 — CKM texture and charged leptons of the 12 survivors

Producer: `analysis/w33_pass11104_ckm_and_leptons.py` (needs the field dumps; frozen scan
`data/w33_pass11104_ckm_leptons_12.json`)
Certificate: `data/w33_pass11104_ckm_and_leptons.json`
Regression: `tests/test_w33_pass11104_ckm_and_leptons.py`

## Method

Each model has three localized light-Higgs candidates for each of H_u and H_d, one at each fixed point of the family
torus. That gives 9 pairs (H_u, H_d) per model. For every pair:
* the up and down Yukawa matrices are built to all orders with the Pass 11102 engine;
* the CKM matrix is V = U_u† U_d, and the exponents of |V_ij| in ε are read off;
* the charged leptons use the same H_d;
* the neutrino Dirac couplings L N H_u use the same H_u.

## Results (12 models, 108 pairs)

| | result |
|---|---|
| top and bottom in the same doublet | **36 / 108**: only when H_u and H_d sit at the same fixed point, and each H_u has exactly one such partner; the other 72 give V_tb ~ ε², so they are excluded |
| CKM, aligned pairs | V_tb = O(1); \|V_cb\|, \|V_ub\|, \|V_ts\|, \|V_td\| ~ ε²; Cabibbo block O(1) (24/36) or ~ε (12/36) |
| charged leptons, same H_d | a single heavy τ in 36/36; **μ and e at the same order** ((0,1,1) or (0,2,2)) |
| neutrinos | 39 or 51 SM-singlet fermions; tree-level Dirac couplings in 12/12 |

## Comparison with data

* In the large-area regime m_c/m_t ~ ε². The texture then gives **V_cb ~ m_c/m_t**, but the data give
  V_cb/(m_c/m_t) ≈ 11 at M_Z. That is an order-of-magnitude tension that O(1) coefficients must absorb.
* The texture gives V_ub ~ V_cb, but the data give V_ub/V_cb ≈ 0.09. That is another factor of 10.
* The Δ(54) degeneracy of Passes 11102 and 11103 reaches the leptons: m_μ/m_e ≈ 207 is not generated.
* Light neutrinos need a seesaw: Majorana masses for N from the scalar VEVs, not computed here.

## Reading

The survivors impose one sharp requirement on the vacuum: the light H_u and H_d must be localized at the **same** fixed
point. Given that, the third generation separates correctly (V_tb ≈ 1, a single heavy t, b and τ). The two light
generations stay degenerate in every sector, and that degeneracy is the common obstruction. The texture is qualitatively
right in the heavy sector and off by about 10 in V_cb and V_ub.

Scope: parametric exponents with random O(1) coefficients. The vector-like states are treated only through the
rectangular Yukawa blocks. The choice of light Higgs pair is assumed, not derived.
