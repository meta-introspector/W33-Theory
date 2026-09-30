# Pass 11191 — polygamy: prior art, and the qutrit bulge is not special

Producer: `analysis/w33_pass11191_polygamy_across_dimension.py`
Certificate: `data/w33_pass11191_polygamy_across_dimension.json`
Regression: `tests/test_w33_pass11191_polygamy_across_dimension.py`

**Prior art (found by the Pass 11188–11192 sweep; it corrects Passes 11148–11187).**
* The optimum ψ* = √(7/15)|000⟩ + √(2/15) Σ_{a=1,2} |a⟩(|0a⟩ + |a0⟩) lies in the family d|000⟩ + a Σ|j0j⟩ + b Σ|jj0⟩
  of **Allen & Meyer, "Polynomial monogamy relations for entanglement negativity", PRL 118, 080402 (2017),
  arXiv:1502.04807**. They block-diagonalise its partial transposes, give its marginal spectra, and **conjecture from
  numerics that the family traces the whole achievable negativity region for D > 2** (proved for qubits).
* The "global optimality open" of Passes 11162, 11167 and 11187 is that conjecture specialised to the sum, not a new
  open problem. Pass 11172's frontier claim is its (N_AB, N_AC) projection.
* **He & Vidal, PRA 91, 012339 (2015)** show numerically that negativity violates linear monogamy and conjecture the
  squared version. ψ* obeys the squared inequality: 8/15 ≤ N²_{A|BC} ≈ 0.945.
* **Not found in the literature:** the value 4/√15, the certified local maximum (11162), the unconditional bound 25/18
  (11167) and the torus-class branch and bound (11187). Those remain the repo's contributions. The branch and bound
  used floating point with a 10⁻⁹ margin, not directed rounding.

**Computation.** The negativity is normalised so that a maximally entangled pair has N = 1:
N = (‖ρ^T‖₁ − 1)/(D − 1), the repo's N for D = 3. The maximum of N_AB + N_AC on the Allen–Meyer family:

| D | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| max | 1.044815 | 1.032796 | 1.025909 | 1.021430 | 1.018277 | 1.015937 |

* D = 2 reproduces (8√2 − 4)/7, and a full-state-space search with 40 random restarts reaches the same value.
* D = 3 reproduces 4/√15.
* The excess over 1 decreases with D, roughly like 1/D. This is consistent with Allen–Meyer's remark that linear
  monogamy is recovered up to O(1/D).

**Correction (failure mode 2, over-read).** "The sum exceeds 1" is **not a qutrit effect**. Three qubits exceed 1 by
more (0.0448) than three qutrits (0.0328). Earlier wording that presented the qutrit bulge as special, or as a feature
of W(3,3), is withdrawn. The qutrit statements that stand are the exact value on the family, its certification and
the partial-class proofs.
