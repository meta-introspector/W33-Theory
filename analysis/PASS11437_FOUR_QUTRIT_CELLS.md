# Pass 11437 — four qutrits cell by cell: the n = 3 excess over 1/8 lives in one cell, and it does not persist

Producer: `analysis/w33_pass11437_four_qutrit_cells.py`
Certificate: `data/w33_pass11437_four_qutrit_cells.json`
Regression: `tests/test_w33_pass11433_11437.py`

**Method.** The one-gate bad-class fraction is an exact average over the magic-axis cells:
P_n = Σ_cells share × (bad fraction in the cell).
* The shares come from vector counts among the 3^{2n} − 1 images Mz₁:
  * collinear vectors (other than ±z₁): 3^{2n−1} − 3;
  * non-collinear vectors: 3^{2n} − 3^{2n−1}.
* Each cell is sampled directly, M = h·g:
  * g is a long random word;
  * h is a transvection map sending gz₁ to a uniformly random vector of the target cell (Pass 11352);
  * each class is decided by Pass 11421's sparse decider.
* **Control at n = 3** against the exact cells of Pass 11373 (1500 classes per cell type):

| | sampled | exact |
|---|---|---|
| collinear cells together | 0.152 | 65/432 = 0.1505 |
| non-collinear | 0.1227 | 10/81 = 0.1235 |
| P₃ | 0.1334 ± 0.0064 | 437/3276 = 0.13339 |

  The split of the collinear sample between "different lines" and "same line" is 1017 : 474. The exact vector ratio
  is 162 : 76.

## Four qutrits (3000 classes per cell type)

| cell | n = 3 (exact) | **n = 4 (sampled)** |
|---|---|---|
| same line, other | 77/342 = 0.225 | **140/1001 = 0.140** |
| collinear, different lines | 1/9 = 0.111 | **250/1997 = 0.125** |
| all collinear cells | 65/432 = 0.1505 | **0.130 ± 0.006** |
| non-collinear | 10/81 = 0.1235 | **0.127 ± 0.006** |
| **P_n** | **0.13339** | **0.1284 ± 0.0045** |

The shares are exact: collinear 273/820 and non-collinear 2187/3280. Two classes were undecided.

**Reading.**
* **Where the excess sits.** The n = 3 excess of P₃ over 1/8 comes from the **"same line, other"** cell, at 0.225 for
  n = 3. That cell falls to 0.140 at n = 4. The collinear-different-lines and non-collinear cells both sit near 1/8
  at n = 4.
* **The size of the drop.** This estimate and Pass 11421's uniform sample (0.118 ± 0.006) agree within 1.4σ. Their
  combination, about 0.125 ± 0.004, is consistent with 1/8.
* **Not a monotone trend.** The apparent "drop below 1/8" of Pass 11421 is mostly sampling noise around 1/8.
* **A conjecture, not a result:** the cell fractions approach 1/8 and P_n → 1/8. n = 3 is the outlier because of the
  same-line cell.
* **Refuted by Pass 11460.** At five qutrits P₅ = 0.1382 ± 0.0041, 3.3σ above 1/8, with every cell near 0.14. The
  sequence 1/8, 1/8, 0.1334, ≈ 0.125, 0.138 does not settle at 1/8.
