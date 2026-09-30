# Pass 11204 — a second law by counting: arrow-free reversible dynamics becomes rare

Producer: `analysis/w33_pass11204_second_law_counting.py`
Certificate: `data/w33_pass11204_second_law_counting.json`
Regression: `tests/test_w33_pass11199_to_11206_frozen.py`

**Question.** A Clifford tick has no intrinsic arrow (A = 0, Pass 11183) iff some choice of subsystems makes it a
product of independent one-qutrit ticks. Equivalently, F₃^{2n} is an orthogonal sum of n invariant nondegenerate
planes. How common is arrow-free dynamics as the register grows?

**Method.**
* The decomposition is tested exactly for each sampled tick: invariant planes are found (cyclic span(v, Sv) or pairs
  of eigenvectors), and a backtracking search looks for n pairwise orthogonal ones.
* The test agrees element by element with the arrow computation on **all 51 840 elements of Sp(4,3)** (19 152
  arrow-free = 133/360 exactly). This is the known-positive control.
* Samples are uniform random elements of Sp(2n,3), built as random symplectic bases.

**Result.**

| n qutrits | arrow-free fraction | source |
|---|---|---|
| 1 | 1 | trivial |
| 2 | 133/360 = 0.3694 | exact (Pass 11188); sample 0.3715 ± 0.0034 |
| 3 | 55241/884520 = 0.0625 | exact (Pass 11188); sample 0.0600 ± 0.0017 |
| 4 | 0.0058 ± 0.0005 | 116 / 20 000 |
| 5 | < 0.001 (95%) | 0 / 3 000 |

The fraction falls by factors of about 2.7, 5.9, 11 and then more than 6, faster than exponentially.

**Reading.** On the substrate an arrow of time is generic. Reversible dynamics that some choice of subsystems renders
arrow-free is a vanishing corner once there are a few qutrits: a second law obtained by counting, not by coarse-graining
or a low-entropy initial state. The counting is over Clifford ticks with the uniform measure. No claim is made about
thermodynamic entropy or physical initial conditions.

**Cross-track note (added with Passes 11207–11208).** Arrow-free is A = 0, and by Theorem A = n − c(S) (Pass 11207)
that holds exactly when c = n, i.e. the tick preserves n mutually orthogonal nondegenerate planes. Summing over the
conjugacy classes turns the sampled rows into exact fractions:

| n | exact arrow-free fraction | source | sample above |
|---|---|---|---|
| 4 | 43 815 911 / 8 703 676 800 = 0.005034 | all 278 classes of PSp(8,3), `data/w33_pass11208_arrow_universality.json` | 0.0058 ± 0.0005 (1.5 σ) |
| 5 | 36 659 605 018 759 / 180 638 872 346 275 200 = 0.000203 | all 940 classes of PSp(10,3), same file | 0 / 3000 |

The successive ratios are 2.71, 5.91, 12.4 and 24.8, so the decay is super-exponential and the ratio roughly doubles
with each qutrit, as the counting here predicted. The mean number of protected qutrits converges instead: E[c] =
133/180, 106927/147420, 142080247/195832728 for n = 2, 3, 4, which tends to 0.7255… (Pass 11208).
