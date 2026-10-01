# Pass 11233: the qubit protection limit from the char-2 cycle index

Producer: `analysis/w33_pass11233_qubit_limit_cycle_index.py`
Certificate: `data/w33_pass11233_qubit_limit_cycle_index.json`
Regression: `tests/test_w33_pass11233_qubit_limit_cycle_index.py`

**Question.** Pass 11220 computed E_n[c], the mean number of qubits a random Clifford tick protects, exactly for
n ≤ 6 from GAP's conjugacy classes. It then extrapolated geometrically to lim E_n[c] ≈ 0.665166. Master's Pass 11215
had obtained the qutrit limit from Fulman's cycle index. For qubits the unipotent factor needs Hesselink's
characteristic-2 classes.

## What is exact

* **Additivity.** c = m₁(1)/2 + χ₂ m₂(1) + m₁(x²+x+1) (Pass 11217) is additive over primary parts. So
  Σ_n uⁿ E_n[c] = (1/(1−u)) · (F_unip/P_unip + F_U/P_U).
* **Unitary part (x²+x+1, Q = 2).** Pass 11215's unitary sums with q = 2 pass Steinberg's check. Its limiting
  contribution is **0.27261186563783750005…**, with partial sums stable to 40 digits.
* **Unipotent normalisation.** P(u) = Σ_m w_m u^m with w_m = 2^{2m²}/|Sp(2m,2)| (Steinberg). By the q-binomial theorem,
  P(u) = ∏_{i≥0} 1/(1 − u·2^{−1−2i}). So 1/P is entire, and the unipotent contribution to the limit is exactly
  F(1)/P(1), where F_m is the mean of m₁/2 + χ₂m₂ over the unipotent elements of Sp(2m,2).
* **F₁ … F₆ exactly.** From the class data:
  1, 61/128, 0.286720…, 0.265019…, 0.264576…, **0.2645298…**.
* **Validation.** Steinberg's counts hold for m ≤ 6. E[2^{dim ker(u−1)}] = 3 − 2^{1−2m}, proved by counting unipotents
  in a vector stabiliser and checked on the classes. The factors reproduce **all five exact values of Pass 11220**
  (n = 2…6) identically.

## The limit, and a correction

| | value |
|---|---|
| unitary contribution | 0.27261186563783750 |
| unipotent contribution through m = 6 | 0.39003115463 |
| weight of m ≥ 7 in the unipotent factor | 0.0095179 |
| rigorous bracket (F_m ≤ log₂3 for m ≥ 7) | [0.66264, 0.67773] |
| point estimate with F_m = F₆ for m ≥ 7 | **0.66516078** |

* **Correction to Pass 11220.** E_n[c] is **not monotone**. The unitary series has a negative term at n = 7
  (b₇ = −4.8·10⁻⁷), so E₇ < E₆ = 0.6651612. The limit therefore lies **below** E₆. Pass 11220's geometric extrapolation
  (0.665166, from three positive differences) was wrong in the sixth decimal.
* **Sharpened in Pass 11244.** There, closed forms for the char-2 unipotent factor (verified on all classes with
  m ≤ 6) give F_m for every m and the limit **0.66516061457…**, conditional on those identities.
* **Not completed.** The Monte Carlo run of F₇ and F₈ (transvection walks) was stopped in favour of Pass 11244, which
  carries the m = 7 Monte Carlo test.

**Prior art.**
* Master's Pass 11215 supplies the unitary sums and the odd-q method; this pass is the char-2 analogue for the unitary
  and normalisation parts.
* The other track's Passes 11234–11235 (pushed meanwhile) confirm the n = 2, 3 qubit means 29/40 and 5981/8960 by
  brute force over all of Sp(4,2) and Sp(6,2). Those are the first two values reproduced here.
