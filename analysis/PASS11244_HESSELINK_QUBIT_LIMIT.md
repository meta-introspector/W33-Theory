# Pass 11244: the char-2 unipotent factor in closed form — the qubit limit 0.66516065…

Producer: `analysis/w33_pass11244_hesselink_qubit_limit.py`
Certificate: `data/w33_pass11244_hesselink_qubit_limit.json`
Regression: `tests/test_w33_pass11244_hesselink_qubit_limit.py`

**Starting point.** Pass 11233 reduced the qubit protection limit to F_m, the mean of m₁/2 + χ₂m₂ over the unipotent
elements of Sp(2m,2), and knew F_m only for m ≤ 6. Two identities, compared exactly with every unipotent class of
Sp(2m,2) for m ≤ 6, close the gap.

1. **Jordan types.** The total weight of the unipotent elements of each Jordan type λ equals the odd-q symplectic
   weight of Pass 11215 evaluated at q = 2, summed over the signs of the even parts. Checked on all **92** Jordan
   types, which are all the symplectic partitions with m ≤ 6.
2. **Hesselink index at 2.** Among the unipotents of type λ, the fraction with χ₂ = 0 is **2^{−m₂}** when m₂ is even
   and **0** when m₂ is odd, independent of the other parts. Checked on all **53** types with a part 2. For even m₂
   this is the chance that a random nondegenerate symmetric F₂-form of rank m₂ is alternating (Pass 11245 checks
   that count).

**Consequences.**
* F_m is now given for every m by a finite sum over partitions. It reproduces F₁…F₆ exactly and gives
  F₇ = 0.264518475, F₈ = 0.264515644, …, F_∞ ≈ 0.2645147.
* A Monte Carlo of F₇ (240 000 random elements, 2797 unipotent) gives 0.2639 ± 0.0090 (z = −0.07). That is consistent,
  but only a weak test.
* **Qubit limit, conditional on identities 1 and 2: lim E_n[c] = 0.66516065**, with bracket
  [0.6651604993, 0.6651614201] from F_m ≤ log₂3 on the tail m > 20.
* The unconditional bracket remains Pass 11233's [0.6626, 0.6777].
* **Correction.** Pass 11233's "F_m = F₆" estimate 0.66516078 was slightly high. The 0.66516061 quoted mid-session
  omitted the m > 22 tail; the value with the tail is **0.66516065**.
* Qubits protect less than qutrits: 0.665 against 0.7255 (Pass 11215).
