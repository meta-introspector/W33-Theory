# Pass 11266 — how often does magic break the arrow? Exact counts, and the one-qutrit rule

Producer: `analysis/w33_pass11266_one_magic_gate.py` (decisions by Pass 11252's exact criterion)
Certificate: `data/w33_pass11266_one_magic_gate.json`
Regression: `tests/test_w33_pass11265_11269.py`

## One qutrit, exhaustive

| cubic gates k | words | violating | probability |
|---|---|---|---|
| 1 | 216 | 18 | **1/12** |
| 2 | 216² = 46 656 | 18 954 | **13/32** |
| 3 | 216³ = 10 077 696 | 5 524 362 | **421/768** |

* For k ≥ 2 the first Clifford is absorbed by conjugation, so the words C_k T ⋯ C_1 T cover every case.
* **The rule for one gate (exact):** C·T violates iff both of these hold:
  * the symplectic part of C is a unit shear [[1,0],[c,1]], i.e. C fixes the Z axis, so C is a diagonal Clifford times
    a Pauli;
  * the Pauli part shifts X.

  That is 3 shears × 6 shifts = 18.
* The parity-twisted shears [[2,0],[c,2]] **never** violate.
* This refines Pass 11213's "cubic phase plus nonzero shift".

## Two qutrits: all 51 840 × 81 = 4 199 040 Cliffords C, U = C·(T⊗I)

A uniform sample of 200,000 two-qutrit Clifford cosets gives **18,248 violations**, probability **0.09124 ± 0.00064** (one standard error). This is sampled, not exhaustive; the producer estimates about 37 CPU-hours for the full 4,199,040-case census.

* The one-qutrit rule does **not** generalise. Most violators have a Clifford that moves Z₁, and the predicate
  "C Z₁ C† fails to commute with Z₁ ⇒ reversible" is false (268 counterexamples among 5000 samples).
* One sub-rule survives in every sample: **if C maps Z₁ to Z₁⁻¹, C·T₁ is reversible** (59/59). This is the two-qutrit
  echo of the parity-twisted shears.
* The full rule for n ≥ 2 is OPEN.

## Three qutrits (2000 random Cliffords)

A uniform sample of 2,000 three-qutrit Cliffords gives **190 violations**, probability **0.095 ± 0.0066** (one standard error). This is a reconnaissance sample only.

**Reading.**
* A single magic gate breaks substrate time reversal with probability **1/12** on one qutrit, rising slowly with the
  number of qutrits.
* Each further gate raises the violating fraction sharply: 1/12, 13/32, 421/768 on one qutrit.
* Compare Pass 11268: perfect scramblers are all-or-nothing in where the magic sits.
