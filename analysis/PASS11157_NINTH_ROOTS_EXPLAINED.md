# Pass 11157 — why ninth roots of unity: one identity explains both cubics

Producer: `analysis/w33_pass11157_ninth_roots_explained.py`
Regression: `tests/test_w33_pass11157_ninth_roots_explained.py`

Pass 11152 found char(B) = (x³ − 6x − 2)(x³ − 3x + 1)² for the best Pauli-only Bell operator (the Z and X bases on both
qutrits), but left the ninth roots of unity unexplained.

* **Weyl form.** B is a sum of eight two-qutrit Weyl operators: Z⊗Z², Z⊗X², X⊗Z², X⊗X² and their adjoints. Each has
  coefficient of modulus 1/√3 and phase ±π/6 or ±5π/6.
* **Symmetry.** The only nontrivial Pauli symmetry is **U = (XZ²)⊗(XZ²)**, a joint relabelling of both parties' outcomes
  with U³ = I. B splits into three 3-dimensional charge sectors.
* **The identity** (40-digit check; exact, since every entry lies in (1/9)ℤ[ω]):

      B³ − 3B + 1 = (B + 1)(I + U + U²).

  * **Charged sectors (U = ω, ω²):** the right side vanishes, so B³ − 3B + 1 = 0. With B = 2cos θ this reads
    2cos 3θ = −1 = 2cos(2π/3): the Chebyshev **trisection** of the symmetry's angle 2π/3. So θ = 2πm/9 for m = 1, 2, 4,
    which is where the ninth roots come from.
  * **Neutral sector (U = 1):** B³ − 3B + 1 = 3(B + 1), i.e. B³ − 6B − 2 = 0. Its largest root 2.6017 is the Pauli-only
    maximum.

Both cubics are one identity read in different charge sectors of the relabelling symmetry. This closes the question left
open in 11152.
