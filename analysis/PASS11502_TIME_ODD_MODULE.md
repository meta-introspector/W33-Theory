# Pass 11502 — Δ⁽²⁾ = Δ₆/108 fully proved, and the generator degrees of the time-even ring and time-odd module of a qutrit state

Producer: `analysis/w33_pass11502_time_odd_module.py`
Certificate: `data/w33_pass11502_time_odd_module.json`
Regression: `tests/test_w33_pass11498_11505.py`

## (A) The last step of the 1/108 derivation is now exact

* **What Pass 11492 had.** It derived Δ⁽²⁾(ψ⊗s) = (1/120 + (27/40)·3⁻⁶)·Δ₆(ψ) = Δ₆(ψ)/108 from the compressions
  A_C = (1 ⊗ ⟨0|) C (1 ⊗ |0⟩) of all 4 199 040 two-qutrit Cliffords. The vanishing of the rank-one terms' odd part was only
  numerical (3·10⁻¹⁶).
* **The census.** Every rank-one compression is |x⟩⟨y|, with unit singular value, and **x, y are stabiliser states**.
  The 1 259 712 of them are spread **exactly uniformly over the 144 ordered stabiliser pairs**, 8748 each.
* **Why that suffices.** Their total odd contribution is
  8748·(Σ_x p_x⁶)·(Σ_y p_y⁶ − Σ_y p_{ȳ}⁶) = 0,
  because complex conjugation permutes the 12 stabiliser states (s ↦ s̄).
* **So the identity is proved:**

> **Δ⁽²⁾(ψ ⊗ s) = Δ₆(ψ)/108 exactly, with 1/108 = 1/120 + 1/1080.**

* **What each piece is:**
  * 1/120 of the two-qutrit Clifford group compresses to a one-qutrit Clifford;
  * 27/40 compresses to 3^{−1/2} times one, uniformly over the one-qutrit Clifford group;
  * the rest compresses to rank one, uniformly over stabiliser pairs, and contributes nothing odd.

## (B) Module structure (numerical ranks, exact dimensions)

**Method.**
* The 12 stabiliser probabilities span all Hermitian quadratic forms. So the Reynolds averages
  R_m = Σ_g sgn^ε(g) m(gψ) of p-monomials span every even (ε = 0) and odd (ε = 1) invariant.
* Bases of E_k and O_k were built at 160 random unit rays. A bidegree-(k, k) form is fixed by its values on the sphere.
  The basis sizes are checked against Pass 11491's proved Molien series.
* **Generators in degree k** are counted as
  * E_k − rank{e_i e_j}, the products of even invariants of positive degrees summing to k;
  * O_k − rank{e o}, with e even of degree ≥ 1 and o odd.
* **Gaps.** At every degree the kept singular values are ≥ 10⁻⁷ and the dropped ones ≤ 10⁻¹³ (relative).

| k | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E_k | 1 | 1 | 2 | 3 | 4 | 6 | 8 | 11 | 15 | 19 | 24 | 32 | 40 |
| O_k | 0 | 0 | 0 | 0 | 0 | 1 | 2 | 3 | 5 | 8 | 11 | 16 | 22 |
| new even generators | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 0 | 0 |
| new odd generators | — | — | — | — | — | **1** | **1** | **1** | **1** | **1** | 0 | 0 | 0 |

**Reading.**
* **The even ring** (the time-symmetric observables of a qutrit state) needs, through degree 13, generators in degrees
  **1, 3, 4, 5, 6, 7, 8, 9**: eight in all, one per degree, with none in degree 2.
  * The degree-1 generator is the norm |ψ|².
  * D₂ = 1 means the norm squared is the only degree-2 invariant, as for any 2-design.
* **The odd module** (the time arrows) needs exactly **one new generator in each degree 6, 7, 8, 9, 10**, and none in 11–13.
  * h₆ is the (3, 2, 1) MUB chirality (Pass 11434).
  * g₇ is the (4, 2, 1) chirality (Pass 11491).
  * The degree-8 generator needs two stabiliser states from one MUB (Pass 11491).
  * The degree-9 and degree-10 generators are not written down.
* So the time arrows of a qutrit state form a module on **five generators of consecutive degrees 6–10**, as far as computed.
* This corrects the free-module reading suggested by the rational form t⁶(1+t)/∏(1−t^d). That reading would have only
  two generators (6, 7), and Pass 11491 had already shown its primary degrees cannot exist.

**Scope.**
* The ranks are numerical, but the gaps are about six decades and the dimensions are exact.
* Generator degrees found up to degree K do not exclude generators above K. No a-priori degree bound is used: the
  acting group includes the U(1) phase, so the finite-group Noether bound does not apply as stated.
* Only the counts are claimed; the generators are not written down.

## Correction (Pass 11515)

* **The error.** The statement that the new degree-8 odd invariant **needs** two stabiliser states from the same MUB is
  wrong: it was an over-read of a largest-residual search.
* **The fact.** The (5, 2, 1) chirality R(p_a⁵ p_b² p_c) over three *different* MUBs is already a new generator in
  degree 8.
* **The full picture.** All five odd generators are the (k, 2, 1) chiralities for k = 3 … 7. The module is free
  through degree 12, with one relation in degree 13 (Pass 11515).
