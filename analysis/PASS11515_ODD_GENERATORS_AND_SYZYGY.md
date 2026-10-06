# Pass 11515 — the time arrows of a qutrit state are the (k, 2, 1) MUB chiralities, k = 3 … 7, with one relation in degree 13

Producer: `analysis/w33_pass11515_odd_generators_and_syzygy.py`
Certificate: `data/w33_pass11515_odd_generators_and_syzygy.json`
Regression: `tests/test_w33_pass11511_11515.py`

## Explicit generators

Pass 11502 counted one new time-odd generator in each degree 6 … 10 over the even ring. Here they are written down.
* **Search.** For each degree d, the stabiliser-probability monomials were searched in order: fewest distinct states
  first, then exponent pattern. The first monomial m whose odd Reynolds average
  R_m = (1/432) Σ_g sgn(g) m(gψ) lies outside the span of {even invariants × lower generators} was taken.
* **Result.** Every one of the five has the same shape:

> **g_{k+3} = Σ_g sgn(g) · p_{ga}^k p_{gb}² p_{gc},  k = 3, 4, 5, 6, 7,**
>
> with a, b, c stabiliser states in three different mutually unbiased bases (Z, X, XZ).

| degree | monomial | exponents | new-direction gap |
|---|---|---|---|
| 6 | p_{Z0} p_{X0}² p_{Y0}³ | (3, 2, 1) | 1 (first) |
| 7 | p_{Z0} p_{X0}² p_{Y0}⁴ | (4, 2, 1) | 1.3·10⁻² |
| 8 | p_{Z0} p_{X0}² p_{Y0}⁵ | (5, 2, 1) | 3.0·10⁻³ |
| 9 | p_{Z0} p_{X0}² p_{Y0}⁶ | (6, 2, 1) | 2.2·10⁻⁴ |
| 10 | p_{Z0} p_{X0}² p_{Y0}⁷ | (7, 2, 1) | 1.2·10⁻⁵ |

Here Y = XZ. The gap column is the relative singular value of the new direction; each dropped value is 0 to machine
precision.

**Reading.**
* The time-odd invariants of a qutrit state are generated, through degree 13, by **one family**: the chiralities of an
  oriented triple of MUBs, weighted p^k·p²·p.
* h₆ (Pass 11434) is k = 3, and g₇ (Pass 11491) is k = 4.

## The first relation

**Prediction from the Molien numbers alone.**
* If the odd module were free on generators of degrees 6 … 10 over the even ring E, degree k would hold
  Σ_i E_{k−d_i} products.
* With Pass 11491's exact E_k and O_k this count is 1, 2, 3, 5, 8, 11, 16 for k = 6 … 12, equal to O_k.
* At k = 13 it is **23 against O₁₃ = 22**.

**Verified** with the explicit generators and an even basis:

| k | 6 | 7 | 8 | 9 | 10 | 11 | 12 | **13** |
|---|---|---|---|---|---|---|---|---|
| products e·g_i | 1 | 2 | 3 | 5 | 8 | 11 | 16 | **23** |
| rank | 1 | 2 | 3 | 5 | 8 | 11 | 16 | **22** |
| relations | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** |

* The rank gap at degree 13 runs from 2.3·10⁻⁷ (kept) to 3·10⁻¹⁶ (dropped).

> **The module of time arrows is free on the five chiralities through degree 12; its first syzygy is a single relation in
> degree 13.**

**Scope.**
* Ranks are numerical, with gaps of six or more decades. The dimensions are exact (Pass 11491).
* "Generated through degree 13" does not exclude generators above 13 (Pass 11502).
* The relation is located but not written out.

## Correction (applies to Passes 11491 and 11502)

* **The error.** Pass 11491 said the new degree-8 odd invariant "needs two stabiliser states from the same MUB", and Pass
  11502 repeated it.
* **Why it was wrong.** 11491's search returned the candidate with the largest residual, and that candidate happened to
  use a same-MUB pair. "Needs" was an over-read.
* **The fact.** The (5, 2, 1) chirality over three *different* MUBs is already a new generator in degree 8.
