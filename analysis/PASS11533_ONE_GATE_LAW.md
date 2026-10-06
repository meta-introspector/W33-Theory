# Pass 11533 — the one-gate law: violating frames by linear algebra, a rigorous lower bound on Pₙ, and the 1/8 conjecture refuted

Producer: `analysis/w33_pass11533_one_gate_law.py` (stages `validate`, `est n samples`)
Certificate: `data/w33_pass11533_one_gate_law.json`
Regression: `tests/test_w33_pass11531_11538.py`

## The law

**Definition.** For the tick U = W(a)V_M T₁, call v **universally anti-fixed** if three conditions hold:
* Mv = v;
* ω(v, z₁) = 0;
* QJv = −v for every Q in every solution space of (S) (k = 0, 1, 2).

These conditions are linear in v, so the universally anti-fixed vectors form a **subspace W**.

**Theorem (sufficiency, every n).** Every frame with ω(v, a) ≠ 0 for some v ∈ W violates. If (S) has no solution at all,
every frame violates.
* The proof is Pass 11498's argument with z₁ + Mz₁ replaced by any v ∈ W.
* ω(v, z₁) = 0 makes the shear fix v and removes the T-twist; r_Q = 0; pairing (F) with ω(v, ·) leaves 2ω(v, a) = 0.

**Predictions.**
* The law predicts 3^{2n} − 3^{2n − dim W} violating frames.
* It predicts the class is bad iff W ≠ 0 or (S) has no solution.
* Computing W is a few linear solves, about 7 ms per class at n = 4. The exact decider needs seconds to minutes per class.

## Validation against the exact deciders

**n = 2 (all 51 840 classes; 2 beyond the decider's cap).**
* The class verdict agrees on **51 838/51 838**.
* The violating-frame count agrees on 51 814.
* The 24 exceptions are the fixed-axis transvections (Pass 11486). There the true violating set is a union over a
  2-dimensional space containing W.
* **Law P₂ = 6480/51840 = 1/8, exactly.**

**n = 3 (all 2308 orbits of Pass 11373, with exact weights).**
* The verdict agrees on 2190 of the 2199 decided orbits.
* **Law P₃ = 97/728**, against the exact 437/3276. The difference is exactly **1/6552**.
* The deficit is entirely **9 orbits of the "same line, other" cell**. Each has 162 violating frames forming an affine
  piece: {ω(w, a) = 0, ω(v, a) ≠ 0} for two M-fixed vectors w, v ⊥ z₁.
* Pass 11537's union law accounts for them.
* In every other cell the law's fraction is exact: non-collinear 10/81, collinear-different-lines 1/9, and the four rule
  cells.

## Precise bad-class fractions (the law gives a rigorous lower bound)

Every class the law calls bad is bad, so on any sample the law's fraction is at most the true one.

Samples: 1.2·10⁶ for n = 4 and 5, 4·10⁵ for n = 6. Each entry is the law's fraction ± standard error.

| n | Pₙ (law, lower bound) | non-collinear | collinear, different lines | same line, other |
|---|---|---|---|---|
| 2 | **1/8** (exact, all classes) | 1/9 | 0 | 1/2 |
| 3 | **97/728 = 0.13324** (exact; true 437/3276) | 10/81 | 1/9 | 17/76 |
| 4 | 0.13371 ± 0.00031 | 0.13267 ± 0.00038 | 0.12264 ± 0.00064 | 0.16036 ± 0.00101 |
| 5 | 0.13417 ± 0.00031 | 0.13358 ± 0.00038 | 0.13291 ± 0.00066 | 0.13998 ± 0.00095 |
| 6 | 0.13468 ± 0.00054 | 0.13495 ± 0.00066 | 0.13280 ± 0.00114 | 0.13680 ± 0.00163 |

* **The bad-class fraction does not settle at 1/8.** The law's lower bound alone puts P₄, P₅ and P₆ 28σ, 29σ and 18σ
  above 0.125.
* **It rises slowly:** 0.1332, 0.1337, 0.1342, 0.1347.
* **The cells converge.** The three bulk cells, far apart at n = 2 and 3, approach each other near 0.133 to 0.137.
* **Earlier decider samples, now superseded.** They were P₄ = 0.128 ± 0.005, P₅ = 0.138 ± 0.004 and P₆ = 0.131 ± 0.006.
  All are consistent with these bounds within 1.1σ; the law gives errors 10 to 15 times smaller.
* **Scope.** These are lower bounds. The W-law misses exceptional classes, of mass 1/6552 at n = 3. The union law
  (Pass 11537) is exact where it has been checked, but it is too slow to sample at this scale.

## The refuted conjecture

* The exact non-collinear fractions 1/9 (n = 2) and 10/81 (n = 3) fit (1 − 9^{1−n})/8 = Σ_{j=1}^{n−1} 9^{−j}, whose limit
  is 1/8.
* At n = 4 the law's non-collinear fraction is **0.13267 ± 0.00038**, against the predicted 0.12483: **z = +20.7**.
* The law is exact in that cell at n = 2 and 3. Its value is also a lower bound on the true fraction, so the true
  fraction is at least this.
* **Refuted.** The non-collinear fraction does not converge to 1/8; it *increases* from n = 3 to n = 4.
* The structure also differs from what the formula suggested. At n = 3 the non-collinear bad mass splits as
  8/81 from W ≠ 0 plus 2/81 from classes with no solution of (S). That is not the 1/9 + 1/81 the geometric sum implied.
