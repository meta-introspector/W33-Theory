# Pass 11373 — one magic gate on three qutrits, exactly: 240857/2388204 of all ticks and 437/3276 of all classes break the arrow

Producer: `analysis/w33_pass11373_three_qutrit_exact_fraction.py`
Data: `data/w33_pass11373_orbits_n{2,3}.txt` (GAP orbit lists), `data/w33_pass11373_centralisers_n{2,3}.txt`
Certificate: `data/w33_pass11373_three_qutrit_exact_fraction.json`
Regression: `tests/test_w33_pass11369_11373.py`

**Question.**
* Pass 11352 sampled the three-qutrit one-gate violating class fraction as 0.1333 ± 0.0019 and excluded the
  two-qutrit 1/8.
* Pass 11358 found which mechanisms drift.
* Is the fraction exact, and which rules survive?

## Method (exact)

* **The ticks.** U = W(a) V_M (T ⊗ I ⊗ I), with M ∈ Sp(6,3) (9 170 703 360 classes) and 729 frames a.
* **Orbit reduction.** Let H = Stab(z₁), of order 12 597 120.
  * Every N ∈ H lifts to a Clifford commuting with T₁, so the number of violating frames is constant on
    H-conjugation orbits of Sp(6,3).
  * **GAP** enumerates the orbits as C_G(g)\G/H double cosets over the 141 conjugacy classes: **2308 orbits**. Their
    sizes sum to |Sp(6,3)|.
* **Fast classes (2199).** The F₃-linear decider of Pass 11350.
* **Slow classes (109).** These have high symmetry, and their (S) solution spaces have affine dimension 14–30.
  * For N in the centraliser C_H(M), V_N(W(a)V_M T₁)V_N† = W(Na) V_M T₁ up to phase. So the violating frames are
    unions of orbits of the linear action a ↦ Na.
  * GAP gives generators of C_H(M), checked against |H : C_H(M)| = orbit size.
  * Each class splits into 6 to 63 frame orbits. One frame per orbit is decided by the exact Weyl criterion
    (Pass 11252).
* **Controls.**
  * The whole pipeline at n = 2 (200 orbits) reproduces **223/2430** and **1/8** exactly (Passes 11309/11330).
  * For every slow class, two random non-representative frames were decided directly. All **218** agree with their
    orbit's verdict.

## Results

> **Violating fraction over all ticks: 240857/2388204 ≈ 0.100853.**
> **Fraction of classes with at least one violating frame: 437/3276 ≈ 0.133394.**

The class fraction agrees with Pass 11352's sample (0.1333 ± 0.0019) at z = 0.05. Two qutrits give 1/8 and 223/2430.

**By the W(5,3)-type geometry of the magic axis** (bad-class fraction):

| cell (where M moves z₁) | share of Sp(6,3) | n = 2 | **n = 3** | violating fraction, n = 3 |
|---|---|---|---|---|
| fixes it: Mz₁ = z₁ | 1/728 | 1 | **1** | 4408/6561 |
| reverses it: Mz₁ = −z₁ | 1/728 | 0 | **0** | 0 |
| same line, M²z₁ = z₁ | — | 1 | **1** | 164/243 |
| same line, M²z₁ = −z₁ | — | 0 | **0** | 0 |
| same line, other | — | 1/2 | **77/342** | 649/3078 |
| collinear, different lines | — | **0** | **1/9** | 2/27 |
| non-collinear | 486/728 | **1/9** | **10/81** | 1784/19683 |

**Geometric decomposition (exact).**
* Every cell is a set of images Mz₁ among the 3^{2n} − 1 nonzero vectors, so each cell's share is a vector count
  over 80 (n = 2) or 728 (n = 3).

| cell | vectors, n = 2 | vectors, n = 3 |
|---|---|---|
| fixed / reversed / M²z₁ = ±z₁ | 1 each | 1 each |
| same line, other | 4 | 76 |
| collinear, different lines | 18 | 162 |
| non-collinear | 54 | 486 |

* Hence:

> 1/8 = [2 + 4·½ + 18·0 + 54·⅑] / 80,  437/3276 = [2 + 76·(77/342) + 162·⅑ + 486·(10/81)] / 728 = (80 + 154/9)/728.

* The "2" is the two always-bad singleton cells; it is universal. At n = 3 the collinear and non-collinear cells
  contribute 18 + 60 = 78 where n = 2 had 0 + 6.

**Reading.**
* **Four rules are universal, now exactly at n = 3:**
  * fixing the magic axis always breaks the arrow;
  * reversing it never does;
  * M²z₁ = z₁ on the same line always does;
  * M²z₁ = −z₁ never does.
* **The 1/9 migrates.**
  * At n = 2 it sat on the non-collinear cell, and the collinear-different-lines cell was clean.
  * At n = 3 the collinear cell carries exactly 1/9, and the non-collinear cell carries 10/81 = 1/9 + 1/81.
  * The same-line half-rule (1/2) becomes 77/342.
* The bad-frame counts per class are {0, 162, 486, 540, 648, 729}. The good-frame counts 567 and 189 are not powers
  of 3, so the affine law fails exactly where Pass 11352 found it failing.
* Pass 11358's sampled rules ("no (S) solution ⇒ violating", "one solution ⇒ reversible") are consistent with this
  census. Their proof is still open.

## A proof route for "reversing the magic axis never breaks the arrow" (`--route`)

* At k = 0, equation (S) of Pass 11350 reads Q⁻¹MQ = J M⁻¹ J, with Q z₁ = z₁.
* Any anti-symplectic involution σ with σMσ = M⁻¹ and **σz₁ = −z₁** gives the solution Q = σJ.
  * Wonenburger: every symplectic matrix is reversed by an anti-symplectic involution. The constraint σz₁ = −z₁ is
    the new part.
* **n = 2, the whole reversed cell (computed):**
  * all **648/648** classes have such a σ;
  * for **621** of them a single σ already makes all 81 frames reversible;
  * for the other 27 the best single σ covers 27 frames, and the union over solutions covers the rest (census).
* So the rule splits into an **existence lemma**, which held in every class checked, and a **covering lemma** for the
  frame condition (F). Neither is proved here.

**Open.**
* A closed formula for 437/3276 = 437/(2²·3²·7·13) from the geometry of W(5,3).
* The general-n law: the fixed-point rules look universal, while the 1/9 cell moves.
