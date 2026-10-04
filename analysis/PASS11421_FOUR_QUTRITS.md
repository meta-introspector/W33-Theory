# Pass 11421 — four qutrits, one magic gate: the four magic-axis rules hold at n = 4, and the bad-class fraction is not monotone in n

Producer: `analysis/w33_pass11421_four_qutrits.py`
Certificate: `data/w33_pass11421_four_qutrits.json`
Regression: `tests/test_w33_pass11418_11422.py`

**Method.**
* The F₃-linear decider of Pass 11350 is exact for every n. Its operator table (6561 Weyl operators of size 81)
  would cost 690 MB at n = 4.
* `SparseDecider` replaces the table by the closed form W(p)|j⟩ = ω^{2x·z + z·j}|j + x⟩.
  * The Weil unitary is assembled from 6561 phased rank-one entries.
  * Frames are read off one column and n phase ratios.
* **Validated** against the table-based decider: identical good-frame sets on all 60 classes checked at n = 2 and all
  25 at n = 3.
* Sp(8,3) has about 6×10¹⁷ elements, so this is a sample, not a census. Classes whose (S) solution space exceeds the
  cap are reported as undecided, never guessed.

## The four rules at n = 4 (targeted classes)

| cell | construction | bad | good | undecided |
|---|---|---|---|---|
| Mz₁ = z₁ | random symplectic composed with a transvection map onto z₁ | **59** | 0 | 1 |
| Mz₁ = −z₁ | the same, onto −z₁ | 0 | **56** | 4 |
| same line, M²z₁ = z₁ | an exact n = 2 class A₂ ⊕ B, B ∈ Sp(4,3), conjugated by random N ∈ Stab(z₁) | **46** | 0 | 14 |
| same line, M²z₁ = −z₁ | the same | 0 | **45** | 15 |

**No exception on any decided class.** The four rules of Pass 11373 hold at n = 1, 2 and 3 exactly, and at n = 4 on
every decided targeted class.

## Uniform sample (3000 classes, all decided)

* **Bad-class fraction 0.1183 ± 0.0059.**
* By cell:

| cell | bad / sampled |
|---|---|
| collinear, different lines | 78/663 = 0.118 |
| non-collinear | 233/1998 = 0.117 |
| same line, other | 44/339 = 0.130 |

* Bad frames per class: 4374 (= ⅔ of 6561) in 270 classes and all 6561 in 83 classes. 5832 and 1458 occur once each.

**Reading.**
* The bad-class fraction runs 1/8, 1/8, 437/3276 = 0.1334, and 0.118 ± 0.006 for n = 1, 2, 3, 4.
* n = 4 is consistent with 1/8 and lies 2.6σ below the exact n = 3 value. The n = 3 excess over 1/8 does **not**
  grow monotonically. A closed form in n is not in sight.
* The cell structure keeps the n = 3 pattern: collinear and non-collinear cells both carry about 1/9 to 1/8.
