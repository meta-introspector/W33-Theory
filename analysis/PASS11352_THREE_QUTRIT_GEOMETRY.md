# Pass 11352 — one magic gate on three qutrits: the bad fraction is 0.1333 ± 0.0019 (not 1/8), and the affine law fails

Producers: `analysis/w33_pass11352_three_qutrit_geometry.py`, `analysis/w33_pass11352_affine_failure.py`
Certificates: `data/w33_pass11352_three_qutrit_geometry.json`, `data/w33_pass11352_affine_failure.json`
Regression: `tests/test_w33_pass11350_11354.py`

**Method.**
* Pass 11350's exact linear decider is applied to uniformly sampled Sp(6,3) classes (long random words in the H, S, SUM
  generators), in two independent batches of 8000 and 24 000 classes.
* The rare geometric cells are also targeted: M is composed with transvections carrying Mz₁ to ±z₁.
* **Independent check:** 60 classes (10 of them bad) × 11 frames are decided with the Weyl criterion of Pass 11252.
  **660/660 agree.**

**The bad fraction is not 1/8.**

| | |
|---|---|
| decided classes | 31 994 (+ 6 beyond the enumeration cap, reported, not guessed) |
| bad | 4265 |
| **fraction** | **0.1333 ± 0.0019** |
| against 1/8 | **+4.4σ (excluded)** |
| against 1/9 | +11.7σ (excluded) |

* So the 1/8 shared by one and two qutrits does not persist to three.
* The value is consistent with 2/15, but other simple fractions (e.g. 21/160) lie within about 1σ. No closed form is
  claimed.
* Pass 11333's slower estimate 0.111 ± 0.009 was low by 2.2σ.

**Which W(3,3) rules carry over** (pooled uniform and targeted counts).

| cell | n = 2 (Pass 11331) | n = 3 |
|---|---|---|
| Mz₁ = z₁ | all bad | **187/187 bad** |
| Mz₁ = −z₁ | all good | **193/193 good** (2 undecided) |
| same line, M²z₁ = z₁ | all bad | **42/42 bad** |
| same line, M²z₁ = −z₁ | all good | **36/36 good** |
| collinear, different lines | all good | **806/6990 bad**: the rule does not carry over |
| same line, other | 1/2 | 709/3373 |
| non-collinear | 1/9 | 2671/21 472 |

**The affine law fails at n = 3.**
* Two classes (seeds 619187 and 622316) have 567 reversible frames. 567 is not a power of 3, so these frames do not
  form an affine subspace.
* Their reversible frames are **exactly a union of affine hyperplanes** (5 hyperplanes contained).
* Each was confirmed with the Weyl criterion on 12 frames (6 inside, 6 outside): **12/12 agree**.
* This is the "union of affine pieces" structure predicted by Pass 11350. At n = 2 the union always collapsed to a
  single piece; at n = 3 it does not.
* Pass 11333's single-subspace finding (two classes, 486 violating frames) was a special case.
