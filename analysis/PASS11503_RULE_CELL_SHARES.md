# Pass 11503 — the exact shares of the four rule cells, for every n

Producer: `analysis/w33_pass11503_rule_cell_shares.py`
Certificate: `data/w33_pass11503_rule_cell_shares.json`
Regression: `tests/test_w33_pass11498_11505.py`

## Statement

With N = 3^{2n} − 1:

| rule cell | share of Sp(2n, 3) | n = 2 classes | rule |
|---|---|---|---|
| Mz₁ = z₁ (inside the next) | 1/N | 648 | bad |
| M²z₁ = z₁, v ≠ 0 | **2/N** | 1296 | bad: the magic-axis law (Pass 11498) |
| Mz₁ = −z₁ | 1/N | 648 | good (Pass 11499) |
| M²z₁ = −z₁ | **3/N** | 1944 | good (Pass 11499) |

**Derivation.**
* Mz₁ is uniform on the N nonzero vectors.
* By Witt's theorem, Sp(2n, 3) is transitive on independent ordered pairs (u, w) with a given ω(u, w).
* For y = Mz₁ independent of z₁, M²z₁ = ±z₁ requires My = ±z₁. Then ω(z₁, y) = ω(Mz₁, My) = ±ω(y, z₁).
  * **Plus sign:** this forces ω(z₁, y) = 0. The isotropic y contribute 1/N in total.
  * **Minus sign:** ω is left free. The two non-isotropic values contribute 2/N, and the isotropic one contributes 1/N.
* Adding y = z₁ (to the + cell) gives 2/N and 3/N.

**Checks.**
* n = 2: counted over all 51 840 elements.
* n = 3: Pass 11373's exact orbit sizes over all 9 170 703 360 elements.
* Every share matches.

## What is new and what is not

* **Already in the corpus (cite).** Pass 11373's table gives each of the singleton cells (fixed, reversed, same-line
  M²z₁ = ±z₁) one vector out of N. The bad share 2/N = fixed + same-line follows from it.
* **New.** The good rule cell **M²z₁ = −z₁ straddles Pass 11373's cells.**
  * 1/N lies in "same line".
  * **2/N lies inside the non-collinear cell.**
  * So the non-collinear cell, whose bad fraction is 1/9 at n = 2 and 10/81 at n = 3, contains a sub-cell of share 2/N
    that is **always good**.
  * Pass 11499 proves that sub-cell good at n = 2 on every class.
* **Scope.** The four rule cells together have share 6/N, which shrinks with n. They govern none of the bulk of the
  ≈ 1/8 bad-class fraction.
