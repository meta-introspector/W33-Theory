# Pass 11536 — the good rules for the regular unipotent blocks of dimension 8

Producer: `analysis/w33_pass11536_good_rules_dim8.py` (the method of Pass 11532 in Sp(8, 3); flag `--no-x2`)
Certificate: `data/w33_pass11536_good_rules_dim8.json`
Regression: `tests/test_w33_pass11531_11538.py`

* **What is new in dimension 8.** The only new indecomposable unipotent class is the regular one, J₈. It comes in two
  classes, told apart by the square class of ω(v, N⁷v).
  * J₄ ⊕ J₄ and all smaller combinations split into blocks already settled by Passes 11511 and 11532.
  * There is no odd pair J_k ⊕ J_k with 2k = 8, since k must be odd.
* **The witnesses.** Random elements of Sp(8, 3) were drawn and their unipotent parts classified. Each J₈ representative
  found was given witnesses, exactly as in Pass 11532.

| class | (U) | (E) for −u |
|---|---|---|
| J₈, form class 1 | 4/4 proved | 4/4 proved |
| J₈, form class 2 | 4/4 proved | 4/4 proved |

Both classes were found within 287 samples.

* **Not done: the x² + 1 type in dimension 8.** This is the regular unipotent of U(4, 9). It did not occur in 800 random
  samples, because random elements of the centraliser rarely have a regular unipotent part. A direct construction is
  needed, so this case is **open**.
* **Status of the good rules.** They now hold for every n whenever:
  * the unipotent part at eigenvalue 1 has indecomposable pieces of dimension ≤ 6, or is J₈;
  * the z₁-block has pieces of dimension ≤ 6.
