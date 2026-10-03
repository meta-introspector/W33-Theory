# Pass 11358 — the three-qutrit one-gate mechanism: two universal rules and two that drift

Producer: `analysis/w33_pass11358_three_qutrit_mechanism.py`
Certificate: `data/w33_pass11358_three_qutrit_mechanism.json`
Regression: `tests/test_w33_pass11355_11360.py`

**Data.** 3000 uniformly sampled Sp(6,3) classes, 2998 decided by the F₃-linear decider of Pass 11350. Each class was
tabulated by the cell of the magic axis z₁ and by the number of symplectic solutions of
(S): M s^{−k} Q (JMJ) = Q, Q z₁ = z₁. Bad fraction 0.1341, consistent with Pass 11352's 0.1333 ± 0.0019.

**A remark (exact).** (S) is unchanged under M ↦ −M, since both factors flip sign. So the parity −I acts only on the
frame condition (F).

**Rules that hold at n = 2 and n = 3.**
* **(S) has no symplectic solution ⇒ violating in every frame.** This is a theorem (Pass 11350): no Clifford E exists.
  * At n = 2 such classes occur only in the single-line cell (1296).
  * At n = 3 they occur also in the **non-collinear** cell (43 of 2998 sampled). This is a source of badness that
    n = 2 lacks.
* **(S) has exactly one solution ⇒ reversible in every frame.** Verified on the whole n = 2 non-collinear cell (23 328 classes, Pass 11351) and on 2021/2021 sampled
  n = 3 classes (1379 non-collinear, 441 collinear on different lines, 201 single-line). Unproved.

**Rules that drift.**

| | n = 2 | n = 3 |
|---|---|---|
| three solutions | exactly half bad (3888/7776) | non-collinear 147/351 = 0.42; collinear, different lines 78/145 = 0.54 |
| −I swap (bad M ⇔ good −M) on the non-collinear cell | perfect | **66 both-bad pairs** among 2029 |
| −I swap on the collinear-different-lines cell | n/a (all good) | perfectly complementary (78 + 67, no both-bad) |

**Reading.** The excess of 0.1333 over 1/8 comes from:
* symplectically unsolvable classes appearing off the single-line cell; and
* the three-solution half-rule and the −I exchange loosening.

The n = 2 mechanism is the skeleton, and n = 3 adds new obstruction types. A closed count of the n = 3 fraction is open.
