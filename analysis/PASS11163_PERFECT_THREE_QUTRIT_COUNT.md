# Pass 11163 — exactly 286 654 464 = 2¹⁷·3⁷ perfect three-qutrit Clifford gates (128/4095 of Sp(6,3))

Producer: `analysis/w33_pass11163_perfect_three_qutrit_count.py`
Scan: `analysis/w33_pass11163_scan_perfect_count.py`, frozen as `data/w33_pass11163_perfect_count.json`
Regression: `tests/test_w33_pass11163_perfect_three_qutrit_count.py`

**Criterion (now proved).** A three-qutrit Clifford gate S is perfect, meaning its Choi state is AME(6,3), iff all nine
2×2 party blocks are invertible.
* A 3|3 cut with k inputs needs the 2k×2k submatrix on k output and k input parties to be invertible.
* For k = 2, Jacobi's complementary-minor identity, with det S = 1 and S⁻¹ = −JSᵀJ, gives
  det S[{i,i′},{j,j′}] = ± det S_{i″j″}, the complementary block.
* Pass 11158 had stated this step as empirical; it is now proved.

**Count.** S is an ordered symplectic basis (s1,s2 | s3,s4 | s5,s6), and the determinant of block i of a column pair
(a,b) is the local form w_i(a,b).
* There are 176 904 symplectic pairs (a,b); 41 472 of them have all three local forms nonzero.
* The admissible second pairs (c,d) lie in the complement and satisfy two conditions: w_i(c,d) ≠ 0, and the third
  column's determinants 1 − w_i(a,b) − w_i(c,d) (fixed by the row law) are nonzero.
* That gives 11 943 936 good (first, second) pairs.
* The third pair is then any symplectic basis of the remaining plane: 24 choices, block determinants unchanged.

    **#perfect = 24 × 11 943 936 = 286 654 464 = 2¹⁷·3⁷ = 24³·144²,   fraction 128/4095 = 2⁷/(2¹²−1).**

The two-qutrit fraction is 4/15 = 2²/(2⁴−1).

**Checks.**
* The Pass 11158 sampling estimate was 0.03108 ± 0.0004; the exact value is 0.0312576.
* **Known-positive control:** the same method for two qutrits returns 13 824, the Pass 11156 enumeration.

**Prior art.** Pahari (arXiv:2607.00210) classifies *two-qudit* Clifford dual-unitary gates over F_q. It finds q − 2
perfect-tensor cores, which is one at q = 3 and matches Pass 11156's single local orbit, but gives no three-qudit count.
We found no exact count of perfect three-qutrit Clifford gates in the literature we checked.
