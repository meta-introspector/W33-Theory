# Pass 11158 — perfect three-qutrit gates exist, and the determinant law generalises

Producer: `analysis/w33_pass11158_perfect_three_qutrit_gates.py`
Scan: `analysis/w33_pass11158_scan_sp63.py`
Frozen: `data/w33_pass11158_sp63_sample.json`
Regression: `tests/test_w33_pass11158_perfect_three_qutrit_gates.py`

**General law (proved).** For S ∈ Sp(2n, 3) with 2×2 blocks S_ij (output party i, input party j):

    Σᵢ det S_ij = 1 (mod 3) for every column j,   and   Σⱼ det S_ij = 1 for every row i.

This is the symplectic condition restricted to one party, since ω(Mu, Mu′) = det(M)·ω(u, u′) on F₃². For n = 2 it is Pass
11156's det S_AA + det S_BA = 1.

**Three qutrits** (Sp(6,3), 9 170 703 360 elements; 200 000 uniform samples built as products of 80 transvections):
* The column and row laws hold in 100% of samples.
* Perfect means the Choi state is maximal across all 3|3 cuts (AME(6,3)). This holds exactly when all nine blocks are
  invertible; the 2×2 block minors then follow automatically.
* **Perfect gates exist:** a fraction of 0.03108 ± 0.0004, about 2.85·10⁸ gates.
* Every perfect gate has determinants (1, 1, −1) in each column, the only way three nonzero values can sum to 1 mod 3.

The two-qutrit rule "both −1" becomes "exactly one −1 in each column". The existence of AME(6,3) is known in general; new
here are its realisation in W(3,3)'s three-qutrit Clifford group and the determinant law. The exact count is open.
