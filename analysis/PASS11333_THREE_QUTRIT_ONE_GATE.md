# Pass 11333 — one magic gate on three qutrits: the bad fraction is 0.111 ± 0.009 (consistent with 1/8)

Producer: `analysis/w33_pass11333_three_qutrit_one_gate.py`
Certificate: `data/w33_pass11333_three_qutrit_one_gate.json`
Regression: `tests/test_w33_pass11330_11334.py`

**Question.** The fraction of symplectic classes M for which some Pauli frame of C = W(a)·V_M makes C·(T⊗I⊗I) violate
substrate time reversal is:
* exactly 1/8 on one qutrit (3/24, Pass 11266);
* exactly 1/8 on two qutrits (6480/51 840, Pass 11330, exhaustive).

Does 1/8 persist on three qutrits?

**Method.**
* |Sp(6,3)| ≈ 9.2×10⁹, so orbit reduction is out of reach. 1100 classes M are sampled from long random words in the
  Sp(6,3) generators.
* The canonical Weil unitary V_M is built by twirl. Each class is decided exactly (Pass 11252) at the 7 frames
  0, e₁, …, e₆.
* At n = 2, Pass 11330 proved this five-frame screen exact by exhausting every frame. At n = 3 it is exact **if** the
  affine law holds there, which is checked on all 729 frames for two flagged classes.
* Uniformity control: the fraction of sampled M fixing z₁ is 1/1100, against 1/728 expected.

**Results.**

| | |
|---|---|
| flagged bad classes | **122 / 1100 = 0.111 ± 0.009** |
| z-score against 1/8 | −1.4 |
| full 729-frame checks of 2 flagged classes | 486 = (2/3)·729 violating frames each; the non-violating frames form an affine subspace (codimension 1) |

**Reading.**
* **Consistent with 1/8.** The sample also cannot exclude 1/9.
* **The affine law survives at n = 3** on the two classes checked exhaustively. Its general n = 3 status is open.
* **What a sharper test needs.** Mid-range sampling cannot separate the two candidates. A test of the
  W(3,3)-geometry mechanism of Pass 11331 can: at n = 3 it predicts that every class fixing z₁ is bad and every class
  reversing z₁ is good. That needs targeted sampling of those cells, not uniform sampling.
