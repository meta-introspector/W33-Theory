# Pass 11309 — one magic gate on two qutrits: an exact affine law and the exact probability 223/2430

Producer: `analysis/w33_pass11309_one_gate_law.py` (decisions by Pass 11252's exact criterion)
Certificate: `data/w33_pass11309_one_gate_law.json`
Regression: `tests/test_w33_pass11308_11312.py`

**Starting point.**
* Pass 11266, one qutrit: C·T violates substrate time reversal iff C's symplectic part is a unit shear fixing Z and its
  Pauli frame shifts X. That gives (3/24)·(2/3) = 1/12.
* For two qutrits it sampled 0.0912, and the naive fixed-axis rule failed.

## The law (two qutrits, exhaustive over all 4 199 040 Cliffords mod phase)

Write C = P_a·V_M, with a ∈ F₃⁴ the Pauli frame and M one of the 51 840 symplectic classes. Then:

> **C·(T⊗I) violates iff the frame a avoids an affine subspace A_M ⊆ F₃⁴ attached to the class.** (Good classes:
> A_M = F₃⁴.)

* So the frame enters only through one affine subspace, and a class is decided by the 5 frames 0, e₁, …, e₄.
* Every frame of every bad class was then decided.

| | count |
|---|---|
| bad symplectic classes | **6480 = 51 840 / 8** |
| …with A_M of codimension 1 (54 violating frames) | 5160 |
| …with A_M of codimension 2 (72 violating frames) | 24 |
| …with A_M empty (all 81 frames violate) | 1296 |
| violating Cliffords | **385 344** |
| **probability** | **223/2430 = 0.091770** (sampled in Pass 11266: 0.09124 ± 0.00064) |

**Structure.**
* **The bad fraction is exactly 1/8 for one qutrit (3/24) and for two (6480/51840).**
* All 648 classes fixing the magic axis (M z₁ = z₁) are bad.
* Every all-frames class satisfies ⟨z₁, M z₁⟩ = 0.
* The tested axis predicates do not capture the bad set. Its intrinsic description, and whether the 1/8 persists for
  n ≥ 3, are OPEN. The three-qutrit sample of Pass 11266 (0.095 ± 0.007) is consistent with it.
