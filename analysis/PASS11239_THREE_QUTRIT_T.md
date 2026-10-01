# Pass 11239 — three qutrits: the minimal T-violation embeds, and generic search cannot see three-qutrit reversals

Producer: `analysis/w33_pass11239_three_qutrit_t.py`
Certificate: `data/w33_pass11239_three_qutrit_t.json`
Regression: `tests/test_w33_pass11236_11239.py`

**Question.** The two-qutrit statements (Passes 11227, 11235, 11238) are exact: they search every anti-unitary two-qutrit
Clifford. Does the 2π/9 minimal violation, F = (1 + 2cos 2π/9)/3, persist on three qutrits?

**Why only lower bounds.**
* The three-qutrit anti-unitary Clifford group has |Sp(6,3)| × 729 ≈ 3.3 × 10¹² elements. Exhaustive search is out of
  reach.
* Magnitude-only certificates of violation are incomplete: they miss the minimal violator (Pass 11228).
* So every number below is a lower bound on the best reversal fidelity F_T = max_V |tr(V U* V† U)|/27. It certifies
  reversibility when it reaches 1 and certifies nothing when it does not.

**Method.**
* *Product reversals (exact for products).* If U = A ⊗ B, then V = V_A ⊗ V_B gives a trace that factorises, so
  F(A ⊗ B) ≥ F(A)·F(B). A diagonal one-qutrit tick has F = 1 with V = I. The two-qutrit factors use the exhaustive optima.
* *Generic search.* 3000 random three-qutrit Cliffords (words of length 40 in H, S, SUM), with the Pauli part maximised
  exactly over all 729 Paulis. Then greedy hill-climbing from the six best.

**Results.**

| tick | factorisation | product reversal | generic search |
|---|---|---|---|
| [(T⊗T)·SUM] ⊗ I | minimal violator ⊗ I | **0.8440296** = F_min | 0.601 |
| (T⊗T⊗T)·SUM₁₂ | [(T⊗T)·SUM] ⊗ T | **0.8440296** | 0.458 |
| (I⊗I⊗T)·SUM₂₃ | I ⊗ [(I⊗T)·SUM] | **1 (reversible, certified)** | 0.411 |
| (T⊗T⊗T)·SUM₂₃·SUM₁₂ | none | — | 0.237 |
| (T⊗I⊗T)·SUM₂₃·SUM₁₂ | none | — | 0.200 |

**What this shows.**
* **The minimal violation embeds.** A two-qutrit violator tensored with an idle or diagonally ticking third qutrit keeps
  a reversal at the 2π/9 level. Whether a genuinely three-qutrit reversal does better (F → 1) is **open**.
* **The generic search is blind.** On a tick that is exactly reversible (row 3, F = 1), it reaches only 0.41. On the
  product case (row 2) it reaches 0.46 against the exact 0.844. The low numbers for the two entangled circuits (rows 4,
  5) are therefore **no evidence** that they violate T, or of how strongly.
* **Over-read avoided.** A first reading of rows 4 and 5 as "three-qutrit violation is stronger" would have been
  wrong. Row 3 is the control that shows it: the search finds 0.41 where the truth is 1.

**Open.** A sound three-qutrit statement needs either the orbit structure of Sp(6,3) acting on the reversal problem
(an analogue of the 51 840-element enumeration), or a complete invariant that is not magnitude-only (Pass 11228).
