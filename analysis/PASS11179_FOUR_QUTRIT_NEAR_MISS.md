# Pass 11179 — how close four qutrits come to a perfect tick: two of 35 cuts fail, each by two trits, in a rigid paired pattern

Producer: `analysis/w33_pass11179_four_qutrit_near_miss.py`
Scan: `analysis/w33_pass11179_scan_four_qutrit_nearmiss.py`, frozen as `data/w33_pass11179_nearmiss_scan.json` and
`data/w33_pass11179_nearmiss_deficit_scan.json`
Regression: `tests/test_w33_pass11179_four_qutrit_near_miss.py`

**Background.** A perfect four-qutrit tick is impossible: its Choi state would be AME(8,3), which the shadow
inequalities exclude (Huber, Eltschka, Siewert and Gühne, arXiv:1708.06298; Pass 11170). So how close can one get? Up
to local Cliffords, every 8-qutrit stabilizer state is a weighted qutrit graph. Its balanced cut S|S^c is maximal
exactly when Γ[S, S^c] is invertible over F₃.

**Search (tabu).**
* Minimising the number of deficient balanced cuts: **120 of 120 restarts reach exactly 2**, never 1.
* Minimising the total deficit Σ(4 − rank): **36 of 36 restarts reach exactly 4**.

**Structure.** Every optimum found looks the same (20 checked):
* both deficient cuts have **rank 2**: two trits short, not one;
* the state is **3-uniform**: all 56 of the 3|5 cuts are maximal;
* the 8 qutrits split into **four pairs** P₀, P₁, P₂, P₃:
  * the deficient cuts are P₀+P₁ | P₂+P₃ and P₀+P₂ | P₁+P₃;
  * the third pairing, P₀+P₃ | P₁+P₂, is maximal.

**As a tick.** Take the maximal cut P₀+P₃ | P₁+P₂ as input | output. The result is a four-qutrit Clifford gate
(symplectic, checked) that is perfect except in one place:
* every single-qutrit block is invertible;
* the four **aligned pair-to-pair blocks** (input pair P₀ or P₃ to output pair P₁ or P₂) have rank 2 instead of 4, so
  each carries two trits instead of four.

**Reading.** The perfect-tick ladder is 1, 2, 3, 5 qutrits, never 4. What four qutrits miss is structured, not random
noise: the register pairs itself up, and information fails to spread fully only between aligned pairs.

**Scope.** The minimum of 2 (total deficit 4) is what the search reaches every time. It is a numerical statement, not a
proof that a single deficient cut is impossible; the shadow bound excludes only 0.
