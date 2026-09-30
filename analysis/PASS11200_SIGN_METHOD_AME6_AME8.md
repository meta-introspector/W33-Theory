# Pass 11200 — the sign method across the AME(2m,3) family: projective lines over F_{2m-1}, and a hypothetical AG(3,2)

Producer: `analysis/w33_pass11200_sign_method_ame6_ame8.py`
Certificate: `data/w33_pass11200_sign_method_ame6_ame8.json`
Regression: `tests/test_w33_pass11200_sign_method_ame6_ame8.py`

**Method (Pass 11190, general m).**
* In an AME(2m,3) stabilizer state, for an (m−1)-set K the 2-dim space of stabilisers trivial on K restricts
  isomorphically to each of the m+1 other parties. The symplectic signs sum to 0 mod 3.
* The possible splits of the (m+1)-sets are therefore 3|0 (m = 2), 2|2 (m = 3), 4|1 (m = 4) and 3|3 or 6|0 (m = 5).
* Only m = 2 and m = 5 allow uniform sets.

**AME(6,3), three-qutrit perfect gates.**
* Every local graph (on the 5 parties outside one party) is a **pentagon**, in 72 of 72 cases over 12 random AME(6,3)
  graph states.
* The sign structure is unique: all 12 are isomorphic.
* Its automorphism group has order **120**, with 25 involutions, 30 elements of order 4, 24 of order 5 and 20 each of
  orders 3 and 6. That is S5 ≅ **PGL(2,5)** acting on the **projective line over F5**.
* This is the three-qutrit analogue of Pass 11190's PGL(2,9) on PG(1,9). It explains the permutation patterns of
  perfect three-qutrit gates, where each 4-set pairs a party with one partner.

**AME(8,3), four qutrits.**
* AME(8,3) does not exist, by the shadow inequalities (Huber, Eltschka, Siewert and Gühne, arXiv:1708.06298; Pass
  11170).
* The sign method forces every local graph on 6 parties to be a **perfect matching** (degrees in {1,4}).
* Even so, CP-SAT with the isotropy duality is **feasible**. It has exactly **30** solutions, all isomorphic, each with
  automorphism group of order **1344 = |AGL(3,2)|**.
* In every solution, removing the singleton of each 5-set leaves the unique block it contains in a Steiner system
  **S(3,4,8)**, the affine planes of AG(3,2).
* So a hypothetical AME(8,3) stabilizer state would carry the affine geometry AG(3,2). Its non-existence is **not** a
  sign obstruction: the sign method is consistent there, and the shadow inequalities see more.

**The family.**

| m | parties | sign structure | symmetry |
|---|---|---|---|
| 2 | 4 | trivial | S4 = PGL(2,3) on PG(1,3) |
| 3 | 6 | pentagons | PGL(2,5) on PG(1,5) |
| 4 | 8 | S(3,4,8) = AG(3,2) (not realised) | AGL(3,2) |
| 5 | 10 | S(3,4,10) | PGL(2,9) on PG(1,9) (Pass 11190) |

**Scope.** Stabilizer AME states. The AME(6,3) uniqueness rests on 12 random samples and the automorphism
computation. It is not a proof over all AME(6,3) states.
