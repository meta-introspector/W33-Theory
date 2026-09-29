# Pass 11169 — six local orbits of perfect three-qutrit gates, labelled by which future receives each past time-reversed; unique up to relabelling

Producer: `analysis/w33_pass11169_perfect_gate_local_orbits.py`
Regression: `tests/test_w33_pass11169_perfect_gate_local_orbits.py`

**The invariant π.** For a perfect S ∈ Sp(6,3), every column of block determinants is (1, 1, −1) (Pass 11158), and so is
every row. So the −1's form a permutation matrix: det S_ij = −1 exactly when i = π(j).
* A block of determinant −1 reverses the orientation of the qutrit phase space: it is a local time reversal, as in
  Pass 11156's anti-unitary law.
* So **π ∈ S₃ says in which future each past appears time-reversed.**
* Local Clifford operations S ↦ LSR, with L, R ∈ SL(2,3)³, preserve every block determinant, so π is a local invariant.

**Orbits.**

| | value |
|---|---|
| stabiliser order of a representative of each of the 6 classes | 4 |
| orbit size under SL(2,3)⁶ (order 24⁶ = 191 102 976) | 47 775 744 |
| class size = 286 654 464 / 6 (classes are equal because output permutations permute them) | 47 775 744 |

* **Each π-class is one local orbit, so π is a complete local invariant.**
* Relabelling parties (π ↦ σπτ⁻¹) joins the six classes. Up to local operations and relabelling, **the perfect
  three-qutrit gate is unique**, as the perfect two-qutrit gate is: stabiliser 24, orbit 13 824, recomputed here.
* The F₉-linear perfect gates (Pass 11164) are spread evenly: **2048 in each class** (12 288 in total).

**Prior art.** Tan (arXiv:2601.19677) proved local-unitary uniqueness of AME(4,3). The three-qutrit statement here is
Clifford-level: local Clifford orbits of perfect Clifford gates. We did not find the π-labelling in the literature.
