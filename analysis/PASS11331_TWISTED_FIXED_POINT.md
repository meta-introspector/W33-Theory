# Pass 11331 — one magic gate is a twisted fixed-point problem on the stabilizer of the magic axis

Producer: `analysis/w33_pass11331_twisted_fixed_point.py`
Certificate: `data/w33_pass11331_twisted_fixed_point.json`
Regression: `tests/test_w33_pass11330_11334.py`

## Exact reformulation

Let U = C·T₁ with C Clifford and T₁ = T⊗I. An anti-unitary Clifford reversal D·K exists iff D U* D† = λU⁻¹. Put
E = D·C̄. Then

> **U is substrate-time-reversible ⇔ some Clifford E solves T₁ E T₁⁻¹ = λ·C⁻¹ E Cᵀ.**  (*)

Conversely, a solution gives the reversal D = E·C̄⁻¹. The right side of (*) is Clifford, so E must lie in
N_T = {E Clifford : T₁ E T₁⁻¹ Clifford}.

## Computed facts

| | |
|---|---|
| symplectic parts of N_T | **exactly Stab(z₁)**, the stabilizer of the magic axis (648 classes; one qutrit: the 3 unit shears) |
| Pauli frames in N_T | all 81 (T₁ is in the third Clifford-hierarchy level) |
| distinct symplectic parts of T₁ E T₁⁻¹ | 1944 = 648 × 3: they depend only on Q and on the X₁ component of e |
| census by (*), one representative per Pass 11330 orbit | **385 344 violators, identical to Pass 11330 for every one of the 51 840 classes** |

The second derivation never calls the Weyl-coefficient criterion of Pass 11252. It is an independent confirmation of
223/2430.

## Where the obstructions sit

* **1296 classes:** (*) has no solution even at the symplectic level. These are exactly the classes in which every
  Pauli frame violates.
* **5184 classes (5160 + 24):** (*) is solvable symplectically, but the Pauli-frame part fails off an affine subspace
  of codimension 1 or 2.
* **45 360 classes:** solvable at both levels for every frame.

**Mechanism of the affine law (sketch, not a proof).** For a fixed symplectic solution Q, the Pauli-level part of (*)
is linear over F₃ in the frames (a, e). So the frames admitting a solution form an affine set. A union over several Q
need not be affine; Pass 11330 verifies the law on every class.

## The bad set is W(3,3) geometry of the magic axis

The magic axis z₁ is a point of W(3,3). Classify each symplectic class M by where it sends z₁ (and z₁'s preimage):

| where M moves the magic axis | classes | bad |
|---|---|---|
| fixes it: Mz₁ = z₁ (= the symplectic part of N_T) | 648 | **648 (all)** |
| reverses it: Mz₁ = −z₁ | 648 | **0** |
| to a collinear point, M⁻¹z₁ on a different line through z₁ | 11 664 | **0** |
| along one line through z₁, M²z₁ = z₁ | 648 | **648 (all)** |
| along one line through z₁, M²z₁ = −z₁ | 648 | **0** |
| along one line through z₁, otherwise | 2592 | **1296**: exactly the symplectically unsolvable, all-81-frame classes |
| to a non-collinear point | 34 992 | **3888 = 1/9** |

* Total: 648 + 648 + 1296 + 3888 = 6480 = |Sp(4,3)|/8.
* So the 1/8 decomposes over the collinearity graph of W(3,3) around the magic axis.
* Inside the non-collinear cell, the 1/9 is uniform across ω(z₁, Mz₁), ω(Mz₁, M⁻¹z₁) and the position of M²z₁. Its
  finer description is open.

**Open.** A closed count of the 1296 symplectically unsolvable classes (1/40 of Sp(4,3)) and of the 5184
frame-obstructed ones. Together they give the bad fraction 1/8, which equals the one-qutrit value 3/24.
