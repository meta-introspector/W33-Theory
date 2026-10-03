# Pass 11330 — the one-magic-gate census made unconditional: 223/2430 over all 4 199 040 two-qutrit Cliffords

Producer: `analysis/w33_pass11330_orbit_census.py`
Certificate: `data/w33_pass11330_orbit_census.json` (+ `data/w33_pass11330_class_counts.npy`,
`data/w33_pass11330_bad_classes.npy`)
Regression: `tests/test_w33_pass11330_11334.py`

**The gap (Codex's audit, PASS11303_11307, publication-time integration correction).**
* Pass 11309 decided all 81 Pauli frames only for the 6480 classes flagged by the frames 0, e₁, …, e₄.
* The other 45 360 classes were certified only through the unproved affine law, so 223/2430 was conditional.

**Method: symmetry, not brute force.**
* Write every Clifford (mod phase) uniquely as C = W(a)·V_M, with V_M the canonical Weil unitary
  (V_M W(p) V_M† = W(Mp) exactly).
* Conjugation by any Clifford D with D T₁ D† ∝ T₁ preserves the verdict for U = C·T₁, because D U D† = (D C D†)·T₁.
* These D form a group K of order 27·72 = 1944: Paulis with no X₁ part, the Weil unitaries of the unit shears on
  qutrit 1, and Sp(2,3) on qutrit 2. Every generator was checked numerically to commute with T₁.
* K acts by (M, a) ↦ (N M N⁻¹, d + N a − N M N⁻¹ d). This formula was checked against matrix conjugation.
* The 4 199 040 elements fall into **4110 orbits** (sizes from 1 to 1944). One representative per orbit was decided
  exactly with the Pass 11252 criterion.

**Results.**

| | |
|---|---|
| orbit verdicts vs direct decisions on 1000 random elements | **1000/1000 agree** |
| violating Cliffords | **385 344** |
| probability | **223/2430** (exhaustive, unconditional) |
| classes by number of violating frames | 0: 45 360; 54: 5160; 72: 24; 81: 1296 |
| classes violating the affine law | **0 of 51 840** |
| five-frame screen 0, e₁…e₄ exact for every class | **yes** |
| agrees with Pass 11309 | yes |

**Status.**
* The affine law (C·T₁ violates iff the frame avoids an affine subspace A_M) is now verified on **every** class.
* The census is exhaustive without any shortcut.
* Pass 11331 confirms it by an independent algebraic method that never calls the Weyl-coefficient criterion.
