# Pass 11183 — an arrow of time that no choice of subsystems removes: quantised, never less than two trits

Producer: `analysis/w33_pass11183_intrinsic_arrow.py`
Regression: `tests/test_w33_pass11183_intrinsic_arrow.py`

**Definition.** Theorem 4.7 and Pass 11161: a tick followed by forgetting the partner loses whatever the tick exported
about a qutrit's past. For an n-qutrit Clifford tick S and a split F (Pass 11180), the entropy formula gives
I(q_in : q_out) = 2 − rank S_F[others, q] (checked against Pass 11165). So forgetting everything but q destroys
export_q = rank S_F[others, q] trits.
* The **arrow strength** of S in F is E(S,F) = Σ_q export_q.
* The **intrinsic arrow** is A(S) = min over all splits of E(S,F): the loss that no choice of subsystems avoids.
* A(S) = 0 iff S is local and keeps every subsystem in place in some split. A swap moves a whole past into the partner.

**Two qutrits (exact, all of Sp(4,3)).**
* **A = 0 for 19 152 ticks and exactly 2 trits for 32 688 (227/360). It is never 1.** Symplecticity forces
  rank S_BA = rank S_AB, so the exported information is even.
* A = 2 for every tick of order 5 or 9. It is also 2 for the order-4 and order-6 ticks whose only local splits swap the
  qutrits.

**Three qutrits** (a sample of 240, plus the paper's ticks).
* **A ∈ {0, 2, 3}, never 1.** If two qutrits' columns are local and non-permuting, symplectic orthogonality makes the
  third plane invariant too. So information never leaks out of a single subsystem alone.
* In the sample, A = 3 occurs at orders 3, 6, 7, 9, 13, 14, and A = 2 at orders 5, 6, 8, 9, 10, 12, 15, 18, 20, 24, 30,
  36. A = 0 occurs only at orders 6 and 12.

**The paper's ticks.**

| tick | A | where |
|---|---|---|
| clock (Theorem 4.3) | **0** | its 108 local splits are all arrow-free |
| perfect interacting tick V K V | **2** | attained in 162 splits |
| F₉ gate K = P² + i𝟙 | **3** | local in 27 splits but cycles the qutrits in every one |

**Reading.** An objective arrow of time, present in every subsystem description, exists exactly for these ticks. It is
quantised: its smallest nonzero value is two trits, one full qutrit record, which is the unit of Theorem 4.7.
