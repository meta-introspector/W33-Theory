# Pass 11228 — measuring T-violation: the minimal violation is a 2π/9 phase, and Pauli polynomial witnesses fail

Producer: `analysis/w33_pass11228_t_witness.py`
Certificate: `data/w33_pass11228_t_witness.json`
Regression: `tests/test_w33_pass11227_11228.py`

**Question (from Pass 11213).** The cubic layer breaks substrate time-reversal symmetry: some ticks are inverted by no
anti-unitary Clifford Θ = V K. How can that be detected, and how large is it?

**1. The natural cubic witness vanishes identically.**
* J₃(U) = Im Σ_{P,Q} tr(U P U† Q)³ is invariant under Cliffords on both sides, and satisfies J₃(U⁻¹) = J₃(U) and
  J₃(U*) = −J₃(U). So it vanishes on every reversible tick.
* It also vanishes on every tick: max |J₃| = 3·10⁻¹³ over 240 words. Pairing each Pauli P with P† makes the sum real.
* Consequence: a T-witness built from the Pauli transfer matrix R(P,Q) must be **odd under exchanging rows and
  columns** (U ↔ U†). A reversible tick has U* and U† in the same two-sided Clifford coset.

**2. Transpose-odd magnitude moments are sound but incomplete.**
* Compare the sorted row and column moments of |R|⁴ and |R|⁶.
* They never fire on a reversible tick (0 of 4319), so a non-zero value certifies T-violation.
* But they detect only 3199 of 7781 violators (41%), and they miss the minimal violator T·X.
* So any complete witness must be phase-sensitive.

**3. A complete, operational measure: the best time-reversal fidelity.**

    F_T(U) = max over the 216 Cliffords V of |tr(V U* V† U)| / d.

* Here V U* V† U is unitary, and |tr| = d exactly when it is a scalar. So **F_T = 1 iff some substrate time reversal
  inverts U**. This was checked against the exact test on all 12 100 sampled words.
* 1 − F_T measures T-violation as an interferometric overlap: the best attempt to undo the dynamics by a time
  reversal.
* **All 18 minimal violators (the cubic phase times a cyclic shift) have F_T = (1 + 2cos 2π/9)/3 = |1 + ζ₉ + ζ₉⁻¹|/3
  = 0.8440296287…** The irreducible T-odd content of the minimal violating tick is a **phase of 2π/9 = 40°**,
  quantised by the cubic gate.
* Deeper circuits reach lower fidelities (minimum 0.7257876540… by depth 3). Mean violator fidelity by depth 1–8:
  0.844, 0.811, 0.867, 0.889, 0.904, 0.910, 0.913, 0.919. As depth grows, violations become more common (Pass 11213)
  but individually milder.

**Reading.** On one qutrit, neither the natural cubic Pauli invariant nor magnitude moments of the transfer matrix
detect T-violation completely. Whether some phase-sensitive polynomial does is left open. The failure of
time-reversal interferometry, 1 − F_T > 0, detects it exactly. Its minimal quantum is
the phase 2π/9 carried by ζ₉^{x³}. No identification with a measured CP phase is made or implied (δ_CKM ≈ 65–70°).
