# Pass 11308 — the shared cyclotomic level is |tr T/3|²; closeness to Clifford bounds time-reversal violation

Producer: `analysis/w33_pass11308_shared_level_and_clifford_bound.py`
Certificate: `data/w33_pass11308_shared_level_and_clifford_bound.json`
Regression: `tests/test_w33_pass11308_11312.py`

**The open coincidence.** Codex's balanced-G26 note (`PASS20261001_BALANCED_G26_TMAGIC_VACUUM.md`) records that
F_min² = ((1 + 2cos 2π/9)/3)² = 0.712386… is both:
* the largest qutrit-stabilizer probability of the T-magic vacuum ray |T3⟩ = T|+⟩; and
* the exact two-qutrit reversal-fidelity level of Passes 11238/11251.

It calls this "a shared algebraic invariant, not yet a derived physical identity".

## Derivation (exact)

T = diag(ζ^{x³}) has spectrum ζ^{0, 1, −1}.
1. **Vacuum side.** ⟨+|T|+⟩ = tr T/3 = (1 + ζ + ζ⁻¹)/3 = F_min. So |⟨+|T3⟩|² = |tr T/3|², and |+⟩ attains the
   maximum over all 12 stabilizer states.
2. **Dynamics side.** The optimal residual of (I⊗T)·SUM·(I⊗T²) has the spectrum of T⊗T: the sumset {−1,0,1}+{−1,0,1}
   of Pass 11251, checked exactly in Z[ζ]. So its reversal fidelity is |tr(T⊗T)|/9 = |tr T/3|².
3. Both numbers are **the squared normalised trace of the cubic gate**. The vacuum number is that trace seen on a
   stabilizer state; the dynamical number is the same trace carried twice by the optimal residual.

## Clifford fidelity bounds T-violation

Let F_Cl(U) = max_C |tr(C†U)|/d over Cliffords C.

* **Theorem (any n): F_T(U) ≥ 2F_Cl(U)² − 1.**
  * Write U = C·R with C the nearest Clifford. C has an anti-unitary Clifford reversal V C* V† = λC† (no Clifford
    T-violation).
  * Then tr(V U* V† U) = λ·tr(A R) with A = C†V R* V†C and |tr A| = |tr R| = d·F_Cl.
  * Splitting off the traceless parts and applying Cauchy–Schwarz gives the bound.
* **Numerics** (one-qutrit stream of Pass 11237, 1344 violating words; two-qutrit stream of Pass 11235, 159 words):
  * **F_T ≥ F_Cl² holds on every word.** The minimum margin is 0.0134 on one qutrit and 0.096 on two. Unproved.
  * **F_T ≥ F_Cl is false.** There are counterexamples on one qutrit and three on two, e.g. (F_T, F_Cl) =
    (0.594859, 0.601275), (0.641015, 0.712386), (0.808257, 0.844030).
  * A first look at 60 two-qutrit words had suggested F_T ≥ F_Cl; the full sample refutes it.
* **Reading.** A tick whose magic is a small perturbation of a Clifford cannot break substrate time reversal by much.
  A large T-violation needs a tick far from every Clifford.
