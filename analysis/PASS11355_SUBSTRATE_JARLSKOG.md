# Pass 11355 — a substrate Jarlskog invariant: one degree-six formula detects every one-qutrit T-violation checked

Producer: `analysis/w33_pass11355_substrate_jarlskog.py`
Certificate: `data/w33_pass11355_substrate_jarlskog.json`
Regression: `tests/test_w33_pass11355_11360.py`

## Theorem 1 (exact)

> A tick U is substrate-time-reversible **iff its transpose is Clifford-conjugate to it** (up to a phase).

Proof: V Ū V† = λU⁻¹ ⇔ Ū = λ V† U⁻¹ V ⇔ Uᵀ = Ū⁻¹ = λ̄ V† U V.

**Corollary.** Every phase-invariant, Clifford-conjugation-invariant function f satisfies f(Uᵀ) = f(U) on reversible
ticks, so f(U) ≠ f(Uᵀ) **witnesses** T-violation.

## The witness in closed form

* The degree-(t,t) invariants are tr(X·π_t(U)), with π_t(U) = U^{⊗t} ⊗ Ū^{⊗t} and X in the Clifford commutant.
* All of them are captured by A_t(U) = avg_C π_t(CUC†), the orthogonal projection of π_t(U) onto the commutant.
* Because ⟨π_t(X), π_t(Y)⟩ = |tr(X†Y)|^{2t}:

> **J_{2t}(U) = avg_C |tr(U†·CUC†)|^{2t} − avg_C |tr(U†·CUᵀC†)|^{2t} = ½‖A_t(U) − A_t(Uᵀ)‖² ≥ 0.**

* It is an average of 216 one-qutrit traces. The identity was checked on explicit 729×729 matrices to 4×10⁻¹⁴.

## Theorem 2 (corrected) — the lowest qutrit witness has degree (3,3)

* **J₂ ≡ 0 for any unitary 2-design.** The commutant of C ⊗ C̄ is spanned by the identity and the maximally entangled
  projector, and both are transpose-symmetric.
* **J₄ ≡ 0 for qutrits, as a verified polynomial identity.** It is not a consequence of a design property.
* **Correction.** An earlier draft "derived" J₄ = 0 from the 2-design property. That derivation is wrong:
  * the degree-(t,t) invariants come from the commutant of C^{⊗t} ⊗ C̄^{⊗t}, whose dimension is the frame potential at
    2t;
  * for d = 5 the Clifford group is a 2-design, yet J₄ is not identically zero.

  See Pass 11357, which decides vanishing by polynomial identity testing on Haar-random unitaries.
* So for qutrits the lowest witness degree is **(3,3)**. It is cubic, like the quark-sector Jarlskog invariant
  Im tr[H_u, H_d]³.

## Finding (computer-verified): J₆ is complete on every one-qutrit word checked

| cubic gates | words (coset-reduced, Pass 11312) | J₆ > 0 | exact violators | per-word check against the Weyl criterion |
|---|---|---|---|---|
| 1 | 216 | 18 | 18 | 0 mismatches |
| 2 | 5184 | 2106 | 2106 | 0 mismatches |
| 3 | 124 416 | 68 202 | 68 202 | 0 mismatches |
| 4 | 2 985 984 | **2 077 650** | **2 077 650** | count equality |

* The smallest positive values are **9/8** (one gate, numerically 1.124999999999961), 0.368 (two gates), 0.0198
  (three) and 1.5×10⁻⁴ (four). All are far above the 10⁻⁹ tolerance.
* **Scope.** "Complete" is verified on these words, not proved for all one-qutrit unitaries. The invariant ring of a
  finite group separates orbits, so *some* degree always suffices. That degree 6 already suffices on everything checked
  is the finding.
