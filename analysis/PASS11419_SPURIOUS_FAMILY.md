# Pass 11419 — anatomy of J₆'s blind spot: two 4-dimensional components, one of them explained by a unique T-odd ray invariant of degree 6

Producer: `analysis/w33_pass11419_spurious_family.py` (main run; `--fact` for the identities and counts)
Certificate: `data/w33_pass11419_spurious_family.json`
Regression: `tests/test_w33_pass11418_11422.py`

**Background.** Pass 11369 found non-reversible zeros of J₆ on PU(3) in two populations. About 65% have a generic
spectrum (the refined point is smooth with rank dF = 4); about 35% have a nearly degenerate spectrum.

## The pseudo-reflection component (theorem, given one computed count)

* **Setup.** Up to phase, a unitary with a repeated eigenvalue is a pseudo-reflection U = I + λP, with P = |ψ⟩⟨ψ|
  and λ = e^{iβ} − 1.
* **Reversibility.** Uᵀ = I + λP̄, so U is reversible iff ψ̄ lies in the Clifford orbit of ψ. Call such rays
  "real-type"; they form a 2-dimensional set in CP².
* **The traces are real.**
  * tr(U†CUC†) = 3 − |λ|²(1 − x), with x = |⟨ψ|Cψ⟩|².
  * tr(U†CUᵀC†) = 3 − |λ|²(1 − y), with y = |⟨ψ|Cψ̄⟩|².
* **Hence the expansion.** Writing Δ_k = avg_C x^k − avg_C y^k:

> J₆(U) = Σ_k C(6,k) 3^{6−k} (−|λ|²)^k Δ_k.

* **Which Δ_k survive.**
  * avg x^k = ⟨P^{⊗k}, B_k(P)⟩, with B_k(P) = avg_C (CPC†)^{⊗k}.
  * The coefficients of B_k(P) are the Clifford-invariant polynomials of bidegree (k,k) in ψ, and B_k(Pᵀ) uses them
    at ψ̄.
  * **Computed: τ-odd invariants first appear at k = 6, with exactly one there.**

| k | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| Clifford invariants of bidegree (k,k) | 1 | 1 | 2 | 3 | 4 | 7 | 10 |
| τ-odd | 0 | 0 | 0 | 0 | 0 | **1** | 2 |

  * So Δ₁ = ⋯ = Δ₅ = 0 (measured: ≤ 10⁻¹⁶) and Δ₆ = 2|h₆(ψ)|², where h₆ is the unique T-odd Clifford-invariant
    polynomial of bidegree (6,6).
  * Δ₁–Δ₄ = 0 also follow independently from J₄ ≡ 0.

> **J₆(I + λ|ψ⟩⟨ψ|) = |λ|¹² · Δ₆(ψ) = |λ|¹² · (avg_C |⟨ψ|Cψ⟩|¹² − avg_C |⟨ψ|Cψ̄⟩|¹²).**

* **Checked numerically.**
  * The β-dependence is exactly |1 − e^{iβ}|¹² (deviation 1.7×10⁻⁹).
  * The matrix J₆[ray, β] has rank 1: its second singular value is 4×10⁻¹³ of the first.
  * J₆(U)/|λ|¹² = J₆(P) to rounding level (1.3×10⁻⁷ of the scale; absolute errors about 10⁻¹³).
  * J₆(P) equals the moment difference Δ₆ to 7×10⁻¹⁰.
* **Consequence.** The zero set {h₆ = 0} is a real hypersurface in CP²:
  * 200 minimisations gave 182 zeros, 125 of them on non-real-type rays (Clifford overlap down to 0.960);
  * the Hessian has rank 1 at every zero tested.
* With β free, {I + λ|ψ⟩⟨ψ| : h₆(ψ) = 0} is a **4-dimensional spurious component** of J₆ (3 + 1 dimensions). It sits
  inside the 5-dimensional degenerate-spectrum locus, while the reversible part of that locus is only 3-dimensional.

**The near-degenerate descents land on it.**
* All 17 near-degenerate certified spurious points from fresh descents (of 49 certified) snap onto the component:
  minimising Δ₆ from their isolated eigenvector reaches zero after moving the ray by at most 8×10⁻⁹.
* 16 of the 17 are off real-type.
* The other 32 points have a generic spectrum (smallest gap 0.040·2π). They form the second component, whose
  structure stays open.

**A wrong guess, recorded.** The natural cubic Clifford invariant S₃(ψ) = Σ_p ⟨ψ|W(p)|ψ⟩³ is **real** for every ray
(imaginary part ≤ 4×10⁻¹⁵), hence τ-even. It is not the odd invariant: that one first appears at degree 6.

**Why J₈ does not have this component** (see Pass 11418).
* On pseudo-reflections, J₈ = |λ|¹²(252 Δ₆ − 24|λ|² Δ₇ + |λ|⁴ Δ₈). It is not rank 1 in (ray, β).
* Its zeros at a fixed β need a combination of three moment gaps to vanish, and degree 7 already carries two odd
  invariants.
* Every J₈ zero found on the locus is real-type.

**Open.** An explicit formula for h₆; the structure of the generic-spectrum component; whether any Clifford+T word
lies exactly on either component.
