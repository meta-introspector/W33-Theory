# Pass 11370 — why the time-reversal witness has degree 10, 6, 4: a twisted Frobenius–Schur count

Producer: `analysis/w33_pass11370_twisted_indicator_degree_law.py`
Certificate: `data/w33_pass11370_twisted_indicator_degree_law.json`
Regression: `tests/test_w33_pass11369_11373.py`

**Question (from Pass 11357).** The first witness J_{2t} that is not identically zero has degree 10 for qubits, 6 for
qutrits and 4 for d = 5. Why these numbers?

## The reduction (exact)

* Pass 11355: J_{2t} ≡ 0 exactly when every Clifford-conjugation-invariant polynomial of degree (t,t) is unchanged by
  τ: f(U) ↦ f(Uᵀ).
* By Peter–Weyl, the degree-(t,t) polynomials on PU(d) are the matrix coefficients U ↦ tr(M λ(U)) of the irreps
  λ = (α, β) with |α| = |β| ≤ t (mixed Young diagrams).
* In a real basis λ(Uᵀ) = λ(U)ᵀ, so τ acts by M ↦ Mᵀ. The trace of M ↦ P Mᵀ Q on End(V) is tr(P Qᵀ). Hence:

> **#invariants in λ = (1/|Cl|) Σ_C |χ_λ(C)|²,  #τ-even − #τ-odd = (1/|Cl|) Σ_C χ_λ(C̄ C).**

* The second sum is the **twisted Frobenius–Schur indicator** for the automorphism C ↦ C̄ (Kawanaka–Matsuyama 1990).
  The Clifford group is closed under complex conjugation.
* Write λ restricted to Cl as Σ m_ρ ρ. Then the number of τ-odd invariants is
  n⁻_λ = ½ Σ_ρ m_ρ (m_ρ − ν_K(ρ)), with ν_K ∈ {−1, 0, 1}.

> **Theorem.** J_{2t} is not identically zero ⟺ some irrep λ with |α| = |β| ≤ t has n⁻_λ > 0. That happens exactly when
> λ, restricted to the Clifford group, repeats an irrep or contains one with twisted indicator ≠ 1.

## Qubits: the octahedral group explains the 10

* For SU(2), C̄ = Y C Y up to phase, with Y a Clifford. So the twisted indicator is the ordinary Frobenius–Schur
  indicator of the octahedral rotation group, whose irreps are all real.
* Hence n⁻ = Σ_ρ m_ρ(m_ρ − 1)/2: an odd invariant exists exactly when spin j restricted to O **repeats an irrep**.

| j | spin j restricted to O | n⁻ |
|---|---|---|
| 1 | T₁ | 0 |
| 2 | E + T₂ | 0 |
| 3 | A₂ + T₁ + T₂ | 0 |
| 4 | A₁ + E + T₁ + T₂ | 0 |
| **5** | **E + 2T₁ + T₂** | **1** |
| 6 | A₁ + A₂ + E + T₁ + 2T₂ | 1 |

* The first repeat is at j = 5, which first occurs in degree (5,5). So the witness degree is **10**.

## The table (computed from characters; every count is checked to be an integer)

| d | \|Cl\| | first t with an odd invariant | witness degree | odd invariants at that degree |
|---|---|---|---|---|
| 2 | 24 | 5 | **10** (Pass 11357: 10 ✓) | 1, in spin 5 |
| 3 | 216 | 3 | **6** (Pass 11357: 6 ✓) | 7: five in (3,3), one each in (4,1) and (1,4) |
| 5 | 3000 | 2 | **4** (Pass 11357: 4 ✓) | 3, in (2,0,0,2) |
| 7 | 16 464 | 2 | **4** (**predicted**) | 18 |

* **Prediction tested.** For d = 7 the count was computed first. Polynomial identity testing then gave J₂ ≤ 2×10⁻¹⁶
  and J₄ = 0.162 on Haar-random unitaries, so degree 4 is confirmed.

**Reading.**
* The dimension dependence of the arrow's witness is group theory: how complex conjugation acts on the irreps of the
  Clifford group inside the irreps of PU(d).
* Qubits need degree 10 because the octahedral group is "too real". Spin 5 is the first spin that repeats an
  octahedral irrep.
* Qutrits have exactly **7** independent T-odd invariants of the lowest degree. Pass 11369 shows that 7 is not
  enough to cut out the reversible set on all of PU(3).

**Prior art.** Twisted Frobenius–Schur indicators: Kawanaka–Matsuyama, Hokkaido Math. J. 19 (1990). Restricting
SO(3) irreps to the octahedral group is textbook crystal-field theory (Bethe 1929). In this repository, Pass 455 already
computed twisted indicators, for the order-27 Heisenberg groups under a centre-inverting involution; that is a different
group and a different twist. We found no statement linking twisted indicators to the degree of a time-reversal witness.
The search covered the corpus for "Frobenius-Schur", "twisted indicator", "Kawanaka" and "crystal field".
