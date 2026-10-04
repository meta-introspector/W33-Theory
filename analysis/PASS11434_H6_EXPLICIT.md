# Pass 11434 — the T-odd invariant of a qutrit state is a chirality of the MUB geometry: h₆ = signed sum of p_a³ p_b² p_c over oriented triples of MUBs

Producer: `analysis/w33_pass11434_h6_explicit.py`
Certificate: `data/w33_pass11434_h6_explicit.json`
Regression: `tests/test_w33_pass11433_11437.py`

**Background.**
* Pass 11419 proved that the first T-odd Clifford-invariant polynomial of a qutrit ray has bidegree (6,6), and that it
  is unique there.
* Call it h₆. Then J₆(1 + λ|ψ⟩⟨ψ|) = |λ|¹² Δ₆(ψ) with Δ₆ = 2|h₆|².
* Two natural guesses are wrong: Σ_p χ_p³ and Σ_p χ_p⁶ (χ_p = ⟨ψ|W(p)|ψ⟩). Both are **real** for every ray, hence
  T-even, because χ_{−p} = χ̄_p pairs the terms.

## The explicit form (computed, exact up to rounding)

* **Coordinates.** The 12 stabiliser states s of the four MUBs (Z, X, XZ, XZ²) determine a ray through
  p_s(ψ) = |⟨s|ψ⟩|², and the extended Clifford group permutes them. That group has 216 unitaries with sign +1 and 216
  antiunitaries with sign −1, i.e. 432 distinct permutations of the 12 states.
* **The search.** For each degree-6 monomial m in the p_s, the signed average
  R_m(ψ) = (1/432) Σ_g sgn(g) m(gψ)
  is T-odd and Clifford-invariant of bidegree (6,6), hence a multiple of h₆. The degree-6 monomials fall into **74**
  orbits. Exactly **8** give a nonzero R_m.
* **All 8 have the same shape:** exponents 3, 2, 1 on three stabiliser states from three **different** MUBs. All are
  proportional to h₆; the constant ratio Δ₆ / R_m² is 349 920 or 4 × 349 920.

> **h₆(ψ) = (1/432) Σ_{g ∈ ext. Clifford} sgn(g) · p_{g·a}(ψ)³ p_{g·b}(ψ)² p_{g·c}(ψ),**
> a, b, c stabiliser states of three distinct MUBs, and **Δ₆ = 349 920 · h₆²** (= 2⁵·3⁷·5 · h₆²).

**Reading.**
* Reversing time conjugates ψ, which reverses the orientation of the MUB triple. h₆ counts the weight of a state on
  oriented MUB triples with an orientation sign: it is a **chirality** of the state relative to the four mutually
  unbiased bases.
* The pseudo-reflection component of J₆'s blind spot is exactly the achiral hypersurface {h₆ = 0}, minus the
  real-type rays.
* **Heuristic for degree 6, not a proof.** The Clifford group acts on the four MUBs as A₄. The antiunitaries act by
  odd permutations; K fixes the Z and X bases and swaps XZ ↔ XZ².
  * A₄ is 2-transitive, so ordered pairs of MUBs carry no orientation.
  * Ordered triples split into two A₄-orbits, which are chiral.
  * Telling a triple from its reversal with one monomial needs three distinct exponents, the smallest being
    1 + 2 + 3 = 6.
  * The proof that nothing odd exists below degree 6 is the invariant count of Pass 11419, not this argument.

## Generic component (observation)

| | generic J₆-spurious points (24) | Haar-random unitaries (300) |
|---|---|---|
| some eigenvector with Δ₆ < 10⁻⁸ | **79%** | 13% |
| median of the smallest eigenvector Δ₆ | **3.3×10⁻¹⁰** | 5.8×10⁻⁷ |

* Generic-spectrum spurious points tend strongly to have an **achiral eigenvector**: one whose h₆ ≈ 0.
* That condition alone defines a 7-dimensional set, while the component is 4-dimensional, so further conditions must
  hold.
* The characterisation of the generic component remains open.

**Prior art.** MUBs and stabiliser states of a qutrit are standard (Wootters–Fields 1989; Gross 2006). We found no
earlier identification of the lowest T-odd Clifford invariant of a ray with a chirality of MUB triples. The search
covered the corpus for "chirality", "MUB", "oriented triple" and "349920".
* The integer 349 920 does occur in `analysis/2026-05-18_minimal_logical_dual_visibility_scheme.md`, as the count
  1620 · 432 / 2. There, 432 equals the order of the extended qutrit Clifford group used here. The match is recorded,
  not explained.
* "Oriented triples" in Pass 11222 (product Lagrangians) and in the 2026-05-30 Fano-triangle note are unrelated
  objects.
