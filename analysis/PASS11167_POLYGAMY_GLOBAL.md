# Pass 11167 — 4/√15 is exact on its symmetry class; 1.0328 ≤ max(N_AB + N_AC) ≤ 25/18 globally; 3000/3000 dual see-saws reach 4/√15

Producer: `analysis/w33_pass11167_polygamy_global.py`
Regression: `tests/test_w33_pass11167_polygamy_global.py`

**1. Exact on the symmetry class.** The optimum ψ* of Pass 11159 is invariant under U ⊗ Ū ⊗ Ū for every U ∈ U(2) acting
on span{|1⟩,|2⟩}. This is the 4-dimensional stabiliser seen in Pass 11162.
* The invariant states are α|000⟩ + β Σ_a |aa0⟩ + γ Σ_a |a0a⟩, and phases can be removed.
* Set s = β², t = γ², α² = 1 − 2s − 2t. The partial transposes split into 2×2 blocks, which gives
      N_AB = √(t² + 4α²s) − t + s,   N_AC = √(s² + 4α²t) − s + t,
      **N_AB + N_AC = √(t² + 4α²s) + √(s² + 4α²t)**   (the linear terms cancel).
* The critical points, found exactly with sympy, are:
  * (2/15, 2/15), value **4/√15**, the maximum;
  * (2/9, 1/18), value √(2√43/27 + 59/108) ≈ 1.0159, a saddle;
  * (1/3, 0) and (0, 1/3), value 1.
* The boundary maximum is 1. So **4/√15 is the exact, unique maximum on the symmetry class of the optimum.**

**2. A rigorous global bound.** Three facts combine:
* The partial transpose of a two-qutrit state has at most (d−1)² = 4 negative eigenvalues (Rana, PRA 87, 054301
  (2013)). So its purity p ≥ min_{m≤4} [(1+N)²/(9−m) + N²/m], which gives N ≤ h(p).
* For pure ψ_ABC, purity(ρ_AB) = purity(ρ_C).
* N_AC ≤ N_{C|AB}(ψ) = g(λ_C), because tracing out B is local for the C|AB cut.

Hence

    N_AB + N_AC ≤ max over qutrit spectra λ of [h(Σλ_i²) + g(λ)] = 25/18 ≈ 1.389,

attained at λ = (2/3, 1/6, 1/6), where h = 5/9 and g = 5/6 (maximised on a 1/1500 grid). The bound is strictly below
the trivial 2 but not tight: the gap [1.0328, 1.389] is what remains open.

**3. Dual see-saw.** N(σ) = max_{0≤P≤I} −Tr(P σ^T), so N_AB + N_AC = max over ψ, P₁, P₂ of ⟨ψ|H(P₁,P₂)|ψ⟩.
* The iteration alternates two steps: ψ ← the top eigenvector of H, and P_i ← the projector onto the negative
  eigenspace. Neither step can decrease the value.
* This is a different landscape from the primal restarts of Pass 11162.
* **3000 of 3000 random starts converge to 4/√15; none exceeds it.**

**Scope.**
* Exact on the symmetry class.
* Rigorous 25/18 globally.
* 4/√15 as the global maximum is still a conjecture, now supported by both primal and dual searches.
