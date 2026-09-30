# Pass 11186 — perfect five-qutrit gates beyond F₉: 17 892 more, the same two patterns; the determinant law holds at every order

Producer: `analysis/w33_pass11186_five_qutrit_beyond_f9.py`
Scan: `analysis/w33_pass11186_scan_ame10_tabu.py`, frozen as `data/w33_pass11186_ame10_tabu.json`
Regression: `tests/test_w33_pass11186_five_qutrit_beyond_f9.py`

**AME(10,3) graph states by tabu search.** Plain hill climbing found none (Pass 11175). Tabu search found **71 in 300
restarts**; the other 229 stalled at 8 singular blocks. With every choice of five input legs, the 71 graphs give
**17 892 perfect five-qutrit gates**, all verified symplectic and perfect.

**Patterns.** They realise exactly **the same two orientation patterns as the 2 642 411 520 F₉-linear gates**: the
10-cycle (5112) and a full row and column plus a permutation (12 780). Neither source produced the 4+6-cycle,
double-cross or all-(−1) pattern.
**Conjecture:** only these two orientation patterns occur for perfect five-qutrit gates.

**The determinant law of every order.** The symplectic form is a sum of local forms, ω = Σᵢ ωᵢ. Since S preserves
ωᵏ/k! = Σ_{|I|=k} ∧_{i∈I} ωᵢ, evaluating on the 2k-dimensional input space of a party set J gives

    Σ_{|I| = k} det S[I, J] = 1  (mod 3)   for every k and every set J of k input parties.

* k = 1 is Pass 11158's law, and k = n is det S = 1.
* Jacobi's identity with S⁻¹ = −JSᵀJ gives the complement symmetry det S[I,J] = det S[Iᶜ,Jᶜ] **exactly**.
* Both were verified on random symplectic matrices for 2 to 5 qutrits.
* These laws constrain the higher-order patterns, but on their own they do not exclude the three unrealised first-order
  patterns. That remains open.

**Prior art.** That the compound of a symplectic matrix preserves ωᵏ is standard linear algebra. The mod-3 sum law and
its use to constrain perfect-gate patterns are what this pass adds.
