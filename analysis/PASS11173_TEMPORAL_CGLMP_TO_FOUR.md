# Pass 11173 — one qudit's temporal Bell value tends to 4: an explicit strategy with I_d = 4 − 6/(d−1) − O(d^−3/2)

Producer: `analysis/w33_pass11173_temporal_cglmp_to_four.py`
Regression: `tests/test_w33_pass11173_temporal_cglmp_to_four.py`

**The ramp identity.** With d outcomes, every CGLMP term collapses to one linear ramp 1 − 2δ/(d−1). Here δ ∈ {0,…,d−1}
is the cyclic lag of Bob's outcome from the target relation:
* δ₀₀ = a − b, δ₁₀ = b − a − 1, δ₁₁ = a − b, δ₀₁ = b − a (mod d).

So

    I_d = 4 − (2/(d−1)) · L,   L = Σ_{x,y} E[δ_xy].

A local model has the four lags summing to −1 mod d, so L ≥ d − 1 and I_d ≤ 2. Approaching 4 means keeping the total
expected lag bounded as d grows. The identity is checked against direct CGLMP evaluation for every d, odd d included.

**The strategy,** read off the Pass 11171 optimisers.
* **Rest vector.** ψ = |0⟩ is the 0-th vector of A₀, of B₀ and of B₁. A₀ = B₀ = the standard basis.
* **Shift register.** On the complement of |0⟩, B₁ is B₀ shifted by one: B₁[a] = |a+1⟩ for a = 1…d−2, and B₁[d−1] = |1⟩.
* **Spreading.** A₁ = Q·B₁, where Q is the Householder reflection exchanging |0⟩ with the uniform vector.

With s = √d the lag is exactly

    L(d) = 3 + 2(s − 2)/(s(s − 1)),

checked with sympy for d = 4 … 64 (perfect squares) and in floating point to d = 1024, where I₁₀₂₄ = 3.99402. The 3
splits into about 1 from the spread vector A₁[0], about 1 from landing on the rest vector (cost d − 1 with probability
1/d² per outcome), and about 1 from the Householder spreading.

**Theorem.** I_d ≥ 4 − 6/(d−1) − O(d^−3/2). With Pass 11171 (I_d < 4 at every finite d), **sup_d I_d = 4, and it is
never attained.**

**The optimisers.** They do slightly better, L ≈ 1.7 (deficit ≈ 3.4/d), with the same architecture. At d = 8, p(a|0) is
0.978 on one vector shared by B₀ and B₁, and B₁ is B₀ shifted by one on its complement.

**Reading.** The qudit's only memory between the two times is its post-measurement vector. That vector still carries
the setting as well as the outcome: setting 0 parks it on the rest vector, setting 1 writes into the shift register.
The cost is O(1/d). This is why time can approach the algebraic maximum while space (CGLMP → 2.9696) cannot.
