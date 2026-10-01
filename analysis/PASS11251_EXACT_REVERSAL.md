# Pass 11251 — F_T((I⊗T)·SUM·(I⊗T²)) = F_min² exactly, and why: two quanta multiply

Producer: `analysis/w33_pass11251_exact_reversal.py`
Certificate: `data/w33_pass11251_exact_reversal.json`
Regression: `tests/test_w33_pass11250_11254.py`

**Question.** Pass 11238 found that the cheapest single-sector T-violator U = (I⊗T)·SUM·(I⊗T²) has best reversal
fidelity 0.7123860142010856. Pass 11237 matched this to F_min² = ((1 + 2cos 2π/9)/3)², but only to double precision.

## Proof (computer-assisted, exact)

1. **Exact arithmetic.** After removing a global phase, every qutrit Clifford has entries in 3^{−r/2}·μ₁₈, and
   T = diag(ζ^{x³}) has entries in μ₉ (ζ = e^{2πi/9}). The global phase cancels in V U* V†. So the trace
   t = tr(V U* V† U) is computed exactly in Z[ζ] = Z[x]/(x⁶ + x³ + 1).
2. **The attained value.** The exhaustive float search over all 51 840 × 81 anti-unitary two-qutrit Cliffords finds
   243 maximisers. All are monomial (r = 0), and 81 of them are tensor products. For one maximiser, exactly:
   * t = 3 + ζ − ζ² − ζ⁴ − 2ζ⁵;
   * t·t̄ = (1 + ζ + ζ⁸)⁴, which is the identity F_T² = F_min⁴ in Z[ζ].
3. **No larger value (Galois gap).**
   * Every two-qutrit Clifford has r ≤ 2: the smallest nonzero modulus is 1/3. So for every V, Y = 81|t|² is an
     algebraic integer of the cubic field Q(ζ)⁺.
   * Galois automorphisms commute with complex conjugation here, so every Galois conjugate of V is unitary up to
     3^{r/2}. All conjugates of Y therefore lie in [0, 81²].
   * Two distinct values of Y differ by at least 1/(2·81²)² = 5.8×10⁻⁹, which is 6.2×10⁻¹³ in F. That is far above
     the float error (≤ 10⁻¹³).
   * So the float maximum is the exact maximum, and **F_T = F_min² is proved.**

## The mechanism

* At a product maximiser, V U* V† U is **diagonal**.
* Its nine phases, in units of 2π/9, are {0, 0, 0, ±1, ±1, ±2}. That is the sumset {−1, 0, 1} + {−1, 0, 1}.
* So the trace is (1 + ζ + ζ⁻¹)²: the product of two copies of the one-qutrit minimal residual trace 1 + ζ + ζ⁻¹
  = 3F_min.
* The "two 2π/9 quanta composed" reading of Pass 11237 is therefore a theorem for this tick, not just a numerical
  coincidence.
* The next-best reversal of this tick has fidelity 0.666666667, numerically 2/3 (not proved exactly here).

**Scope.** This proves the identity for this one tick. Whether every F_min² level factorises in the same way is open.
