# Pass 11491 — how many ways a qutrit state can know the direction of time: the time-odd Hilbert–Molien series in closed form

Producer: `analysis/w33_pass11491_time_odd_molien.py`
Certificate: `data/w33_pass11491_time_odd_molien.json`
Regression: `tests/test_w33_pass11486_11492.py`

## Setting

* **Ray invariant of degree k:** a polynomial p(ψ, ψ̄) of bidegree (k, k) that is fixed by the 216 projective Cliffords
  of one qutrit.
* **Time reversal:** complex conjugation composed with any Clifford. It acts on the invariants as an involution, which
  splits them into even ones E_k and odd ones O_k.
* **What an odd invariant means:** f(ψ̄) = −f(ψ). A state with f(ψ) ≠ 0 is not Clifford-equivalent to its own time
  reverse, so every odd invariant is a witness of the substrate's arrow of time on states.
* **What was already known (cite, do not re-claim):** Pass 11419 proved that the first odd invariant has degree 6 and is
  unique there. Pass 11434 wrote it out explicitly as h₆, a signed sum of p_a³p_b²p_c over oriented triples of MUBs.

## Result (exact, proved)

**The generating functions.** With D_k = (1/216)Σ_g |h_k(eig g)|² and Tw_k = (1/216)Σ_g h_k(eig(ḡg)):

| series | closed form |
|---|---|
| time-odd O(t) | **t⁶(1+t) / ((1−t)(1−t²)(1−t³)(1−t⁴)(1−t⁶))** |
| time-even E(t) | (1 − t² + t⁵ + t⁸ − t¹¹ + t¹³) / ((1−t)(1−t²)(1−t³)(1−t⁴)(1−t⁶)) |
| all invariants D(t) = E + O | (1 − t² + t⁵ + t⁶ + t⁷ + t⁸ − t¹¹ + t¹³) / ((1−t)(1−t²)(1−t³)(1−t⁴)(1−t⁶)) |
| trace of time reversal Tw(t) = E − O | **(1+t⁵) / ((1−t)(1−t³)(1−t⁴))** |

**The coefficients:**
* odd: 0,0,0,0,0,0, **1, 2, 3, 5, 8, 11, 16, 22, 29, 38, 49, 61, 77, 95, …**
* even: 1, 1, 1, 2, 3, 4, 6, 8, 11, 15, 19, 24, 32, …
* all: 1, 1, 1, 2, 3, 4, 7, 10, 14, 20, 27, 35, 48, …

Neither the odd sequence nor the total is in the OEIS (searched 2026-10-04).

**Why this is a proof and not a fit.**
1. Every Molien sum is computed exactly. Eigenvalues are identified as 72nd roots of unity, and each sum is reduced
   modulo Φ₇₂ over ℤ; the producer asserts the result is a rational integer divisible by 216.
2. The eigenvalue ratios of g have order dividing 12 (orders 1, 2, 3, 4, 6), and so do the eigenvalues of ḡg (orders
   1, 3, 4). So D_k and Tw_k are quasi-polynomials in k of period 12 and degree ≤ 4.
3. Each closed form is a proper rational function whose poles are 12th roots of unity, so its coefficients are
   quasi-polynomials of the same kind.
4. Two such quasi-polynomials that agree on 60 consecutive values agree everywhere. They are checked on k = 0..72.

## Consequences (each read directly off the closed forms)

1. **Explicit generators, and an independent check of the counts.**
   * The odd Reynolds averages R_m = (1/432) Σ_g sgn(g) m(gψ) of stabiliser-probability monomials m (Pass 11434's
     construction) were evaluated on 14 random rays.
   * They span exactly **1, 2, 3** dimensions in degrees 6, 7, 8 (135, 199 and 280 monomials). The singular values drop
     by about ten decades after the kept ones.
   * This matches O₆, O₇, O₈ from below with explicit polynomials. The Molien series supplies the matching upper bound.
   * **Degree 6:** h₆ = R(p_a³ p_b² p_c), with a, b, c in three different MUBs (Pass 11434).
   * **Degree 7:** h₆·|ψ|² and the new **g₇ = R(p_a⁴ p_b² p_c)**, again over three different MUBs.
     * g₇ is independent of h₆ for **every one of the 24 oriented MUB triples**, with the same residual 0.0241 each
       time, so it is canonical up to sign.
     * The two lowest time arrows of a qutrit are thus the (3,2,1) and (4,2,1) chiralities of an oriented triple of
       MUBs.
   * **Degree 8:** a third direction appears.
     * The products h₆|ψ|⁴ and g₇|ψ|² give only two, because E₁ = E₂ = 1.
     * The new invariant needs **two stabiliser states from the same MUB**, e.g. R(p_a³ p_{a′}² p_b p_c²).
   * **Proved:** the odd invariants are not all multiples of h₆. ℂ[ψ, ψ̄] is a domain, so h₆·E injects into O, and
     O_k − E_{k−6} = 0, 1, 2, 3, 5, 7, 10, 14, … for k = 6, 7, 8, ….
   * **Correction of a first reading (not committed).**
     * The closed form O(t) = t⁶(1+t)/((1−t)(1−t²)(1−t³)(1−t⁴)(1−t⁶)) suggests a free module on generators of degrees
       6 and 7 over primary invariants of degrees 1, 2, 3, 4, 6.
     * That cannot be literal. The only even invariant of degree 2 is |ψ|⁴ (E₂ = 1), which is algebraically dependent
       on the degree-1 invariant |ψ|². So no system of parameters has degrees 1, 2, 3, 4, 6.
     * Indeed a new odd generator appears at degree 8.
     * The rational form is an exact identity of series, not a Hironaka decomposition.
2. **Time reversal is asymptotically balanced.**
   * The leading coefficients (of (1−t)⁻⁵) are 1/72 for O and for E, and 1/36 for D.
   * Tw has only a third-order pole.
   * So O_k/D_k → 1/2: in high degree, half of all state invariants are time-odd.
   * The substrate's time asymmetry is not sparse at high order. It is absent below degree 6 and then **as common as
     symmetry**.
3. **The same reciprocity holds for all four series.** Each satisfies H(1/t) = −t³ H(t).
   * The odd and even numerators (t⁶ + t⁷, and 1 − t² + t⁵ + t⁸ − t¹¹ + t¹³) are both palindromic about 6.5.
   * This is a Stanley-type functional equation. Identifying which general theorem forces the shift 3 (for this
     U(1) × Hessian action on ℂ³ ⊕ ℂ̄³) is left open; the identity itself is exact.
4. **Design check.** D₁ = D₂ = 1 and D₃ = 2 re-derive the known fact that the qutrit Clifford group is a 2-design but
   not a 3-design (classical, external).

**Prior art.**
* Molien's formula and Stanley reciprocity are classical. The Hessian group's own invariant ring (degrees 6, 9, 12, 12,
  Maschke) concerns polynomials in ψ alone, not bidegree (k, k).
* The only other Molien series in the corpus is Pass 11030 (G₁₂ on ℂ², a different action).
* The time-odd / time-even split of the qutrit ray invariants, in closed form, was not found in the corpus or the OEIS.

## Correction (Pass 11515)

* **The error.** The statement that the new degree-8 odd invariant **needs** two stabiliser states from the same MUB is
  wrong: it was an over-read of a largest-residual search.
* **The fact.** The (5, 2, 1) chirality R(p_a⁵ p_b² p_c) over three *different* MUBs is already a new generator in
  degree 8.
* **The full picture.** All five odd generators are the (k, 2, 1) chiralities for k = 3 … 7. The module is free
  through degree 12, with one relation in degree 13 (Pass 11515).
