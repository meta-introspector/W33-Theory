# Pass 11435 — J₈ is positive on both components of J₆'s blind spot; an exact positivity law on the pseudo-reflection component (and a correction)

Producer: `analysis/w33_pass11435_j8_positivity.py`
Certificate: `data/w33_pass11435_j8_positivity.json`
Regression: `tests/test_w33_pass11433_11437.py`

## Correction of Passes 11418 and 11419 (a displayed formula)

* **What was wrong.** For U = 1 + λP, P = |ψ⟩⟨ψ|, every trace is real:
  tr(U†CUC†) = 3 − |λ|²(1 − x) with x = |⟨ψ|Cψ⟩|², and the same with y for ψ̄. So the expansion coefficients are
  E_k = avg(1 − x)^k − avg(1 − y)^k, and **not** the moment gaps Δ_k = avg x^k − avg y^k.
* **The correct formula.** Since Δ₁ = ⋯ = Δ₅ = 0,

> **J₈(1 + λP) = |λ|¹² [252 Δ₆ − 24|λ|²(7Δ₆ − Δ₇) + |λ|⁴(28Δ₆ − 8Δ₇ + Δ₈)]**

* **What changes.** Passes 11418 and 11419 displayed |λ|¹²(252Δ₆ − 24|λ|²Δ₇ + |λ|⁴Δ₈). The leading term, the limit
  J₈/(|λ|¹²Δ₆) → 252, and every conclusion drawn from them are unchanged. Both notes and the 11418 producer now carry
  the corrected formula.
* **How it was caught.** The old formula would have made J₈ negative near the identity on {h₆ = 0}, where Δ₇ > 0.
  That is impossible, since J₈ is a sum of squares.

## Theorem (exact, given two measured signs): J₈ > 0 on the pseudo-reflection component

* On J₆'s component {Δ₆ = 0} (Pass 11419):

> J₈ = |λ|¹⁴ [24Δ₇ + |λ|²(Δ₈ − 8Δ₇)],  |λ|² ∈ (0, 4].

* This is linear in |λ|², so it is positive for every β ≠ 0 **iff Δ₇ > 0 and Δ₈ > 2Δ₇**.
* **Measured** on 195 non-real-type rays of {h₆ = 0}:
  * Δ₇ ≥ 9.5×10⁻¹¹ > 0;
  * Δ₈/Δ₇ ≥ 3.76 > 2;
  * |Δ₆| ≤ 1.8×10⁻¹⁵;
  * the formula agrees with direct evaluation to 1.1×10⁻⁶ at β = 2.0 and 2.9;
  * the smallest J₈ at β = π is 3.0×10⁻⁵.
* Δ₇ becomes small only where the hypersurface approaches real-type rays, i.e. the reversible set.

## The generic component

* At 24 certified generic-spectrum J₆-spurious points, J₈ ranges over 1.3×10⁻⁴ … (median 0.11). It is positive
  everywhere.

## High precision at small distance

* Low-J₈ points were re-found at distance ≥ δ (Nelder–Mead, 40 starts) and evaluated with mpmath at 50 digits:

| δ | distance | J₈ (50 digits) |
|---|---|---|
| 0.1 | 0.104 | **+4.4976×10⁻⁹** |
| 0.2 | 0.200 | **+9.7959×10⁻⁶** |

* Both are certified positive at those points; double precision agrees to 3 or 4 digits.
* With 200 starts, Pass 11418 reached 3.6×10⁻¹² at δ = 0.1, consistent with the dist¹² decay near the identity. That
  smaller value was not re-evaluated here.

**Reading.**
* Both components of J₆'s blind spot are seen by J₈:
  * the pseudo-reflection component through an exact law with two positive moment gaps;
  * the generic component through direct values.
* Away from the identity J₈ is bounded below on every sample.
* A global certified bound on {dist ≥ δ} (8 dimensions, with a Lipschitz covering) was not attempted. t\* = 4 remains a
  conjecture with no counterexample.
