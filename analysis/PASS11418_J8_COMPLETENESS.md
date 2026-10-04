# Pass 11418 — J₈ survives every test that breaks J₆; a ratio bound is impossible for any witness near the identity

Producer: `analysis/w33_pass11418_j8_completeness.py` (main run; `--asym`)
Certificate: `data/w33_pass11418_j8_completeness.json`
Regression: `tests/test_w33_pass11418_11422.py`

**Question.** Pass 11369 conjectured t\* = 4: the degree-8 witness J₈ is complete on PU(3), so U is reversible iff
J₈(U) = 0. Every test below uses J₆, which is known to fail, as the positive control.

## (c) The pseudo-reflection locus, where J₆'s second component lives (Pass 11419)

| | J₆ | J₈ |
|---|---|---|
| rank-1 in (ray, β): second singular value | 6×10⁻¹⁴ (**factorises**) | 0.05 (does not) |
| zeros found by 150 minimisations | 137 | 36 |
| zeros on non-real-type (non-reversible) rays | **94** | **0** (every zero has Clifford overlap ≥ 0.99999999998) |

**J₈ has no spurious zero on the pseudo-reflection locus.**

## (b) Descents started at J₆'s spurious zeros

* 2000 Haar seeds gave 307 certified J₆-spurious points, 80 of them near-degenerate.
* J₈ descent from each reached 303 zeros. **None is certified non-reversible.** The largest distance at a zero is
  0.0039, with the next witness at rounding level: convergence onto the reversible set.

## (a) The ratio test, and why it must fail for every witness

* min J/dist² over {dist ≥ δ}, 200 Nelder–Mead starts each:
  * J₆ reaches −6×10⁻¹⁴ at δ = 0.1 and 0.3: genuine zeros far from the reversible set.
  * J₈ has **minimum value 7.2×10⁻⁷ at δ = 0.3**, attained on the boundary dist = 0.300, and 3.6×10⁻¹² at δ = 0.1.
* **The small δ = 0.1 value is explained, not a counterexample.**
  * On pseudo-reflections, J₈(I + λP) = |λ|¹²(252Δ₆ − 24|λ|²Δ₇ + |λ|⁴Δ₈), where Δ_k are the moment gaps of Pass
    11419.
  * Checked: J₈/(|λ|¹²Δ₆) = 177, 222, 238 → 252 as β = 1.0, 0.6, 0.4 → 0, identically on three rays.
  * The distance to the reversible set is O(|λ|). So **every** witness vanishes like dist¹² approaching the identity
    along these families, and no bound J ≥ c·dist² can hold. The fit 0.36·dist¹¹ reproduces both values
    (6.4×10⁻⁷ at 0.3, 3.6×10⁻¹² at 0.1).
  * The 0.1 value sits at double-precision rounding and is not certified positive.
* The meaningful quantity is min J₈ over {dist ≥ δ}. It is positive on the sample at δ = 0.3 and consistent with a
  high-order vanishing law below that.

**Reading.**
* J₈ passes all three tests that J₆ fails: pseudo-reflections, seeded descents, and positivity away from the
  reversible set. Together with Pass 11369's 1000 clean global descents, t\* = 4 remains the conjecture.
* The new structural fact: near the identity every T-witness is extremely weak, of order dist¹². T-violation of a
  nearly trivial tick is a twelfth-order effect in its distance from the reversible set.

**Open.**
* A proof that J₈ is complete.
* A certified high-precision lower bound for min J₈ on {dist ≥ δ}.
