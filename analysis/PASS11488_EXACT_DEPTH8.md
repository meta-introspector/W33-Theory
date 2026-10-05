# Pass 11488 — the exact reversible fraction at eight cubic gates, by counting operators instead of words

Producer: `analysis/w33_pass11488_exact_depth8.py` (argument: depth, workers)
Certificate: `data/w33_pass11488_exact_depth8.json`
Regression: `tests/test_w33_pass11486_11492.py`

## The idea: deduplicate operators, carry multiplicities

* **The walk.** This is Pass 11436's K-reduced walk: R_{k−1}T ⋯ R₁T·CT, with C over the 216 Cliffords and R over the
  8 coset representatives of K (|K| = 27). It has 216·8^{k−1} words.
* **The observation.** Two words with the same operator (up to phase) have the same extensions. So the walk can run on
  **distinct operators** carrying their word multiplicities.
* **The method.**
  * Deduplicate through depth 6.
  * Then stream the 8 and 64 extensions of each depth-6 operator for depths 7 and 8. No depth-7 set is ever stored.
* **Validation.**
  * The method reproduces the exact P₁ … P₇ of Passes 11312, 11422 and 11436.
  * Deduplication keys are checked at roundings 10⁻⁵, 10⁻⁶ and 10⁻⁹, which give identical class counts.
  * The faster overlap (two batched matmuls) is checked against Pass 11422's routine.
* **Decider.** Theorem 1 (Pass 11355) at the tight cut 3 − 10⁻⁷ (Pass 11459). The gap is reported.

## Results

| k | distinct operators | distinct / 216 | P_k (exact) |
|---|---|---|---|
| 1 | 216 | 1 | 11/12 ✓ |
| 2 | 1 080 | 5 | 19/32 ✓ |
| 3 | 5 616 | 26 | 347/768 ✓ |
| 4 | 29 376 | 136 | 623/2048 ✓ |
| 5 | 166 752 | 772 | 11003/49152 ✓ |
| 6 | 959 040 | 4 440 | 20327/131072 ✓ |
| 7 | (streamed) | — | 353927/3145728 ✓ |
| **8** | (streamed) | — | **665275/8388608 = 0.0793070** |

✓ marks agreement with the earlier full enumerations.

**Reading.**
* **The denominator law holds through k = 8.** 2^{3k−1}·3^{[k odd]} gives 8 388 608 = 2²³ at k = 8.
  * In integers: P_k·216·8^{k−1}/18 = 11, 57, 347, 1869, 11003, 60981, 353927, **1995825**. All are integers, and the
    sequence is not in the OEIS.
  * Neither the integer sequence nor the rationals satisfy a linear recurrence of order ≤ 3.
* **Overlap gap.** Every reversible word reaches 3 − 1.2·10⁻¹⁴. The largest violating overlap through k = 8 is 2.99916, so
  the cut 3 − 10⁻⁷ separates cleanly.
* **Per-gate rates.**
  * One-step: 0.648, 0.761, 0.673, 0.736, 0.693, 0.725, **0.705**. The odd/even alternation persists.
  * Two-step, √(P_{k+2}/P_k): 0.702, 0.716, 0.704, 0.714, 0.709, **0.715**.
  * These rise slowly toward the constant ≈ 0.72 that Pass 11459 sampled at depths 12–30 with the tight cut. They are
    consistent with it, and with ρ₍₈,₈₎ = 0.7228.
  * Exact data still cannot separate those two descriptions; that would need k ≈ 12 or more.
* **Distinct operators per Clifford: 1, 5, 26, 136, 772, 4440** (not in the OEIS, searched 2026-10-04).
  * Ratios 5, 5.2, 5.23, 5.68, 5.75. They grow toward the 8 cosets per gate, less merging.
  * The guessed recurrence a_k = 6a_{k−1} − 4a_{k−2} (growth 3 + √5) fits 1, 5, 26, 136 but predicts 712, not 772.
    **Refuted.**
* **No sqrt law.**
  * Reversible operators are not the sparse fixed points of a random involution on operators.
  * Their distinct count grows ×3.5 per gate at k ≤ 3, against ×5.2 for all operators.
  * An involution model would give about ×√5.2 ≈ ×2.3.

**Prior art.** P₁ … P₇ and the K-reduction are Passes 11312, 11422 and 11436. Multiplicity-weighted deduplication is a
standard transfer technique, new here only as a computational step.
