# Pass 11422 — the exact depth-5 reversible fraction is 11003/49152; the "0.72 per gate" is not a constant rate, and the ρ₍₈,₈₎ coincidence dissolves

Producer: `analysis/w33_pass11422_depth5_exact.py`
Certificate: `data/w33_pass11422_depth5_exact.json`
Regression: `tests/test_w33_pass11418_11422.py`

**Decider (exact).** By Pass 11355, U is reversible iff Uᵀ = λCUC† for some Clifford C, i.e. iff
max_C |tr((CUC†)†Uᵀ)| = 3.
* The 216 Cliffords are checked in one batched contraction.
* The gap is reported, not assumed: every reversible word reaches at least 2.99999999999999, and every violator stays
  at most 2.96531 at all depths 1–5.
* Words are the coset-reduced R_k T ⋯ R_2 T C_1 T of Pass 11312. Their verdicts have the same multiset as the uniform
  walk.
* **Validated:** depths 1–4 reproduce 11/12, 19/32, 347/768 and 623/2048 (2 077 650 violators) exactly.

## Exact values

| cubic gates k | words | reversible fraction | violating fraction |
|---|---|---|---|
| 1 | 216 | 11/12 | 1/12 |
| 2 | 5184 | 19/32 | 13/32 |
| 3 | 124 416 | 347/768 | 421/768 |
| 4 | 2 985 984 | 623/2048 | 1425/2048 |
| **5** | **71 663 616** | **11003/49152 = 0.22386** | **38149/49152** |

* **Denominators:** 12, 32, 768, 2048, 49152 are 2^{3k−1} · 3^{[k odd]} for every k = 1..5.
* **The successive ratios oscillate:** 0.648, 0.761, 0.673, 0.736. That odd/even alternation points to a negative or
  complex subdominant mode in the decay.

## Sampled tail (4 million words per depth, same exact decider)

| k | 6 | 8 | 10 | 12 | 14 | 16 | 20 |
|---|---|---|---|---|---|---|---|
| reversible fraction | 0.15501 | 0.07936 | 0.04075 | 0.02106 | 0.01107 | 0.00586 | 0.00175 |
| standard error | 0.00018 | 0.00014 | 0.00010 | 0.00007 | 0.00005 | 0.00004 | 0.00002 |

Local per-gate rates from two-step ratios:
* 0.7156 (6→8), 0.7166 (8→10), 0.7189 (10→12), 0.7250 (12→14), 0.7274 (14→16);
* 0.740 ± 0.01 (16→20).

**Reading.**
* The "≈ 0.72 per gate" of Pass 11312 is **not a constant decay rate**. The local rate drifts upward with depth, from
  0.716 to above 0.727.
* A pure geometric fit over k = 8–20 gives 0.7218 ± 0.0004. That fit uses a biased model and is not an asymptotic
  estimate; it sits 2.7σ from ρ₍₈,₈₎.
* The **coincidence ρ₍₈,₈₎ = 0.72276 ≈ 0.72** recorded in Pass 11372 therefore **dissolves**: there is no single
  number for it to match.
* Pass 11372 already showed no Fourier rate theorem applies, since the reversible set is Haar-null. The data agree.
* Whether the local rate converges, and to what (0.7228, ρ₍₁₈,₀₎ = 0.781, or 1), is open.

**Prior art.** This extends Pass 11312's exact table by one depth; the exact decider is Pass 11355's Theorem 1.

## Correction (Pass 11459): the sampled tail used too loose a threshold

* **What was wrong.** The sampled values above (k ≥ 6) counted a word as reversible when its best Clifford overlap
  exceeded **2.999**.
  * Exactly reversible words reach 3 to rounding, with a deficit below 10⁻¹⁰.
  * Clifford+T words are dense, so deep violating words come arbitrarily close: their deficits run about 10⁻⁴ to 10⁻³
    by depth 24–30.
  * At the loose cut these near-misses were counted as reversible: 0.7% false at depth 12, **32% at depth 24, 79% at
    depth 30** (Pass 11459, `threshold_sensitivity`).
* **Consequence.** The "upward drift" of the local rate and the conclusion that "0.72 is not a constant rate and the
  ρ₍₈,₈₎ coincidence dissolves" were at least partly this artefact. **Both are withdrawn.** Pass 11459 re-measures
  with the tight cut 3 − 10⁻⁷.
* **Not affected.**
  * The exact values for k ≤ 5: violators stayed below 2.96531.
  * The exact values for k ≤ 7 (Pass 11436): violators stayed below 2.99861.
