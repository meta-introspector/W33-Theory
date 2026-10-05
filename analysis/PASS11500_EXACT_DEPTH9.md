# Pass 11500 — the exact reversible fraction at nine cubic gates, and a cubed-phase test that decides every word

Producer: `analysis/w33_pass11500_exact_depth9.py` (arguments: depth, workers)
Certificate: `data/w33_pass11500_exact_depth9.json`
Regression: `tests/test_w33_pass11498_11505.py`

## Method

**The walk.** This is Pass 11488's weighted operator deduplication.
* Deduplication runs through depth 7: 5 567 184 distinct operators carrying 216·8⁶ word multiplicities.
* Depths 8 and 9 are streamed as the 8 and 64 extensions of each depth-7 operator, about 4·10⁸ words.

**One new ingredient: a cubed-phase prefilter.**
* Pass 11252's exact criterion says U is reversible ⟹ c(Lp) = μ·ω^{⟨b,Lp⟩}·c(p) for an anti-symplectic L of F₃²
  (24 of them), with c(p) = tr(W(p)†U)/3.
* Cubing removes ω^{⟨b,Lp⟩}. So a reversible U must have |c(Lp)| = |c(p)| and (c(Lp)/c(p))³ constant on the support.
* Words that fail this are violating.
* Words that pass are confirmed by Theorem 1 (overlap = 3, cut 3 − 10⁻⁷).
* This makes the stream about 3× faster.

**Checks.**
* A control on 20 000 random depth-8 words agrees with the full Theorem-1 decider.
* Deduplication keys at 10⁻⁵ and 10⁻⁶ agree at every depth.
* At depth 7 the 10⁻⁹ key *splits* 29 operators (floating noise). That is harmless, because split copies keep their own
  weights. Merging would be the danger, and the coarser keys exclude it.
* P₁ … P₈ are reproduced exactly.

## Results

| k | 7 | 8 | **9** |
|---|---|---|---|
| P_k | 353927/3145728 ✓ | 665275/8388608 ✓ | **11519231/201326592 = 0.0572166** |

* **The denominator law 2^{3k−1}·3^{[k odd]} holds for k = 1 … 9.** 201 326 592 = 2²⁶·3.
* **In integers,** P_k·216·8^{k−1}/18 = 11, 57, 347, 1869, 11003, 60981, 353927, 1995825, **11519231**.
* **Per-gate rates.**
  * One-step: …, 0.6928, 0.7255, 0.7049, **0.7215**.
  * Two-step √(P_{k+2}/P_k):
    * from odd k: 0.7021, 0.7039, 0.7089, **0.7131** (rising);
    * from even k: 0.7158, 0.7140, 0.7151.
  * Both chains approach the constant ≈ 0.72 that Pass 11459 sampled at depths 12–30.
  * The exact data are consistent with ρ₍₈,₈₎ = 0.7228 but do not yet reach it.
  * **Extrapolation fails.**
    * A two-mode fit P_k ≈ Aρ^k + B(−σ)^k gives ρ = 0.7119, 0.7126, 0.7134 as the fit starts from k = 3, 4, 5. The
      residuals are 5, 3 and 2 × 10⁻³, so the model is inadequate and ρ is drifting upward.
    * Aitken acceleration of the two two-step chains gives 0.701, 0.733 and 0.715, which is unstable.
    * Nine exact terms do not determine the asymptotic rate. No value is claimed.
* **The cubed-phase test is an exact decider on every word computed.**
  * Across all words to depth 9, **no violator passed the prefilter**. The largest overlap among prefilter-passing
    violators is undefined, because there are none.
  * So on one-qutrit Clifford+T words, "|c(Lp)| = |c(p)| and (c(Lp)/c(p))³ constant for some anti-symplectic L"
    **coincides** with reversibility.
  * Cubing discards the requirement that the residual ω-phase be *affine*, so the equivalence is not automatic. For
    these words the affine condition is never the deciding one.
  * Whether the coincidence holds for all of PU(3) is **open**. On the measure-zero set where the moduli match, a
    non-affine phase could in principle survive.

**Prior art.** P₁ … P₈ come from Passes 11312, 11422, 11436 and 11488. The criterion is Pass 11252's, and the
reversibility theorem is Pass 11355's.
