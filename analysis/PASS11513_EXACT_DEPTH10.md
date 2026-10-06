# Pass 11513 — the exact reversible fraction at ten cubic gates

Producer: `analysis/w33_pass11513_exact_depth10.py` (argument: workers)
Certificate: `data/w33_pass11513_exact_depth10.json`
Regression: `tests/test_w33_pass11511_11515.py`

## Method

This is Pass 11500's method one level deeper.

**Depth 8 by hashes.**
* Depth 8 has 44 537 472 candidates: the 8 extensions of the 5 567 184 distinct depth-7 operators.
* Each candidate is stored as (depth-7 index, extension) with two 64-bit hashes of its phase-fixed entries, rounded at
  10⁻⁶ and at 10⁻⁵.
* Deduplication is np.unique on the 10⁻⁶ hash.
* **Merge check.** The independent 10⁻⁵ hash must give the same number of classes. A hash collision, or two distinct
  operators agreeing to 10⁻⁶, would show up as a mismatch.
* **Result: 32 445 792 distinct depth-8 operators** (hash check equal).

**Stream.**
* Each distinct depth-8 operator is rebuilt from (index, extension).
* Its own reversibility gives P₈ (validation), its 8 extensions give P₉ (validation), and its 64 extensions give P₁₀.
* **Faster prefilter.** Pass 11500's exact cubed-phase prefilter now computes phases only for the (word, L) pairs whose
  moduli already match. Its decisions are identical on 20 000 depth-9 words, and it is 4.3× faster. Theorem 1 confirms
  each survivor.

## Results

| k | 8 | 9 | **10** |
|---|---|---|---|
| P_k | 665275/8388608 ✓ | 11519231/201326592 ✓ | **21873707/536870912 = 0.0407430** |

* **The denominator law 2^{3k−1}·3^{[k odd]} holds for every k = 1 … 10.** 536 870 912 = 2²⁹.
* **In integers,** P_k·216·8^{k−1}/18 = 11, 57, 347, 1869, 11003, 60981, 353927, 1995825, 11519231, **65621121**.
* **Two-step rates √(P_{k+2}/P_k).**
  * From odd k: 0.7021, 0.7039, 0.7089, **0.7131**.
  * From even k: 0.7158, 0.7140, 0.7151, **0.7168**.
  * Both chains are now rising monotonically at their ends. They are still below the ≈ 0.72 sampled at depths 12–30
    (Pass 11459) and below ρ₍₈,₈₎ = 0.7228.
* **Overlap gap.** Reversible words reach 3 − 1.4·10⁻¹⁴. **No violator passes the prefilter at any depth ≤ 10**, so the
  cubed-phase test is exact on every one-qutrit Clifford+T word through depth 10.
  * This contrasts with PU(3), where Pass 11512 found non-reversible unitaries that pass it.
* **Cost.** 22 882 s on 5 workers.

**Prior art.** P₁ … P₉ come from Passes 11312, 11422, 11436, 11488 and 11500.
