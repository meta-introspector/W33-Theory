# Pass 11171 — one qudit's Bell value in time never reaches 4, but climbs to within ~3.4/d of it (d up to 24); d = 6, 7 corrected

Producer: `analysis/w33_pass11171_temporal_cglmp_large_d.py`
Scan: `analysis/w33_pass11171_scan_temporal_cglmp_large_d.py`, frozen as `data/w33_pass11171_temporal_cglmp_large_d.json`
Regression: `tests/test_w33_pass11171_temporal_cglmp_large_d.py`

**Theorem: 4 is never attained.** The scenario is one qudit with rank-1 Lüders measurements, A_x then B_y, so
P(a,b|x,y) = ⟨a_x|ρ|a_x⟩ |⟨b_y|a_x⟩|².
1. I_d = 4 forces all four terms to their maximum. So P is supported on b = a + s_xy, with s₀₀ = s₀₁ = s₁₁ = 0 and
   s₁₀ = 1.
2. For a ∈ S₀ = supp p(·|0), the vector |a₀⟩ equals |b0_a⟩ and also |b1_a⟩. For a ∈ S₁, |a₁⟩ equals |b0_{a+1}⟩ and also
   |b1_a⟩.
3. ⟨a_x|ρ|a_x⟩ = 0 implies ρ|a_x⟩ = 0. So supp ρ ⊂ span{b1_a : a ∈ S₀ ∩ S₁}, and S₀ ∩ S₁ is nonempty.
4. For a in S₀ ∩ S₁: b0_a = b1_a = b0_{a+1}, so two vectors of one orthonormal basis coincide. Contradiction.

So **I_d < 4 for every finite d**. A classical system with memory reaches 4 (Pass 11144). The qudit's only memory
between the two times is the post-measurement vector |a_x⟩, which cannot encode both a and x perfectly.

**Numerics.**
* For fixed bases, the optimal state is the top eigenvector of W = Σ_x A_x diag(w_x) A_x†, so I_d = max over the four
  bases of λ_max(W).
* The optimiser is Riemannian gradient ascent on U(d)⁴. The analytic gradient is checked against finite differences to
  10⁻⁹. There are 24 restarts per d.
* All values are **lower bounds**.

| d | I_d (temporal) | d·(4 − I_d) |
|---|---|---|
| 6 | **3.4552** (was 3.4531) | 3.27 |
| 7 | **3.5209** (was 3.5168) | 3.35 |
| 8 | 3.5762 | 3.39 |
| 10 | 3.6576 | 3.42 |
| 12 | 3.7147 | 3.42 |
| 14 | 3.7562 | 3.41 |
| 16 | 3.7873 | 3.40 |
| 20 | 3.8310 | 3.38 |
| 24 | 3.8528 | 3.53 (probably under-optimised) |

**Corrections.**
* **Pass 11155:** the d = 6, 7 optima (3.4531, 3.5168) were local optima; the improved values are 3.4552 and 3.5209.
* **Passes 11155 and 11160:** the deficit exponents d^−0.65 and d^−0.75, fitted to d ≤ 10, are superseded. For d ≥ 6
  a power law to 4 fits with exponent 0.97, and a free-limit fit gives 4.03. **The deficit falls like ≈ 3.4/d.**
* The d = 4, 5, 8, 9, 10 optima are reproduced to within 2·10⁻⁴.

**Reading.** The data support sup_d I_d = 4, approached like 1/d and never attained. Time lets one qudit come
arbitrarily close to the algebraic Bell maximum, but finite memory keeps it strictly below.

**Prior art.** Budroni–Emary (arXiv:1309.3678) showed that temporal Leggett–Garg values grow with dimension. Fritz
(2010) showed temporal and spatial CHSH coincide for qubits. The finite-d unattainability proof and the CGLMP numbers
to d = 24 are ours, as far as we found.
