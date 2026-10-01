# Pass 11265 — time-reversal fidelity is a Fourier maximum; three-qutrit numerical branch and bound; residual quantisation audit

Producer: `analysis/w33_pass11265_exact_fidelity.py` (uses the Weyl machinery of Pass 11252)
Certificate: `data/w33_pass11265_exact_fidelity.json`
Regression: `tests/test_w33_pass11265_11269.py`

## The formula

Write U = 3⁻ⁿ Σ_p c(p) W(p) in symmetric Weyl operators (Pass 11252). For an anti-symplectic L, the best reversal in the
Weyl coset of V_L has fidelity

> **F(L) = max_b | Σ_p conj(c(p))·c(Lp)·ω^{⟨b,Lp⟩} | / 9ⁿ**,  and F_T(U) = max_L F(L).

* The inner maximum is a Fourier transform over F₃²ⁿ, taken over all b at once.
* Validated against the matrix twirl on 40 random cases: maximum deviation 1.7×10⁻¹⁵.

## Exhaustive numerical maximisation by best-first branch and bound over L

* Build L one basis vector at a time.
* For a partial map on the span S, bound every completion by:
  * the span part, maximised exactly over the 3^k characters of L(S); plus
  * the rest, by the rearrangement inequality on |c| (L maps the complement of S bijectively onto the complement of
    L(S)).
* In exact arithmetic the largest open bound is an upper bound. The implementation uses complex floating arithmetic
  and prunes at 10⁻¹², so an exhausted tree is a numerical certificate at that tolerance, not an interval-arithmetic proof.

**Two qutrits.** 14/14 words (2 named and 12 random, with 1–4 cubic gates) agree with the exhaustive 51 840 × 81
search, in 6 to 19 521 nodes.

**Three qutrits.**

| tick | F_T | status |
|---|---|---|
| [(T⊗T)·SUM] ⊗ I | **F_min = 0.8440296…** | tree exhausted at 10⁻¹² (811 nodes) |
| (T⊗T⊗T)·SUM₁₂ | **F_min** | tree exhausted at 10⁻¹² (830 nodes) |
| (T⊗T⊗T)·SUM₂₃·SUM₁₂ | **∈ [0.237, 0.963]** | bracket only; the bound is weak because \|c\| is nearly flat |

**What this establishes numerically.** At 10⁻¹² pruning tolerance, a third qutrit does not improve reversal of the
embedded minimal violator. Pass 11239's product cases have exhausted search trees. For the genuinely three-qutrit violator, F_T stays open; a sharper bound needs
the phases, not only the magnitudes.

## Is the optimal residual quantised? It appears only at the bottom of the audited sample

* For every level of the Pass 11235 two-qutrit stream, the first 600 optimal reversals at most were examined.
* Only **F_min, F_min² and 2/3** have an optimal residual V U* V† U whose eigenvalues lie in μ₉. For F_min² it is the
  sumset {−1,0,1} + {−1,0,1} of Pass 11251.
* For **every other level** (0.939, 0.810, 0.764, …), none of the examined optimal residuals has eigenvalues in μ₈₁.
* Thus the audit gives no evidence that "quanta multiply" beyond the three simplest levels; unexamined maximisers are
  not excluded.
