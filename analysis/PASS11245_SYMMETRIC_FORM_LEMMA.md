# Pass 11245: the counting lemma behind the Hesselink identity

Producer: `analysis/w33_pass11245_symmetric_form_lemma.py`
Certificate: `data/w33_pass11245_symmetric_form_lemma.json`
Regression: `tests/test_w33_pass11245_symmetric_form_lemma.py`

**Background.** Pass 11244's identity 2 says that among the unipotents of Sp(2m,2) of Jordan type λ, the fraction with
χ₂ = 0 is 2^{−m₂} for m₂ even and 0 for m₂ odd. Here χ₂ = 0 means the symmetric form ω(v,(u−1)w) induced on the
m₂-dimensional space of size-2 blocks is alternating.

**Proved (closed forms, k ≤ 12; brute force for k ≤ 5).** Over F₂,

#{nondegenerate alternating k×k} / #{nondegenerate symmetric k×k} = 2^{−k} for k even, and 0 for k odd.

The alternating count is |GL_k(2)|/|Sp_k(2)|; the symmetric count is MacWilliams' formula. So identity 2 is exactly
the statement that **the induced form is uniformly distributed over the nondegenerate symmetric forms**.

**Status.**
* Identity 2 is reduced to that uniformity statement, which is not proved here.
* Identity 1 (Jordan-type weights equal the odd-q formula at q = 2) remains verified only, on 92 of 92 Jordan types.
* The qubit limit 0.66516065 of Pass 11244 stays conditional.
