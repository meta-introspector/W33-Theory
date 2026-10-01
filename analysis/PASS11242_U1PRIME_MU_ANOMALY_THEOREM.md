# Pass 11242: a U(1)′ that forbids μ always leaves a colored state massless (theorem, checked on 812 cases)

Producer: `analysis/w33_pass11242_u1prime_mu_anomaly_theorem.py`
Certificate: `data/w33_pass11242_u1prime_mu_anomaly_theorem.json`
Regression: `tests/test_w33_pass11242_u1prime_mu_anomaly_theorem.py`

**Theorem.** Let U(1)′ be free of the SU(3)²–U(1)′ anomaly, and suppose t(H_u) + t(H_d) ≠ 0. Then the colored fields
cannot all become massive through U(1)′-allowed mass terms: vector-like masses, the up-type Yukawas q·ū·H_u and
q̄·u·H_d, and the down-type Yukawas q·d̄·H_d and q̄·d·H_u.

**Proof.**
* A complete set of colored masses is a perfect matching of triplet and antitriplet components in each
  electric-charge sector. Every matched pair carries charge −t(H_u), −t(H_d) or 0.
* Count the pair types. Up sector: n₁ pairs (q,ū) using H_u, n₄ pairs (u,q̄) using H_d. Down sector: m₁ pairs (q,d̄)
  using H_d, m₄ pairs (d,q̄) using H_u. Counting q against q̄ gives n₁ − n₄ = 3 and m₁ − m₄ = 3.
* So the numbers of H_u-pairs and H_d-pairs are equal: a = n₁ + m₄ = m₁ + n₄ = b, and a ≥ 3.
* The anomaly equals the sum of the matched pairs' charges, −a(t(H_u) + t(H_d)). It vanishes, so
  t(H_u) + t(H_d) = 0, a contradiction. ∎

No generation universality is assumed. The family-universal special case, that a U(1)′ forbidding μ needs exotics,
is known from U(1)′-extended MSSM model building.

**Checks.**
* **Anomaly input.** The SU(3)²–U(1) anomaly vanishes for every non-anomalous U(1) of all 215 models.
* **Exhaustive matching test.** Over every FI-cancelling vacuum span of Pass 11232 and every μ-protected H_u: **812
  cases, 0 violations.** No choice of H_d ever gives a perfect colored matching. The matching is done at component
  level, respecting hidden-sector multiplicities.
* **Where the failure sits:** up sector only in 3 cases, down sector only in 46, both in 763.

**Consequence.** In these heterotic models every U(1) except the anomalous one is anomaly-free, and the anomalous one
is broken by the FI vacuum. So **no continuous vacuum symmetry can solve the μ problem without leaving a colored state
massless.** This is the general reason behind Pass 11232's pattern and Pass 11240's correction. The discrete
analogue (Z_N and R-symmetries, where Green–Schwarz shifts can intervene) is Pass 11246.
