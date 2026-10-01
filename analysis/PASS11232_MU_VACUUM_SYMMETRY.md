# Pass 11232: can a vacuum symmetry forbid the μ term? An exact U(1) test on all 215 models

> **Correction (Pass 11240).** The single "clean" vacuum (Z6II_06_1204) is **withdrawn**. The clean test used rational
> spans. Integer powers of the vacuum fields leave discrete remnants of the broken U(1)s, and these already keep 4
> v/v̄ pairs massless; the space-group and R rules then keep q/q̄, d/d̄ and extra Higgs pairs massless too. The
> "forbidden" verdicts below (μ protected in 429 vacua of 54 Z6-II models; none in Z6-I) stand. Pass 11242 proves
> why protection costs light colored states.

Producer: `analysis/w33_pass11232_mu_vacuum_symmetry.py`
Input: Pass 10960's exact-rational ledger, `data/w33_pass10960_heterotic_left_chiral_ledger.json.gz`
(sha256 dd7a2c8b…; its D-flat machinery is reused).
Certificate: `data/w33_pass11232_mu_vacuum_symmetry.json`
Regression: `tests/test_w33_pass11232_mu_vacuum_symmetry.py`

**The open item.** The paper's scorecard (sec09) says of the μ term: "Arises at order 3 or 4 in every model; none is
μ-free through order five. Must come from a vacuum symmetry." This pass tests that sentence's continuous part
exactly.

**Setup.**
* **Vacuum.** An FI-cancelling D-flat direction of the singlets that are neutral under the SM and the hidden sector.
  Every such direction is a nonnegative combination of extreme rays of the D-flat cone, and its support contains a
  ray's support. Rays therefore leave the largest unbroken symmetry, so they suffice for an existence question.
* **Rays.** All extreme rays are computed exactly (cdd, gmp) in 125 of the 128 D-flat models. In the three models
  with more than 44 singlet types (Z6II_06_744, Z6II_52_213, Z6II_12_2184), 3000 LP-sampled vertices are used, each
  re-solved exactly on its support.
* **Unbroken U(1)s.** A vacuum S leaves unbroken the directions t with t·q_s = 0 for all s ∈ S. A coupling of total
  charge Q is forbidden to all orders by an unbroken U(1) iff Q ∉ span_ℚ{q_s}.
* **μ-protected up-type Higgs.** bl_i with every μ_ij = bl_i·l_j forbidden, and some top Yukawa q_a·bl_i·bu_b with
  charge in the span.
* **Clean.** In addition:
  * every other bl_k gets a mass with some l, by structural (matching) rank, so exactly one Higgs pair is light;
  * every vector-like pair (q/bq, u/bu, d/bd, e/be, v/bv, x/bx) gets a mass;
  * there is a singlet s′ outside the vacuum with s′·bl_i·l_j allowed. This is the "U(1)′ solution": μ = λ⟨s′⟩
    once s′ condenses and breaks U(1)′.

## Results

| | count |
|---|---|
| models / anomalous / with FI-cancelling D-flat directions | 215 / 211 / 128 (reproduces 23 Z6-I + 105 Z6-II) |
| distinct vacuum spans examined | 2667 |
| vacua with some μ coupling forbidden by a vacuum U(1) | in 127 of 128 models |
| vacua with a **μ-protected H_u and an allowed top Yukawa** | 429 vacua, **54 models**, all **Z6-II**; 53 have a μ-generating s′ |
| **Z6-I** (the class containing the W33 vacuum's 42 Standard Models) | **0 of 23** D-flat models have any μ-protected vacuum (all rays, exact) |
| **clean** (one light Higgs pair, every vector-like exotic massive, s′ exists) | **1 vacuum, 1 model: Z6II_06_1204** |

1. **In the W33 vacuum's own class (Z6-I), μ cannot be forbidden by any continuous vacuum symmetry.** This holds in
   all 23 D-flat models, over all extreme rays, exactly. If μ is to come from a vacuum symmetry there, it must be a
   discrete one: a space-group or R-symmetry, as in Kappl et al. arXiv:0812.2120 and Pass 10968.
2. **In Z6-II the protecting U(1)′ almost always keeps exotics light.** In 428 of the 429 protected vacua, either
   another Higgs pair or a vector-like exotic stays massless. Among the stored examples (up to 3 per model), the light
   states are d/bd pairs (108), extra Higgs pairs (94), v/bv pairs (71), x/bx (28) and q/bq (6). This is the
   continuous counterpart of Pass 10968's finding: the R-symmetry that makes the vacuum F-flat also keeps 3 + 5 + 3
   exotics massless. The symmetry that protects μ tends to protect everything vector-like.
3. **One clean vacuum exists: Z6II_06_1204.**
   * Its support is {n_2, n_13, n_41, n_55, n_64}. Of the model's seven U(1)s, five are broken. What remains is
     hypercharge and one U(1)′, whose generator is given exactly in the certificate.
   * Normalised by 17/4, the U(1)′ charges are:
     * H_u = bl_1: **−7**;
     * the doublets l_j: −3 (six of them), +2 (two), −8 (two). So every μ_1j has charge −10, −5 or −15, and is
       forbidden;
     * the top Yukawa q_1·bl_1·bu: **6 − 7 + 1 = 0**, allowed;
     * every other bl_k pairs off;
     * all 2 bq/q, 8 d/bd and 8 v/bv vector-like pairs are massive;
     * 11 singlets can play s′.
   * So, at the level of gauge U(1) selection rules, this model has a vacuum in which μ is forbidden to all orders by
     an unbroken U(1)′ and arises only when U(1)′ breaks. That is the UMSSM mechanism (Cvetič–Langacker and
     successors, standard; not re-derived here).

## Limits (read before quoting)

* Only the continuous gauge U(1)s are tested. Discrete space-group and R rules can forbid more, never less.
  * The **forbidding** verdicts (μ protected; no protection in Z6-I being *continuous*) are rigorous.
  * The **allowing** verdicts (top Yukawa, exotic masses, the s′ coupling) are necessary conditions only.
  * The clean vacuum must still be checked against the space-group rules, which need the orbifolder dump of model
    Z6II_06_1204 (not in this repo).
* F-flatness is not imposed (the scope of Pass 10960/10962). Pass 10960's Theorem 9.6 still applies: no
  supersymmetric singlet vacuum of these models preserves matter parity. So the clean vacuum does not fix proton
  stability.
* The unbroken U(1)′ is a massless Z′ until s′ condenses. The mechanism therefore needs ⟨s′⟩ near the TeV scale, and
  LHC Z′ bounds then constrain it. That is a phenomenological cost, not a derivation of the scale.
* Three models use sampled vertices, so "not found" there is not proof (Z6II_06_744 has protected vacua anyway).

**Scorecard reading.** The μ row stays **OPEN**, now sharpened:
* in the W33 vacuum's own class a continuous vacuum symmetry cannot do it, so μ needs a discrete symmetry there;
* in the neighbouring Z6-II class a U(1)′ can do it in 54 models, and cleanly, at the gauge level, in exactly one.
