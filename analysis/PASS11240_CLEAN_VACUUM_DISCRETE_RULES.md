# Pass 11240: the "clean" U(1)′ μ-vacuum does not survive the integer and discrete rules (correction to 11232)

Producer: `analysis/w33_pass11240_clean_vacuum_discrete_rules.py`, which uses Pass 11241's exact lattice engine
Data: Pass 10960 ledger; Pass 10967 space-group and R charges
Certificate: `data/w33_pass11240_clean_vacuum_discrete_rules.json`
Regression: `tests/test_w33_pass11240_clean_vacuum_discrete_rules.py`

**Claim under test.** Pass 11232 found one vacuum in which:
* an unbroken U(1)′ forbids μ;
* the top Yukawa is allowed;
* every other Higgs pair and every vector-like exotic can be massive.

The vacuum has support types {n_2, n_13, n_41, n_55, n_64} in model Z6II_06_1204. The test used rational spans of
U(1) charges.

**What is tested here.**
* **Field choices.** Each support entry is a U(1) *type*. Two of the types contain two fields that differ in their
  discrete charges, so there are **4 choices** of condensing fields, all enumerated.
* **Exact criterion.** For each choice, a coupling is allowed at some order iff its charge vector minus the
  superpotential's R vector lies in the **integer** lattice spanned by the vacuum fields and the discrete
  integrality vectors.
* **Rule sets compared.**
  * gauge U(1)s only, with integer powers, which captures their discrete remnants;
  * plus the space group (Z₆ × Z₃ × Z₂ × Z₂);
  * plus the geometric R-symmetries Z₆^R × Z₃^R × Z₂^R.

## Result

| rule set | μ forbidden for H_u = bl_1, top allowed | light states left (per vacuum choice) | clean |
|---|---|---|---|
| U(1)s, integer lattice | 4/4 | **4 of the 8 v/v̄ pairs** | 0/4 |
| + space group | 4/4 | 4 v/v̄ pairs | 0/4 |
| + space group + R | 4/4 | 2 q/q̄, 5 d/d̄, **all 8** v/v̄, extra Higgs pairs | 0/4 |

1. **μ protection itself is robust.** Adding rules only forbids more, so bl_1·l_j stays forbidden in all 4 choices
   under every rule set. With the R-symmetry the vacuum is moreover **F-flat to all orders** by the lattice
   criterion.
2. **"Clean" was an artifact of rational spans.**
   * A U(1) broken by vacuum fields of charge k leaves a discrete Z_k remnant.
   * The rational test of 11232 ignores these remnants. The integer test keeps them, and they already forbid the
     masses of 4 vector-like v/v̄ pairs. The v are SU(2) doublets with hypercharge ∓1/2.
   * The R-symmetry then also keeps q/q̄, d/d̄ and extra Higgs pairs massless.
   * This is the same pattern Pass 10968 found in Z6II_23: the R-symmetry that makes the vacuum F-flat keeps the
     exotics massless.
3. **Correction recorded.** Pass 11232's "1 clean vacuum (Z6II_06_1204)" is **withdrawn**. Its other statements stand:
   * the "forbidden" verdicts (μ protected in 429 vacua of 54 models);
   * the Z6-I statement for continuous U(1)s.

   Pass 11242 then proves the general reason for the light states.

**Scope.** These are lattice verdicts. "Forbidden" is exact. "Allowed" is a necessary condition, since exponent
positivity is not imposed. Matter parity is not re-tested here (Theorem 9.6 of the paper).
