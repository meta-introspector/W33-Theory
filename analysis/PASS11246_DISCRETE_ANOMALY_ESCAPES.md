# Pass 11246: no discrete escape — μ protection leaves a colored state massless in all 16 326 sampled cases

Producer: `analysis/w33_pass11246_discrete_anomaly_escapes.py`. It combines Pass 11241's integer-lattice engine with
Pass 11242's colored perfect-matching test.
Certificate: `data/w33_pass11246_discrete_anomaly_escapes.json`

**Question.** Pass 11242 proved that a continuous anomaly-free U(1)′ forbidding μ leaves a colored state massless.
For discrete Z_N and R-symmetries the anomaly argument holds only as a congruence, and Green–Schwarz shifts can enter.
That is how Z₄^R constructions evade it in the literature. So we ask whether these models' **full** vacuum symmetries
ever forbid μ while letting every colored component be massive.

**Test.**
* Rules: U(1)s with their integer remnants, the space group and, for Z6-II, the geometric R-symmetries.
* Vacua: up to 60 FI-ray supports × 2 field choices per model, in all 128 D-flat models.
* For each μ-protected H_u, test every H_d for a perfect colored matching built from lattice-allowed mass terms:
  vector-like, H_u-Yukawa and H_d-Yukawa.

**Result.**
* **16 326 μ-protected cases, 0 escapes.** In every case some colored state stays massless, whatever H_d is chosen.
* Z6-I has no protected cases at all, consistent with Pass 11241.

**Reading.** The light-colored-state obstruction proved for continuous U(1)′ holds empirically for the discrete and
R-symmetries these orbifolds actually carry. A Green–Schwarz-type evasion does not occur in the sampled vacua.
Together with Pass 11241, μ cannot come from any symmetry these models' vacua carry without a light colored exotic or
a massless quark. **μ stays OPEN.** The remaining routes are:
* R-symmetries in Z6-I, whose charges are absent from the data (Pass 11247);
* non-symmetry mechanisms, such as accidental suppression.

**Scope.** The vacua are sampled and the "allowed" verdicts are necessary conditions. So "0 escapes" is a strong
negative, not a proof.
