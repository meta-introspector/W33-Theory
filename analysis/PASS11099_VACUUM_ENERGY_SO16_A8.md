# Pass 11099 — the one-loop vacuum energy of the SO(16)×SO(16) A8 models is positive and runs away

Producers: `analysis/w33_pass11099_lambda10_so16xso16.py` (the 10D modular integral),
`analysis/w33_pass11099_vacuum_energy_so16_a8.py`
Frozen evidence: `data/w33_pass11099_massless_bose_fermi.json`
Certificates: `data/w33_pass11099_lambda10_so16xso16.json`, `data/w33_pass11099_vacuum_energy_so16_a8.json`
Regression: `tests/test_w33_pass11099_vacuum_energy_so16_a8.py`

The 104 tachyon-free three-generation models of Pass 11095 are non-supersymmetric, so their one-loop vacuum energy Λ
does not cancel. This pass computes what controls it.

## A. The ten-dimensional SO(16)×SO(16) vacuum energy, recomputed

The partition function is written in SO(2n) level-1 characters (Alvarez-Gaumé, Ginsparg, Moore, Vafa 1986) and
integrated over the fundamental domain. For τ₂ > 1 exact q-expansions are used, where only level-matched terms
survive. The thin region below τ₂ = 1 is done by Gauss–Legendre quadrature.

| check | result |
|---|---|
| modular invariance of the integrand (τ → −1/τ, τ → τ+1) | relative error ≤ 10⁻¹⁴ |
| level-matched massless term | **−2112 = n_B − n_F** = 8·(8 + 240) − (8·256 + 8·256) |
| no level-matched tachyon | ✓ |
| I = ∫_F d²τ/τ₂² Z | **−725.96717635**, stable to 12 digits under refinement |

So **Λ₁₀ = −½ M¹⁰ I = +363 M¹⁰ > 0** (M = M_s/2π). It is **positive**, as Alvarez-Gaumé et al. found.

## B. What this means for the A8 models (large volume)

The orbifold group is G = Z₂W × Z3 with |G| = 6, and Z = (1/6) Σ_{g,h∈G} Z(g,h). Only the 4 of 36 sectors with
g, h ∈ {1, β} carry T⁶ zero modes (β, the Witten element, acts trivially on space), and together they give

  (1/6) Σ_{g,h∈{1,β}} Z(g,h) = **⅓ Z[SO(16)×SO(16) on T⁶]**.

Every other sector rotates all three planes and is moduli-independent. Hence, at large T⁶ volume V (string units),

  **Λ₄ = (V/3) Λ₁₀ + O(V⁰) > 0**

for every one of the 104 models. The one-loop potential grows with the volume, and the dilaton tadpole
(∝ Λ) drives the dilaton towards weak coupling. **There is no stabilised vacuum at large volume.** The value
at the orbifold point needs the full twisted-sector computation, which is not done here.

## C. No escape through Bose–Fermi degeneracy

Exponentially suppressed cosmological constants arise in *interpolating* constructions (Itoyama–Taylor;
Abel–Dienes–Mavroudi), and they require equal numbers of massless bosons and fermions. Among the 104 models,
**n_B − n_F ∈ [−500, −64]**, and it is **never 0**; fermions always win. More generally, Groot Nibbelink et al.
(arXiv:1710.09237) proved that no symmetric non-supersymmetric toroidal orbifold has a one-loop Λ that vanishes
through sector-wise Killing spinors.

## Reading

The W(3,3) SO(16)×SO(16) models share the generic disease of the non-supersymmetric heterotic string: a positive
one-loop vacuum energy, linear in the internal volume, and a dilaton tadpole. Nothing specific to the W(3,3) shift
cures or worsens it at large volume. The leading term is model-independent, being ⅓ of the ten-dimensional value
times the volume. A viable vacuum would need an additional ingredient: stabilisation at small volume by the
twisted sectors, fluxes, or non-perturbative effects. None is computed here.
