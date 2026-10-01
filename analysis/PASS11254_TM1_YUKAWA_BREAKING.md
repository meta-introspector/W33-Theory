# Pass 11254 — a Yukawa model on the line stabiliser's S4: TM1 needs its flavon vacua aligned to ~0.1°

Producer: `analysis/w33_pass11254_tm1_yukawa_breaking.py` (reuses Pass 11248's potential and vacua and Pass 11230's
sum rule)
Certificate: `data/w33_pass11254_tm1_yukawa_breaking.json`
Regression: `tests/test_w33_pass11250_11254.py`

**Background.**
* Pass 11243: TM1 = the neutrino symmetry swaps the charged-lepton point X with a second point Q of the same W33 line.
* Pass 11248: the sign of ε in ε(φ·χ)² selects TM1, and the flavons tilt at O(ε).

**Scope first.** In residual-symmetry TM1, **θ13 is not a function of ε**.
* At ε = 0 the Z₂ fixes one column. θ13, θ23 and δ are set by the Yukawa couplings.
* Confirmed here: across 60 fitted models at ε = 0, with every one exactly TM1, sin²θ23 spans [0.29, 0.71] and cos δ
  spans [−1.00, 1.00].
* So the question "θ13 and δ as functions of ε" (Pass 11248's open item) is replaced by the one that is well posed:
  **how fast does ε break the TM1 relations?**

## Model

L, φ and χ are all in the 3′ of S4. Every term is checked to be S4-equivariant.

* **Charged leptons.** H_e = M_e†M_e = α|φ|² + βφφᵀ + γ diag(φ²) + iκ₁[φ]_× + iκ₃[φ³]_×.
  * The antisymmetric terms split m_μ from m_τ under the residual Z₃. α, β, κ₁ are solved for the physical masses,
    with the electron on the Z₃ eigenvector that carries the TM1 weight 2/3.
* **Neutrinos.** M_ν = a + b|χ|² + c χχᵀ + d diag(χ²) + e K(χ), with complex couplings.
  * The couplings are fitted to NuFIT 6.0: sin²θ13, Δm²21/Δm²31, normal ordering, and m₁² ≤ 0.15 Δm²31.
* **An accidental symmetry, found on the way.** With only quadratic terms, M_ν at the chord vacuum commutes with
  diag(−1, 1, 1). That reflection is not in S4, and it freezes the whole PMNS matrix (no fit is possible).
  * The fix is K(χ): K₁₂ = z(y² − x²), K₁₃ = y(x² − z²), K₂₃ = x(z² − y²).
  * K is the **unique** S4-equivariant symmetric cubic: the Reynolds projection of all 60 cubic monomial maps has
    rank 1.
* This is a representative low-degree family, not the most general model.

## Results (60 models, the exact ε-tilted vacua of Pass 11248)

| quantity | value |
|---|---|
| TM1 at ε = 0 | exact for every model (\|U_e1\|² = 2/3, sin²θ12 = 0.3184) |
| vacuum tilt per unit ε (Pass 11248 potential) | φ 84°, χ 33° (reproduces Pass 11248) |
| d\|U_e1\|²/dθ, neutrino side only (per radian of χ tilt) | median **13**, 10–90% range 4–27; linear in 100% of models |
| d\|U_e1\|²/dθ, charged-lepton side only (per radian of φ tilt) | median **17**, 10–90% range 3.5–46; linear in 100% |
| fraction of models with sin²θ12 inside NuFIT 3σ | 83% at ε = −5×10⁻⁴, 57% at −10⁻³, 35% at −2.5×10⁻³, 10% at −0.02 |

**Reading.**
* The TM1 prediction survives the 3σ sin²θ12 window only if each flavon vacuum is aligned to about
  **0.02/15 rad ≈ 0.08°**.
* In Pass 11248's potential that means **|ε| ≲ 10⁻³**. The sign of ε selects TM1, but its size must be small. This is
  a quantified fine-tuning, not a prediction.
* **Hypothesis tested and rejected.** I had expected the charged-lepton side to be amplified by m_τ²/m_μ² ≈ 283. The
  measured median is 17/rad, so it is not.
* **The TM1 sum rule.** The relation tying cos δ to θ13 and θ23 is broken much faster: the median is small, but the
  10–90% range is ±300 per unit ε. That is because cos δ is a ratio of small quantities.

**Open.**
* A dynamical reason for |ε| ≲ 10⁻³, for example ε generated only at loop level.
* A complete model with shaping symmetries, where the charged-lepton rows come from distinct operators.
