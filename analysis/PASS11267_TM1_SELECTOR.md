# Pass 11267 — TM1 needs no tuned sign: a bilinear selects it for either sign, and boson loops fix the quartic's sign and size

Producer: `analysis/w33_pass11267_tm1_selector.py` (reuses the potential and vacua of Pass 11248 (parallel session) and
the mixing rule of Pass 11230)
Certificate: `data/w33_pass11267_tm1_selector.json`
Regression: `tests/test_w33_pass11265_11269.py`

**The problem.**
* In Pass 11248's sequestered S4 model, TM1 needs ε < 0 in ε(φ·χ)²; ε > 0 selects the excluded θ13 = 0 orientation.
* Pass 11254 showed that TM1 also needs |ε| ≲ 10⁻³.
* So both the sign and the size of ε looked like inputs.

## A. A sector-mixing bilinear selects TM1 for either sign (exact)

* On the sequestered vacuum manifold (φ a body diagonal ±d_i, χ a chord ±(d_i − d_j)), the S4 invariant φ·χ takes
  exactly the values **±4 on TM1 pairs and 0 on θ13 = 0 pairs**.
* So κ(φ·χ), **for either sign of κ**, lowers a TM1 pair by −4|κ| and leaves every θ13 = 0 pair at 0.
* Numerical multistart minimisation of Pass 11248's potential plus κ(φ·χ), seeded from every symmetric vacuum:
  * every tested minimum is TM1 for κ = ±0.001, ±0.005, ±0.02;
  * TM1 also wins in this scan against a wrong-sign ε = +5×10⁻⁴ once |κ| ≥ 0.005;
  * the vacuum tilt is 67°·|κ| for φ and 26°·|κ| for χ. With Pass 11254's alignment bound (~0.1°), this needs
    |κ| ≲ 1.5×10⁻³.
* **Price.** The bilinear exists only if the separate parities φ → −φ and χ → −χ are broken.

## B. With the parities exact, a loop fixes the sign by statistics

* Suppose the sectors are sequestered at a cutoff Λ, so ε(Λ) = 0. A messenger triplet with the S4-equivariant mass
  matrix M² = m₀² + g φφᵀ + g χχᵀ then generates ε through its Coleman–Weinberg potential.
* The eigenvalues depend on |φ|, |χ| and (φ·χ)², so the parities are kept.
* Comparing the TM1 and θ13 = 0 configurations at the vacuum norms:

| g | Λ/m₀ | boson messenger: ε_eff | fermion messenger: ε_eff |
|---|---|---|---|
| 0.1 | 10 | −1.4×10⁻⁴ → **TM1** | +1.4×10⁻⁴ → θ13 = 0 |
| 0.1 | 100 | −2.9×10⁻⁴ → **TM1** | +2.9×10⁻⁴ → θ13 = 0 |
| 0.3 | 10 | −1.3×10⁻³ → **TM1** | +1.3×10⁻³ → θ13 = 0 |
| 0.3 | 100 | −2.6×10⁻³ → **TM1** | +2.6×10⁻³ → θ13 = 0 |

**Reading.**
* **In this messenger model the sign of ε tracks statistics.** A bosonic messenger gives ε < 0 and hence TM1; a
  fermionic one gives the excluded orientation.
* **The size is natural.** A one-loop origin gives |ε| ~ g²ln(Λ²/m²)/(16π²)·(vev factors), i.e. 10⁻⁴ to 10⁻³ for
  g ~ 0.1–0.3. That is inside Pass 11254's window without tuning.
* Both the sign and the size that TM1 needs come from one assumption: **sequestered sectors communicating through a
  scalar messenger**.

**Scope.** This is a toy messenger with one equivariant coupling structure. The statistics sign is the standard
Coleman–Weinberg sign (bosons +, fermions −). The messenger is not identified with a W33 object here (OPEN).
