# Pass 10980 — can flavour structure hide the regenerated R-parity violation?

Producer: `analysis/w33_pass10980_rpv_flavour_upper_bound.py`
Certificate: `data/w33_pass10980_rpv_flavour_upper_bound.json`
Regression: `tests/test_w33_pass10980_rpv_flavour_upper_bound.py`

## The idea

Pass 10969 showed that heavy superpartners cannot rescue the class: the Higgs quartic caps them near
3×10¹⁰ GeV. At that mass, though, the tree-level proton bound relaxes to
|λ′₁ⱼₖ λ″₁₁ₖ| ≲ 10⁻²⁷ (m̃/100 GeV)² ≈ **10⁻¹⁰**. Loop bounds on other flavour combinations become
irrelevant at such squark masses. So one escape remains: flavour texture might make the
first-generation projections of the FI-regenerated couplings tiny, even though some of them are of
order ⟨n⟩/M_s.

## Method

The 23 D-flat Z6-I models are analysed, each in the vacuum given by its Pass 10960 D-flat witness,
which contains the forced sneutrino-like singlet.

* **Rules.** Gauge + space group. Z6-I defines no R-rules, so these over-permit couplings, and every
  coupling found is an **upper bound** on its true size.
* **Coupling sizes.** Each coupling is ε^(order−3) (ε = 0.3) times a random O(1) coefficient, with
  orders from exact integer programs.
* **Light states.** Exotic vector-like pairs are removed through their mass matrices. H_u is the
  light bl direction; the first-generation q₁ and u^c₁ come from the up-Yukawa SVD.
* **Figure of merit.** P = max over d^c flavours and light lepton directions of |λ″(u₁dⱼdₖ)|·|λ′(q₁Ldₖ)|.
  This is conservative, because without H_d the down basis is undetermined.

## Result

| | models |
| --- | --- |
| no light H_u in the vacuum (μ-problem) | 7 |
| two massless exotic d-triplets (5 light d^c instead of 3) | 16 |
| P below 10⁻¹⁰ | **0** |

Across the 16 analysable models, P ranges from 2×10⁻³ to 0.34, at least seven orders of magnitude
above what the proton needs. These are upper bounds, so this alone does not prove that no texture could
work. It does show that no suppression mechanism appears. And the vacua are excluded anyway, independently
of proton decay: every one lacks a light Higgs or carries massless exotic triplets.

## Scope

* The analysis uses the Pass 10960 D-flat witness vacua only.
* Weak rules make P an upper bound.
* Coefficients are random O(1) numbers, not computed string couplings.
