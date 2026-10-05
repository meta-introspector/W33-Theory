# Pass 11458 — the generic component of J₆'s blind spot is a genuine accident: no hidden symmetry, and no single achiral-eigenvector description

Producer: `analysis/w33_pass11458_generic_component.py` (main run; `--sym`)
Certificate: `data/w33_pass11458_generic_component.json`
Regression: `tests/test_w33_pass11456_11460.py`

**Background.** J₆'s non-reversible zeros on PU(3) have two 4-dimensional components (Passes 11369, 11419).
* The pseudo-reflection component is explained: it is the achiral hypersurface {h₆ = 0} (Passes 11419, 11434).
* For the generic-spectrum component, Pass 11434 observed that 79% of its points have an eigenvector with
  Δ₆ < 10⁻⁸. For Haar-random unitaries the rate is 13%.

## Refined points (Gauss–Newton on the full 729 × 729 residual A₃(U) − A₃(Uᵀ))

| seed | residual ‖F‖ | rank dF | eigenvector Δ₆ (sorted) | distance to reversible | J₈ |
|---|---|---|---|---|---|
| 113690016 | 8.0×10⁻¹⁵ | **4** | 6.8×10⁻⁸, 8.9×10⁻⁸, 3.2×10⁻⁷ | 0.125 | 0.139 |
| 113690072 | 6.6×10⁻¹⁵ | **4** | 4.2×10⁻⁹, 1.3×10⁻⁷, 3.8×10⁻⁶ | 0.319 | 0.851 |
| 113690092 | 7.1×10⁻¹⁵ | **4** | 1.4×10⁻⁸, 3.0×10⁻⁶, 8.8×10⁻⁶ | 0.399 | 0.088 |
| 113690129 | 6.8×10⁻¹⁵ | **4** | 2.9×10⁻⁹, 1.9×10⁻⁶, 1.2×10⁻⁵ | 0.585 | 1.473 |
| 113690141 | 8.5×10⁻¹⁵ | **4** | 3.3×10⁻¹¹, 3.7×10⁻⁹, 3.0×10⁻⁶ | 0.187 | 0.019 |
| 113690144 | 8.2×10⁻¹⁵ | **4** | **3.6×10⁻¹⁵**, 8.7×10⁻¹¹, 9.5×10⁻⁷ | 0.183 | 0.004 |

* The component is 4-dimensional at every one of the 6 refined points: dF has rank 4, i.e. codimension 4.
* **The eigenvector Δ₆ values do not move under refinement, and they are mixed.**
  * At four points no eigenvector is achiral: the smallest Δ₆ lies between 2.9×10⁻⁹ and 6.8×10⁻⁸.
  * At one point an eigenvector is achiral to working precision (Δ₆ = 3.6×10⁻¹⁵).
  * One point is borderline (3.3×10⁻¹¹).
* So "the generic component consists of unitaries with an achiral eigenvector" is **refuted as a characterisation**.
  It holds at some points and fails at others, which suggests the generic set has more than one piece. The 79%
  association of Pass 11434 is mostly small chirality, not zero chirality.

## No hidden symmetry

* **The idea.** Suppose a group G′ ⊃ Cl had the same degree-(3,3) commutant. Then A₃ would be G′-conjugation
  invariant, and every U with Uᵀ ∝ gUg⁻¹ for g ∈ G′ would have J₆ = 0. That would explain a blind spot as
  reversibility under a symmetry that third moments cannot tell apart from the Clifford group.
* **Continuous G′: excluded.** It would need some H ∈ su(3) with d/dε f_Y(e^{iεH}Ue^{−iεH}) = 0 for every spanning
  invariant f_Y(U) = avg_C|tr(Y†CUC†)|⁶. Over 80 random pairs (Y, U) the derivative map has **8 nonzero singular
  values**; the smallest relative value is 0.45. So there is no continuous hidden symmetry.
* **Finite G′: excluded.** The projective qutrit Clifford group is the Hessian group of order 216. By Blichfeldt's
  classification of finite subgroups of PSU(3), no finite primitive group strictly contains it.

**Reading.**
* The generic component is a genuine accident of the cubic moment map. The seven T-odd invariants of degree 3 vanish
  on a 4-dimensional set with no symmetry behind it and no achiral-eigenvector structure.
* Every such point is caught by J₈ (Pass 11435).
* Its explicit description remains open.
