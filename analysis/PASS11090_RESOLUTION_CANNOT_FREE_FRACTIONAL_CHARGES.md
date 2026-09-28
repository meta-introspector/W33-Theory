# Pass 11090 — blowing up the fixed points cannot free the massless fractional charges

Producer: `analysis/w33_pass11090_resolution_cannot_free_fractional_charges.py`
Certificate: `data/w33_pass11090_resolution_cannot_free_fractional_charges.json`
Regression: `tests/test_w33_pass11090_resolution_cannot_free_fractional_charges.py`
Builds on: Passes 11024 (vacua, exact solver) and 11087 (symmetry blocks).

A vacuum in which twisted singlets take values is the orbifold limit of a resolved compactification: the
singlet VEVs are the blow-up moduli of the fixed points they live at. The question is whether the smooth
geometry could give the massless fractional charges a mass. In the smooth geometry the fixed-point
labels lose their meaning. So the space-group selection rules need not survive, and nor need the
point-group sector labels or the discrete R-rotations. The gauge charges, the U(1)s and the hidden
centres, always do. Dropping selection rules can only **add** couplings. The massless count under a
weaker rule set is therefore a lower bound for every resolution that keeps the rules still imposed.

## Result (exact holomorphic mass ranks, fractionally charged classes)

| massless fractional multiplets | all rules | no space group | no space/point group | **gauge charges only** |
|---|---|---|---|---|
| 55 hidden-condensate parity vacua | 120–166 | 108–156 | 108–156 | **18–130** |
| 28 singlet parity vacua | 93–115 | 81–105 | 81–105 | **15–69** |

* **No resolution frees them.** Even if a resolution destroyed every discrete selection rule, keeping
  only gauge invariance, every one of the 83 vacua would still keep at least 15 fractionally charged
  multiplets massless: at least 18 in the condensate vacua, at least 15 in the singlet vacua.
## Why: charge quantization is a character of the space group

For every state, take the charge class **c = 3Q + t (mod 3)**, where t is the colour triality. It is 0
for every Standard-Model state, and nonzero exactly for the fractionally charged ones. Across all 13
Z3×Z3 Standard Models:
* **c is constant on every fixed point.** This holds on all 946 fixed-point labels (sector and
  translation labels).
* **c is a linear Z₃ character of the space group,** unique in each model.
  * In 9 models it reads a single translation label: a Wilson line along one torus.
  * In c1_2822 it reads the sector labels.
  * In the two c4 models it reads a sector and a translation label together.
* **The untwisted sector has c = 0.** 516 of the 946 fixed points have c ≠ 0.

The consequences follow at once:
* Every state at a fixed point with c ≠ 0 is fractionally charged.
* Such a fixed point hosts no neutral field, so no SM-singlet modulus lives there. The count of SM
  singlets at fractional fixed points is 0 in all 13 models.
* **The fixed points that carry fractional charges can be blown up only by VEVs that break colour or
  electromagnetism.**

This is the orbifold form of the Wen–Witten observation that fractional charges come from Wilson-line
sectors (X.-G. Wen and E. Witten, Nucl. Phys. B261 (1985) 651). Here it is made exact and used as a
no-go for resolution.

Control: the vacuum fields and the fractional states share twisted sectors. Sectors (0,⅓), (⅓,⅓), (⅔,⅓)
and (⅓,⅔) host both. They never share a fixed point, and the character explains why. A first draft of
this pass claimed that the states with no allowed coupling "all sit at unresolved fixed points". The
control showed that statement to be vacuous, since no fractional state sits at a resolved fixed point at
all. It is replaced by the theorem above.

Together with Pass 11087, this answers item 4 of the Pass 11024 follow-up. The space-group selection
rules are necessary for only 16% of the forbidden mass terms. Dropping them, and every other discrete
rule, lowers the massless count but never to zero.

Scope: the rule-dropping bound covers every resolution that preserves the gauge group. It says nothing
about resolutions that also break the gauge group further. Those would need a new vacuum, including
hidden-charged moduli (Passes 11024 D and 11088).
