# Pass 11119 — 21 survivors with a neutral exit, and 6 whose every nearby instability is Standard-Model-neutral

Producer: `analysis/w33_pass11119_neutral_tori_across_491.py`
Data:
* `data/w33_pass11119_torus_resolved_tachyons_491.json`
* `data/w33_pass11119_gauntlet_new_candidates.json`
* `data/w33_pass11119_diagonal_onset_21.json`
* `data/w33_pass11119_extended_moduli_8.json`

Certificate: `data/w33_pass11119_neutral_tori_across_491.json`
Regression: `tests/test_w33_pass11119_neutral_tori_across_491.py`

## 1. Which models put their neutral tachyons on a torus of their own? (all 491)

Each Wilson-line torus is shrunk alone to Im T = 1, with the other at 3i and the family torus at ρ, and the tachyons that
appear are classified:
* **60 models** have a torus that produces only Standard-Model-neutral tachyons;
* of these, **38 are neutral-only** (the other torus gives none or neutral) and **22 are separated** (the other torus
  gives charged ones, like 53 and 57);
* **all 60 have twisted quarks**, the heavy-top class.

## 2. The gauntlet (51 new candidates, full spectra from both engines)

The tests are the survivors' criteria:
* all fractional charges removable (hidden group unbroken, all rules; Pass 11097);
* a single heavy top from a localised Higgs;
* a single heavy bottom.

**19 new models pass. With 53 and 57, that makes 21 survivors whose Wilson-line instability has a neutral exit.**

## 3. The diagonal: both tori shrinking together (B = 0), as the potential drives them

| survivors | first tachyon on the diagonal | also at Im T = 1 |
|---|---|---|
| **8** (36621, 46043, 24165, 40521, 5904, 17224, 19688, 31386) | **neutral only** | neutral only |
| 11 (incl. 53, 57) | neutral and charged, degenerate | both |
| 2 (3818, 5273) | charged | both |

The onset is at Im T = √3 or 1.513.

## 4. How far does neutrality extend? (the 8)

The 8 were checked at the SU(3) point ρ, with the B-field on (0.3 + 1.2i, 0.5 + 1.3i, 0.2 + 0.9i), with unequal radii
(1, 2), (2, 1), (1.2, 0.9), and at smaller radius (0.8, 0.6).
* **6 models (36621, 46043, 24165, 40521, 5904, 17224) have only neutral tachyons along the whole B = 0 diagonal down to
  Im T = 0.6**, and at every nearby point.
* Charged tachyons appear only deep inside: with a B-field at small radius (0.2 + 0.9i) and at the SU(3) point ρ.
* In 19688 and 31386 charged ones also appear at Im T = 0.8.

## Reading

In these 6 survivors the first instability the one-loop potential leads to is **Standard-Model-neutral with a wide
margin**:
* the neutral onset is at Im T ≈ 1.5–1.7, while charged tachyons need Im T ≲ 0.9 together with a B-field;
* the potential slightly prefers B = 0 (Pass 11110).

Unlike model 57's 1.5% knife-edge (Pass 11118), the neutral exit here does not depend on initial conditions near the
diagonal. Its condensate breaks one U(1) and preserves the Standard Model (the structure found in Pass 11117 for 53 and
57). These models share the class-wide problems:
* the family degeneracy (Δ(54), m_c = m_u, m_{c,u}/m_t = ½ at T\* = ρ);
* the positive one-loop energy and the dilaton runaway.

The endpoint of their neutral condensation is the next computation.

Scope: tachyon classes are sampled at the listed moduli points, not on a continuum. The potentials of the 21 were not
recomputed; the direction of flow is taken from the class-wide results of Passes 11110 and 11118.
