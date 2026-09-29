# Pass 11112 — a Standard-Model-neutral corridor: in models 53 and 57 the potential sends the neutral torus to its tachyon first

Producer: `analysis/w33_pass11112_neutral_tachyon_corridor.py`
Data:
* `data/w33_pass11112_bfield_degeneracy.json`
* `data/w33_pass11112_torus_resolved_tachyons.json`
* `data/w33_pass11112_pair_integrands_model57.json`, `..._model53.json`

Certificate: `data/w33_pass11112_neutral_tachyon_corridor.json`
Regression: `tests/test_w33_pass11112_neutral_tachyon_corridor.py`

Pass 11107 found that the winding tachyons are fractionally charged in 10 of the 12 survivors. Only models 53 and 57
also have Standard-Model-neutral ones (Y = 0, colour and SU(2) singlets, hidden singlets), and those are exactly
degenerate with the charged ones. A non-supersymmetric string has no D-terms, so gauge symmetry does not decide which
direction condenses.

## Findings

**A. The degeneracy survives any common B-field.** With both Wilson-line tori at x + iy (x = 0 … ½, three radii per
model, 36 points), the lowest neutral and charged levels agree to all digits.

**B. Its origin: the two kinds live on different tori.** The neutral tachyons wind only in Wilson-line torus 1, with
N = (±1, 0) or (±2, 0). The charged ones wind only in torus 2. Shrinking one torus at a time gives:

| shrink (other torus at 3i) | model 57 | model 53 |
|---|---|---|
| torus 1 to 1.4i / 1.2i / 1.0i | only neutral: Δ = −0.064 / −0.122 / −0.179 | only neutral: −0.096 / −0.154 / −0.211 |
| torus 2 likewise | only charged, same levels | only charged, same levels |

**C. The potential decides which torus reaches its critical radius first.** The kinetic metrics are equal (1/4y²).
Using the one-loop engine:

| model | Λ(torus 1 small) − Λ(torus 2 small) | exact strip + cap part | fitted tail part |
|---|---|---|---|
| 57, (1.6, 1.8) vs (1.8, 1.6) | **−1.29** | −1.21 | −0.08 |
| 57, (1.6, 2.0) vs (2.0, 1.6) | **−2.27** | −2.10 | −0.17 |
| 53, (1.8, 2.0) vs (2.0, 1.8) | **−0.94** | −0.85 | −0.09 |
| 53, (1.8, 2.2) vs (2.2, 1.8) | **−1.67** | −1.48 | −0.18 |

In both models **the potential pulls the neutral torus down faster.** The asymmetry comes from the exact part of the
modular integral, not the fitted tail.

## Reading

In both models that have neutral winding tachyons, the one-loop flow meets the small-radius instability first along a
Standard-Model-neutral direction. Condensation there breaks only extra U(1)s, and electromagnetism and the SM gauge group
survive. This is the first SM-preserving exit from the tachyonic boundary in the class.

The preference is model-specific: model 2, whose tachyons are all charged, prefers its other torus (Pass 11110). What is
not computed:
* the endpoint of the condensation, which needs the tachyon potential beyond quadratic order;
* whether the charged torus stays above its critical radius all the way there.

## Follow-up (Pass 11118)

The full gradient flow on a 6 × 6 grid reaches the neutral boundary from equal starting radii. The preference is
narrow, though: a starting asymmetry of about 1.5% in favour of shrinking the charged torus is enough to reverse it
(separatrix δ\* between 0.010 and 0.020).
