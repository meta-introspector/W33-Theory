# Pass 11108 — tachyon-free over the reachable moduli ⇒ untwisted quarks ⇒ no heavy top (5/5)

Producer: `analysis/w33_pass11108_tachyon_free_moduli_rescan.py` (`rescan <files>` runs the enumerator on scan output)
Data:
* `data/w33_pass11108_tachyon_levels_104.json`
* `data/w33_pass11108_model_shifts_104.json`
* `data/w33_pass11108_rescan_levels.json`
* `data/w33_pass11108_rescan_candidates_gauntlet.json`

Certificate: `data/w33_pass11108_tachyon_free_moduli_rescan.json`
Regression: `tests/test_w33_pass11108_tachyon_free_moduli_rescan.py`

## Question

The one-loop potential shrinks the Wilson-line tori (Passes 11106, 11110), and below a critical radius a fractionally
charged winding state becomes tachyonic (Pass 11107). Which models stay tachyon-free where their moduli are sent? Every
model is run through the explicit enumerator of Pass 11107, with both Wilson-line tori at T_WL = ρ (the SU(3) point),
i, 1.5i and 2i, and the family torus at ρ.

## A. The 104 models of Pass 11095

| T_WL | tachyon-free |
|---|---|
| 2i | 104 |
| 1.5i | 55 |
| i | 55 |
| **ρ** | **1 (model 18)** |

* The lowest level at ρ takes only two values: **Δ = −1/18** (76 models) and **−1/6** (27).
* None of the 12 survivors is tachyon-free at ρ. Model 77 alone survives down to i on the imaginary axis.
* Model 18 fails the earlier tests: 48 fractionally charged states stay massless under all rules, its up-type Yukawas
  are untwisted (top = charm), and it has no down-type Yukawa.
* **Of the 31 models with twisted quark Yukawas, none is tachyon-free at ρ.**

## B. A fresh Wilson-line scan

Setup: the non-SUSY orbifolder, 8 × 50,000 tries on the 4 productive Witten-shift classes (seeds 20260941–20260952;
WSL `~/orb/p1109x/a8/rescan`).

The full scan (400,000 tries) gives 771 SM-like models, tachyon-free at large radius. **387 are inequivalent**:
deduplication by shift vectors agrees exactly with the orbifolder's own count.
* Levels at ρ: −1/18 in 277 models, −1/6 in 106, and **4 are tachyon-free** (A8SM_20260952_1033, _7375,
  A8SM_20260972_10118, A8SM_20260982_8666).
* Full spectra from both engines (patched orbifolder 1.2.1 levdump2; non-SUSY orbifolder nsosm) put all four through
  the gauntlet:

| model | quark doublets | up-type texture | tree-level down Yukawa | light fractional states (hidden unbroken / broken) |
|---|---|---|---|---|
| 1033 | untwisted ×3 | (0, 0, 2): top = charm | none | 48 / 12 |
| 7375 | untwisted ×3 | (0, 0, 2) | none | 24 / 0 |
| 8666 | untwisted ×3 | (0, 0, 2) | none | 48 / 12 |
| 10118 | untwisted ×3 | (0, 0, 2) | none | 48 / 12 |

## Reading

**Every model tachyon-free at the SU(3) point (5 of 5) has untwisted quarks**, and none of the 31 twisted-quark models
is tachyon-free there. The property that protects against the winding tachyon at small radius (untwisted families) is
the property that kills the single heavy top (Pass 11098: the ε-structure makes top = charm) and removes the down-type
Yukawa.

In this sample, tachyon-freedom over the reachable moduli and a realistic quark sector exclude each other. The
correlation is observed (5/5 and 0/31), not proven.

## Correction (Pass 11113)

The "0/31" count was a convenience sample. With the quark sectors of all 387 new models classified, the table is:

|  | twisted up-quarks | untwisted up-quarks |
|---|---|---|
| ρ-safe | 0 | 5 |
| tachyonic at ρ | 152 | 334 |

Fisher's one-sided p = 0.16. The correlation is suggestive but not significant: the ρ-safe rate is about 1% everywhere.
The Reading above ("exclude each other") is an over-read and is withdrawn. What stands is that no model found so far is
both ρ-safe and has twisted quarks.
