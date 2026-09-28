# Pass 11095 — an independent code checks the N = 0 engine, and the W(3,3) A8 shift has tachyon-free non-supersymmetric three-generation completions

Producer: `analysis/w33_pass11095_so16xso16_a8_standard_models.py`
Drivers: `analysis/orbifolder_n0_drivers/` (`nsoscan.cpp`, `nsosm.cpp`, `nsodump.cpp` for the non-SUSY orbifolder;
`levdump.cpp`, `anomdump.cpp` for ours; `nso2orb.py` converts between the two formats)
Frozen evidence: `data/w33_pass11095_nonsusy_orbifolder_crosscheck.json`, `data/w33_pass11095_a8_scan_evidence.json`
(including all 104 model definitions)
Certificate: `data/w33_pass11095_so16xso16_a8_standard_models.json`
Regression: `tests/test_w33_pass11095_so16xso16_a8_standard_models.py`

## A. The non-SUSY orbifolder as an independent check

The non-SUSY orbifolder (Escalante-Notario, Pérez-Martínez, Ramos-Sánchez, Vaudrevange, arXiv:2504.20137;
github.com/StringsIFUNAM/nonSUSYorbifolder) builds orbifolds of the SO(16)×SO(16) string. The Witten twist
v₀ = (0,1,1,1) with its gauge shift V₀ is one Z₂ generator, and each point-group twist is supersymmetric. It
was compiled here without its readline prompt and driven directly.
* **Massless spectra: 19/19 identical.** All 19 Z₂W × Z_N sample models it ships (Z3, Z4, Z6-I, Z6-II, Z7,
  Z8-I, Z8-II, Z12-I, Z12-II) were converted to orbifolder-1.2.1 format and run through our patched engine. The
  gauge groups and the field-by-field massless spectra agree in every case (434–828 fields; graviton fields,
  which orbifolder 1.2.1 does not list, are excluded).
* **The comparison found a sixth fault.** Orbifolder 1.2.1 gives a twist its geometric order: 3 for
  (0, 4/3, 4/3, 1/3). On fermions the order is 6, because 3v has odd sum. The wrong order made the local
  modular-invariance test drop the Witten-twisted sectors. PATCH 6 fixes it, and 6b keeps the pure Witten
  sector, which rotates no plane, from being treated as a twisted sector without oscillators. The Z6-I twins
  are unaffected (6v′ has even sum), and their spectra are byte-identical.
* **Prior art.** Their code already contains our patch 3, the right-moving oscillator sign, with the comment
  "wrong trafo sign for R-moving oscillator excitations in original Orbifolder". Pass 11092's fix is a
  rediscovery and is credited to them.
* **Tachyon control.** SO(16)×E8 on T⁶/Z3 in standard embedding: both codes find exactly one tachyon, a
  (10,1,1) in the Witten sector.
* **Limits of their code.** It cannot load the E8×E8 twins of Passes 11092–11096. Its right-mover lattice is
  hard-coded for the Witten-twist structure, and the load fails with the same "not invariant under
  constructing element" error.

## B. SO(16)×SO(16) completions of the W(3,3) shift (exact)

The W(3,3) twist of T⁶/Z3 is the A8 Kac pair V₁ = (⅙⁷, ⅚; 0⁷, ⅔) (Pass 10979). An SO(16)×SO(16) version needs a
Witten shift V₀ = (a; b) with a and b in ½·{norm-4 vectors of E8}. Those form one Weyl orbit of 2160 vectors, so
every choice gives SO(16)×SO(16) in ten dimensions. The shift must also satisfy **V₀·V₁ ∈ Z**.
* **793,508** of the 4,665,600 pairs are admissible.
* They fall into **six gauge-group classes** before Wilson lines:

  | first E8 | second E8 | Witten shifts |
  |---|---|---|
  | SU(5)×SU(4) | SO(8)×SU(4) | 431,200 |
  | SU(5)×SU(4) | SU(7) | 293,888 |
  | SU(8) | SU(7) | 57,344 |
  | SU(5)×SU(4) | SO(14) | 6,580 |
  | SU(8) | SO(8)×SU(4) | 4,480 |
  | SU(8) | SO(14) | 16 |

  (Each class also has U(1) factors.)
* **The Witten shift always breaks the SU(9)** of the A8 twist, because SU(9) ⊄ SO(16).
* The non-SUSY orbifolder accepts four classes, with exactly the predicted gauge groups (a cross-check of the
  enumeration). It rejects the two classes containing SO(14) as "orbifold group ill-defined".

## C. The scan: 104 tachyon-free three-generation models

20,000 random Wilson-line draws were made on one V₀ from each loadable class, with V₀ and V₁ fixed (80,000
models).
* **Every model is tachyon-free** (80,000/80,000).
* **153 draws pass the Standard-Model test** of the non-SUSY orbifolder, which tests SM-like chiral content, a
  Higgs doublet and vector-like exotics. These are **104 inequivalent models**: 16, 34, 31 and 23 per class.

**Every one of the 104 was re-derived with our independent engine.**

| check | result |
|---|---|
| tachyons, from our mass-level engine at the only tachyonic level (M²/8 = −½, Witten sector) | **0** in all 104 |
| anomalies (our checker) | all 104 consistent: 101 with one Green–Schwarz U(1), 3 with **no** anomalous U(1) |
| net chiral content, **their** hypercharge (Y·Y = ⅚) applied to **our** fermions | **exactly three generations** (3 × [(3,2)_{1/6} + (3̄,1)_{−2/3} + (3̄,1)_{1/3} + (1,2)_{−1/2} + (1,1)₁]) in 104/104; 23 of them appear with 3 ↔ 3̄ relabelled |
| net chiral fractionally charged fermions | **none** in 104/104 |
| Higgs-doublet scalar fields | 6–42 per model |

## What these models are, and are not

**They are non-supersymmetric, tachyon-free, anomaly-free string models with exactly three chiral generations,
built on the W(3,3) twist of T⁶/Z3 and its A8 shift.** Without supersymmetry there are no dimension-4 or
dimension-5 proton-decay operators, the obstruction that closed the supersymmetric routes.

Not shown:
* **Masses for the exotics.** Every model carries **88–230 vector-like fractionally charged fermion states**.
  Vector-like means masses are allowed in principle, but they must be generated by couplings to scalar VEVs.
  In Z3×Z3 the analogous states proved unremovable (Passes 11024, 11087–11090). This is the first thing to
  test.
* **A stable vacuum.** A non-supersymmetric string has a one-loop cosmological constant and dilaton tadpole.
* **Yukawa couplings.**

Prior art. Standard-Model-like tachyon-free SO(16)×SO(16) models are known in large numbers: more than
170,000 over the 138 Abelian geometries (Pérez-Martínez, Ramos-Sánchez, Vaudrevange, arXiv:2105.03460). What is
new here is only that the W(3,3) A8 shift admits them.
