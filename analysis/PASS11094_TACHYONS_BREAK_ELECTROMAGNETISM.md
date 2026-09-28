# Pass 11094 — every tachyon of the Z6-I twins is half-charged or coloured

Producer: `analysis/w33_pass11094_tachyons_break_electromagnetism.py`
Frozen orbifolder evidence: `data/w33_pass11094_tachyon_charges_orbifolder.json`
Certificate: `data/w33_pass11094_tachyons_break_electromagnetism.json`
Regression: `tests/test_w33_pass11094_tachyons_break_electromagnetism.py`

Pass 11092 showed that all 87 non-supersymmetric Z6-I W(3,3) twins are tachyonic. A tachyon is an instability,
not by itself a verdict: a tachyonic vacuum rolls somewhere. This pass asks what the tachyons carry, and
therefore what their condensation must break.

## A second tachyon engine

PATCH 7 of `analysis/w33_pass11092_orbifolder_n0_patches.py` gives orbifolder an optional mass level
μ < 0, read from `ORB_MASS_LEVEL`. Every non-identity sector is then solved at M²/8 = μ, and the full
projection machinery runs: centralisers, γ-phases and level matching. Patches 7b and 7c make a sector with no
right-movers at that level empty rather than fatal. The option is inert when the variable is unset; the
massless spectra of the controls are byte-identical. Three validations:
* **Z6-I twins, μ = −⅙.** The θ and θ⁵ sectors each carry exactly the Pass 11092 Python count, in **87/87**
  twins. There is nothing in any other sector, and every tachyon is a scalar.
* **Supersymmetric parents, μ = −⅙.** There are **0** twisted states.
* **A tachyonic control against an independent code.** The SO(16)×E8 string on T⁶/Z3 in standard embedding
  has a Witten shift of (0,0,0,1,0⁴; 0⁸). Both our mass-level engine and the non-SUSY orbifolder
  (arXiv:2504.20137) find **exactly one** tachyon, a (10,1,1) in the Witten sector.

## Quantum numbers — two independent computations

1. **orbifolder.** The tachyon fields at μ = −⅙ are listed, and each twin's Standard-Model hypercharge is
   carried into the plain orbifold U(1) basis by exact linear algebra. Hypercharge is orthogonal to every
   non-Abelian root, as checked.
2. **Pure Python, weight by weight.** Every tachyonic momentum p_sh of Pass 11092 (½p_sh² = 7/12 with
   N_L = 0, or 5/12 with one α₋₁/₆) is evaluated against the exact 16-dimensional hypercharge vector and the
   colour and SU(2)_L simple roots of the same model.

| | result |
|---|---|
| tachyon representations (orbifolder) | (1,1)_{±½} (263 + 199), (1,2)₀ (101), (3,1)/(3̄,1)_{±1/6} (33) |
| colourless tachyonic weights (Python) | **all** of charge ±½ (351 of +½, 267 of −½) |
| colourless *and* neutral tachyonic weights | **0** (Python), **0** fields with a neutral colourless component (orbifolder) |
| twins with a coloured tachyon | **27** in both engines |

## Reading

**Every tachyon of every Z6-I twin carries electric charge ±½, or colour.** None has a colourless, electrically
neutral component, so any condensate of these fields breaks U(1)_EM, and in 27 twins colour as well. The
twins cannot relax, through their own tachyons, into a vacuum that keeps the Standard-Model gauge group. With
Pass 11093 (their parents' excited exotics v, w, x become the tachyons) the picture is one mechanism. The
same half-charged exotic states that decorate every Z6-I W(3,3) Standard Model are, in the twin, lowered by
exactly ⅙ below the massless level.

Scope. This is a statement about the tachyonic fields at the orbifold point: a condensate of any of them breaks
electromagnetism. The endpoint of condensation (for example a resolved geometry) is not computed.
