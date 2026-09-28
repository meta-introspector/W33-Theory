# Pass 11103 — no alignment of the scalar vacuum splits charm from up

Producer: `analysis/w33_pass11103_reflection_breaking_vacua.py` (needs the field dumps; results frozen)
Frozen scans: `data/w33_pass11103_alignments_hidden_unbroken.json`, `data/w33_pass11103_alignments_hidden_broken.json`
Certificate: `data/w33_pass11103_reflection_breaking_vacua.json`
Regression: `tests/test_w33_pass11103_reflection_breaking_vacua.py`

## Background

Pass 11102 found charm and up at the same order in all 12 surviving SO(16)×SO(16) A8 models.
* **The protecting symmetry.** A point reflection of the family torus protects the degeneracy. It is the −I of the
  flavour group Δ(54) = H27 : ⟨−I⟩ (Pass 11105).
* **The loophole.** Pass 11102 gave every scalar a VEV of the same size. The twisted scalars sit at definite fixed
  points, though, and a vacuum can break the reflection spontaneously by giving them different VEVs.

## Method

1. **VEV weights.** Each scalar gets VEV ε^w. The weight w depends on the scalar's geometric fixed point g = l·c in
   the family torus, taken relative to the light Higgs point p0. Note that θ² class labels are inverted, so the class
   label c is not the point. Untwisted scalars get w = 1.

   | alignment | (w_p0, w_p1, w_p2) |
   |---|---|
   | sym | (1, 1, 1) |
   | p2_small | (1, 1, 2) |
   | split12 | (1, 1.5, 3) |
   | only_p0 | (1, –, –) |
   | only_p1 | (–, 1, –) |
   | p0p1 | (1, 1, –) |

   A dash means zero VEV.
2. **Entries.** Every up and down Yukawa entry gets its minimal weighted order Σ w_s e_s from an integer program under
   all selection rules. Cubic entries are exact and carry the geometric factor ε_geo = ε^g, for:
   * g = 1, a small torus area;
   * g = 4, a large area, where the VEV terms can dominate.
3. **VEV sets.** Two sets are scanned:
   * A: hidden singlets only, 120–126 scalars;
   * B: A plus the hidden-charged SM-neutral scalars, 126–138. B breaks the hidden group.

## Result

The scan covers 12 models × 3 light Higgs × (up, down) × 6 alignments × 2 regimes × 2 VEV sets.

| sector | parametric splitting of the two light generations |
|---|---|
| **up (charm/up)** | **never**: 0 cases in either VEV set |
| down (strange/down) | only for reflection-breaking alignments at large area (gap 1 or 3), in **exactly the 8 models with three extra d̄-like states** (24/36); never in the 4 minimal models |
| any, VEVs aligned at the Higgs point | never, as Schur's lemma requires (Pass 11105 D) |

**Why the up sector is protected.** The leading diagonal charm and up entries each need one scalar at the Higgs point
p0 and one at a non-Higgs point. They are reached at the same weighted cost because the θ and θ² scalars at a point
carry opposite class charges. With VEVs only off p0, the diagonal entries vanish. What remains is the off-diagonal
block, which is exactly degenerate by the reflection.

**The down-sector split** runs through mixing with the vector-like d̄ states: it appears in the 3×6 Q·d̄ block. The
vector-like mass matrix is not included here.

## Reading

The up-quark hierarchy m_c/m_u ≈ 600 cannot come from any alignment of the Standard-Model-neutral scalar vacuum in these
models, with or without the hidden group broken. It needs one of:
* an O(1)-coefficient tuning;
* physics outside the scalar vacuum, such as moduli-dependent Yukawa couplings or non-perturbative effects.

This is the Kobayashi–Nilles–Plöger–Raby–Ratz Δ(54) (hep-ph/0611020) acting on the survivors. The negative result for
spontaneous breaking is specific to their spectra.

Scope: orders and parametric exponents with random O(1) coefficients. The VEVs are assumed, not derived from a
potential.
