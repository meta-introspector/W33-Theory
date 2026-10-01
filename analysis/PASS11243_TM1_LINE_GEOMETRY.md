# Pass 11243: TM₁ is a configuration of points on a line of W33

Producer: `analysis/w33_pass11243_tm1_line_geometry.py` (builds W33 and Sp(4,3) from scratch; reuses Pass 11230's
invariant-potential sampler)
Certificate: `data/w33_pass11243_tm1_line_geometry.json`
Regression: `tests/test_w33_pass11243_tm1_line_geometry.py`

**Question.** Pass 11225 traced the TM₁ column to the S₄ of the line stabiliser 3³:S₄. Pass 11230 showed that no
generic S₄ flavon potential selects TM₁ (0/1500). This pass asks:
* which object of W33 carries the TM₁ residual symmetries;
* why potentials miss it.

## A. The line stabiliser, built from scratch

* W33 is the 40 points of PG(3,3) with the symplectic form.
* Sp(4,3) is enumerated from its 40 transvections and has 51840 elements.
* The stabiliser of a totally isotropic line L has order 1296. It acts on the 4 points of L as the **full S₄**, with
  kernel of order 54 (±1 times 3³). This is 3³:S₄ in PSp(4,3).

## B. TM₁ in the geometry

The triplet **3** comes from the 4 points of L sent to the 4 body diagonals d_i of a cube; **3′** is **3** ⊗ sign.
Columns are computed relative to the charged-lepton Z₃, the 3-cycles fixing one point X of L.

| neutrino involution | column |
|---|---|
| transposition (XQ), which **moves X** | **(2/3, 1/6, 1/6) = TM₁** (3 involutions) |
| transposition fixing X | (1/2, 1/2, 0): θ₁₃ = 0, excluded |
| double transposition | (1/3, 1/3, 1/3) = TM₂ |

* **TM₁ ⟺ the neutrino symmetry swaps the charged-lepton point X with a second point Q of the same line**, and fixes
  the other two points R and S.
* In **3′** the vacuum directions are geometric:
  * the charged-lepton flavon points along the point itself, φ_e ∝ d_X, whose stabiliser is exactly Z₃;
  * the neutrino flavon points along the chord, φ_ν ∝ d_X − d_Q, whose stabiliser is exactly ⟨(XQ)⟩.
* In coordinates, d_X − d_Q = (0, 2, 2). This is the TM₁ direction (0,1,1) of 3′ found in Pass 11230.

## C. Why generic potentials miss it

* Z₃ and ⟨(XQ)⟩ meet trivially. Including the sign flips that every quartic potential respects, the pair
  (d_X, d_X − d_Q) is fixed by a single element, (RS) × (−1, −1).
* That element forces only part of the gradient of an invariant potential to vanish. Over 60 random invariant
  potentials, the gradient component transverse to the two radial directions spans **2 dimensions**.
* So TM₁ is a critical point only on a **codimension-2** subvariety of couplings, which has measure zero. That is why
  Pass 11230 found 0/1500.
* For the aligned pair (d_X, d_X) the joint stabiliser contains Z₃, the transverse gradient vanishes identically
  (rank 0), and aligned vacua are automatic. Pass 11230 found exactly those.

## D. The sequestered limit and the tie-breaker

Suppose the two sectors couple only through |φ|²|χ|² (a sequestered limit). Then the vacuum manifold is
orbit(d_X) × orbit(d_X − d_Q), with 8 × 12 = 96 points:
* 48 of them are TM₁;
* 48 are the excluded θ₁₃ = 0 column.

The two orientation-sensitive quartic cross invariants on these two classes are:

| class | (φ·χ)² | Σφᵢ²χᵢ² |
|---|---|---|
| TM₁ | **16** | 8 |
| θ₁₃ = 0 | **0** | 8 |

So at first order in a small cross-coupling ε(φ·χ)², **the sign of ε alone decides**:
* **ε < 0 selects TM₁**, the chord sharing the charged-lepton point;
* ε > 0 selects the excluded θ₁₃ = 0 orientation.

The term also pulls φ toward χ at O(ε). This is the codimension-2 obstruction of C showing up as a small deviation
from exact TM₁ at that order. Here the deviation is identified, not computed.

## Reading

* The geometry does supply the alignment, as an incidence statement. The charged leptons see a point of a line (a
  measurement context, Pass 11225), and the neutrinos see the chord from that point to a second point of the same
  line.
* What it does not supply is the potential that prefers this configuration.
* The minimal extra input is one sign, ε < 0 for (φ·χ)², in a sequestered two-flavon potential. With it, TM₁ is
  selected over its only degenerate competitor, up to O(ε) corrections.
* Status of the mixing angles: still **OPEN**. The obstruction is reduced from "no potential selects TM₁" to "one
  sign plus sequestering".

**Prior art.**
* Trimaximal TM₁ from S₄ with a Z₃ / Z₂ split is standard (Albright–Rodejohann; Luhn; King and collaborators, with
  the (0,1,−1) alignment of CSD-type models).
* The vacuum-alignment problem and its driving-field solutions are standard (Altarelli–Feruglio).
* New here:
  * the incidence reading on a W33 line (point plus chord);
  * the codimension-2 count with its explicit stabiliser;
  * the sequestered 48/48 tie and its single-sign breaker.
