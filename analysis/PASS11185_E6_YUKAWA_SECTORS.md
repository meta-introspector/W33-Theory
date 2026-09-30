# Pass 11185 — the E₆ Yukawa reading of splits and gate classes (a counting dictionary; no dynamics claimed)

Producer: `analysis/w33_pass11185_e6_yukawa_sectors.py`
Regression: `tests/test_w33_pass11185_e6_yukawa_sectors.py`

**Known ingredients** (Holotrade dictionary):
* The 27 complete factorisation frames are the 27 lines of the cubic surface, i.e. the 27 weights of E₆.
* Relative to one frame Φ they split 1 + 10 + 16, the SO(10) content of one E₆ generation.
* The 45 two-qutrit splits are the tritangent planes, i.e. the 45 monomials of the E₆ cubic invariant.

**Checked here (exact).**
* The 45 splits are exactly the 45 triangles of the frame graph (valency 10, the line-intersection graph).
* Relative to Φ:
  * the 5 splits through Φ are the **1·10·10** monomials (their other two frames lie in the 10);
  * the other 40 are the **16·16·10** monomials (one frame in the 10, two in the 16).

**Sector transitions by gate class** (per split):

| from | perfect (collinear) | I₃ = −1 (non-collinear) |
|---|---|---|
| 1·10·10 | 4 → 1·10·10, 8 → 16·16·10 | **0** → 1·10·10, 32 → 16·16·10 |
| 16·16·10 | 1 → 1·10·10, 11 → 16·16·10 | 4 → 1·10·10, 28 → 16·16·10 |

* At gate level, starting from a split F₀ in the singlet sector: the 13 824 perfect gates send it to the singlet sector
  4608 times and to the matter sector 9216 times; the 36 864 partially entangling gates always send it to the matter
  sector.
* **The frame stabiliser has order 1920 = |W(D₅)|, the SO(10) Weyl group** (960 in PSp(4,3), as in the Holotrade
  dictionary).
  * Relative to F₀, each of its elements is **either local (384) or perfect (1536), never partially entangling**: it maps
    F₀ to a split of the same frame, and those are pairwise collinear.
  * 768 of its elements fix no split at all: intrinsically entangling ticks inside W(D₅).

**Reading, at the level the computation supports.**
* A merely entangling tick never keeps a singlet-sector split in the singlet sector; only a perfect tick can.
* The ticks preserving an E₆ → SO(10) decomposition (the W(D₅) Weyl group) are, relative to the singlet-sector splits,
  always either local or maximally scrambling.
* No claim is made about fields, charges or Yukawa values.
