# Pass 11114 — the quark-mass determinant is a Hesse cubic; at the stabilised point it degenerates, yet m₂/m_top ≥ ½ for every Higgs direction

Producer: `analysis/w33_pass11114_hesse_mass_determinant.py`
Certificate: `data/w33_pass11114_hesse_mass_determinant.json`
Regression: `tests/test_w33_pass11114_hesse_mass_determinant.py`

## The structure

In the survivors the three Higgs doublets form a Δ(54) triplet, one doublet at each fixed point of the family torus.
For a general light-Higgs direction h, the tree-level up-type matrix is

    M(h) = [[s h1, d h3, d h2], [d h3, s h2, d h1], [d h2, d h1, s h3]],

where s = Y_same(T\*) and d = Y_dist(T\*) are the A2 theta functions of Pass 11109. Couplings with two equal points and a
third different point are not collinear, hence forbidden. The determinant is

    det M(h) = h1h2h3 (s³ + 2d³) − s d² (h1³ + h2³ + h3³),

**a member of the Hesse pencil** x³ + y³ + z³ − 3λxyz with λ(T) = (s³ + 2d³)/(3 s d²).

That membership is forced by symmetry: det M(h) is Δ(27)-invariant up to a phase, and the Δ(27)-invariant plane cubics
are exactly the Hesse pencil (classical; Artebani and Dolgachev, *The Hesse pencil of plane cubic curves*, cited in the
corpus). The content computed here is the **value of λ(T)** and where it lands. The other track's Hesse configuration
(Passes 11061–11066; the September Hesse notes) is the geometry of the up-quark mass determinant over the Higgs
directions.

## At the stabilised point

At T\* = ρ, the one-loop minimum of Pass 11106 and 11110: **d/s = e^{iπ/3}/2 and λ = ω² exactly, so λ³ = 1**. The
determinant is a *singular* Hesse member, a triangle of three lines.

The spectrum on the standard alignments, as singular values relative to the largest:

| alignment | spectrum |
|---|---|
| (1,0,0), (1,1,1), (1,ω,ω²), … | (1, ½, ½) |
| (1,1,ω), (1,ω,1), (1,ω²,ω²) | (1, 1, 0): a massless up quark, but top = charm |

**Two results close the question:**
* Rank 1 (one heavy quark, two massless) would need M_ij² = M_ii M_jj, which forces (d/s)⁶ = 1. That is impossible for
  |d/s| ≠ 1.
* Numerically, min over all Higgs directions of σ₂/σ₁ = |d/s| (60 random starts × Nelder–Mead at ρ, i and 1.5i). The
  minimum is attained by the localised Higgs: ½ at ρ and (√3−1)/2 at i.

**No Higgs alignment beats a localised Higgs, and at the stabilised point the second quark is at least half as heavy
as the top.**

## Beyond one loop

T⁶/Z3 has no N = 2 subsectors, so the gauge threshold corrections do not depend on the Kähler moduli (Dixon,
Kaplunovsky, Louis). Hidden-sector condensates therefore cannot pull T\* away from ρ at this order. The family hierarchy
is closed for these vacua on every route tried: VEV alignment (Pass 11103), a larger family torus (Pass 11109/11110),
and Higgs alignment (here).
