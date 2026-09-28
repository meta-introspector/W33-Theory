# Pass 11088 — baryonic and mixed hidden condensates do not rescue the Z3×Z3 vacua

Producer: `analysis/w33_pass11088_baryonic_hidden_condensates.py`
Engine (shared): `analysis/w33_fast_parity_realizability.py` (`realizable_choice`, new in this pass)
Certificate: `data/w33_pass11088_baryonic_hidden_condensates.json`
Regression: `tests/test_w33_pass11088_baryonic_hidden_condensates.py`

Pass 11024 searched the Z3×Z3 parity models for vacua in which the **curing factor** condenses. The curing
factor is the hidden gauge factor whose breaking alone makes the fractional-charge index vanish: SU(4) in c3
and c4, and SU(2)₁ in c1. That search used only **quadratic** condensates (mesons, SU(2) pairs, 6·6). It named
one residue: **baryonic** invariants, the ε-contractions of N distinct fundamentals. They were left out because
the discrete choice of which baryon condenses made the depth-first search explode. This pass runs the full
census.

## Method

* **Composites.** Every single-factor hidden invariant is in scope: mesons, baryons, antibaryons, SU(2) pairs
  and 6·6. That is 775 composites in c3 (190 mesons, 155 baryons, 430 antibaryons), 802 in c4, and 277 in c1.
* **Search.** Every FI-cancelling extreme ray that contains a condensate of the curing factor is tested for a
  realizable matter parity. The test is exact: a subset-sum dynamic programme over the finite group of
  residue vectors (`realizable_choice`). For fixed discrete group element, the relations that do not touch the
  composite types filter the grid. The rest reduce to reaching a target residue by choosing one variant per
  type. The programme agrees with the depth-first search of Pass 11024 on 99/99 rays (60 in c4, 39 in c3). It
  replaces hours with minutes.
* **Analysis of each realizable support.** Every support with an MSSM-viable parity that contains a baryonic
  condensate is analysed with the curing factor broken. The breaking is **full** (every multiplet splits into
  singlets) whenever a baryon, antibaryon or 6·6 of that factor condenses. It is **aligned**
  (SU(N) → SU(N−1)) if only mesons do. Centre (N-ality) charges enter as extra discrete charges, and mass ranks
  come from the exact monomial orders of Pass 11024.

Both simplifications can only overstate viability and masses. These are composite-even parity (a necessary
condition when a factor is broken by several condensates) and N-ality-allowed couplings under full breaking.
So a vacuum found to keep massless fractional charges is excluded without further assumptions.

## Results

| model | FI rays | realizable curing rays | with a baryonic condensate | massless fractional multiplets | of which symmetry-protected |
|---|---|---|---|---|---|
| c3 | 3069 | 78 | **39** (15 antibaryon only, 24 antibaryon + meson) | **148–180** | 112–162 |
| c4 | 4417 | 39 | **24** (15 antibaryon only, 9 antibaryon + meson) | **142–168** | 90–138 |
| c1 | 68 | 4 | **0**: the curing factor is SU(2), with no baryonic invariant | — | — |

* **Every one of the 63 baryonic vacua has an MSSM-viable parity, and every one keeps at least 142
  fractionally charged multiplets massless.** In every vacuum at least 90 of them are protected by an exact
  unbroken symmetry, since no integer exponents of any sign exist. The rest are forbidden by holomorphy.
* **Mixed condensates are covered.** "Antibaryon + meson" of the same factor is a mixed condensate. The census
  also realizes breaking of **both** SU(4)s at once (c4: 9 vacua, slots 0 and 1). No mixture helps.
* **Only antibaryons ever condense.** In both models no realizable support contains a baryon. This is an
  observation; the mechanism that selects the sign was not isolated.
* **The counts are not comparable with Pass 11024's.** Under full breaking a hidden 4-plet counts as four
  singlet multiplets, so ≥ 142 here and ≥ 120 under aligned breaking there are not a measure of "worse". What
  carries over is that neither count is ever zero.

## Reading

This closes the residue that Pass 11024 named. Across quadratic, baryonic and mixed condensates of one or both
hidden factors, **no parity-preserving vacuum on an FI-cancelling extreme ray with a curing condensate, in
any of the three Z3×Z3 parity models, gives the fractional charges masses**. Together with Pass 11087 (the protecting 3-group has no single removable piece),
Pass 11090 (the charge class is a space-group character, so no blow-up frees them), and Pass 11089 (the
obstruction does not depend on supersymmetry), the Z3×Z3 W(3,3) route is closed on established rules.

Scope. "Condensate" means a single-factor hidden invariant that takes a VEV, together with singlets, along an
**extreme ray** of the FI cone. Two things are not covered:
* A vacuum on a non-extreme face has a larger support. That makes the parity condition stronger, but it can
  also add mass terms.
* Invariants built from several hidden factors at once were not enumerated.
