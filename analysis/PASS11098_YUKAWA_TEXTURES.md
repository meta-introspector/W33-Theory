# Pass 11098 — tree-level Yukawas of the SO(16)×SO(16) A8 models: no single heavy top

Producer: `analysis/w33_pass11098_yukawa_textures.py`
Engine: `analysis/w33_so16_selection_rules.py`
Frozen results (all 104): `data/w33_pass11098_yukawa_104.json`; three sample models are re-run from their field
dumps (`data/w33_pass11097_so16_a8_sample_models.json`)
Certificate: `data/w33_pass11098_yukawa_textures.json`
Regression: `tests/test_w33_pass11098_yukawa_textures.py`

For the 104 tachyon-free three-generation models of Pass 11095 we take every renormalizable coupling ψψφ that the
**exact** tree-level rules allow:
* gauge invariance;
* the Witten Z₂, the point group and the space group;
* H-momentum conservation, with Σ Rᵃ = 0 exactly.

This covers the up-type couplings (q ū H_u), the down-type (q d̄ H_d) and the leptons (l ē H_d). We then ask what
mass pattern they give. Couplings are classified by structure, not just by existence.

* **Untwisted couplings** are the ten-dimensional gauge coupling restricted to the internal space. Under the
  internal SU(4) the fermions are 4's and the scalars 6's (4 ∧ 4). This gives two kinds of coupling:
  * triplet·triplet·scalar carries **ε_{ijk}**;
  * singlet·triplet_i·scalar_j carries **δ_{ij}**. Here the singlet is the fermion with internal momentum
    (−½, −½, −½).

  A matrix Σ_k ε_{ijk} h_k is antisymmetric, so its singular values are **(1, 1, 0) for every Higgs direction h**.
  This is checked on random h.
* **Twisted couplings** are O(1) when the three fields sit at the same fixed point in every torus. Otherwise they are
  suppressed by e^{−area} in each torus where the fixed points differ; only the O(1) ones are kept below.

An earlier estimate that put an independent random coefficient on every allowed triple gave rank 3 for the
untwisted up sector. That over-states it: the ε structure caps the rank at 2. This is the trap recorded as
"structural rank is only an upper bound". A second slip, treating the singlet-type couplings as zero, was caught
by printing the momenta of one model's lepton couplings.

## Results (O(1) couplings, random light-Higgs direction)

| sector | models | O(1) pattern |
|---|---|---|
| up, untwisted | **73** | singular values **(1, 1, 0)**: top and charm degenerate, up massless |
| up, twisted | 31 | rank 3, all O(1) (e.g. 1 : 0.64 : 0.62) |
| down | 77 / 27 | no tree-level coupling / rank 3 |
| lepton | 53 / 24 / 27 | none / **exactly one heavy lepton** / rank 3 |

**No model has a single heavy top at the renormalizable level.** In 73 models the top is exactly degenerate with
the charm at tree level, whatever the Higgs direction. In the other 31 all three up-type quarks are equally heavy.
The top–charm hierarchy (m_c/m_t ≈ 0.007) would have to come from large higher-order effects, or from mixing that
these fields do not have (there are exactly three q and three ū, and no vector-like copies).

What works:
* the down-type Yukawas vanish at tree level in 77 models, so m_b ≪ m_t is natural there;
* in 24 models exactly one charged lepton is heavy at tree level.

Scope. The couplings used are tree-level, at the orbifold point, with O(1) twisted coefficients drawn at random
where the fields share a fixed point. Higher-order couplings, which need the scalar VEVs of Pass 11097, and the
exact twisted-coupling values are not computed.
