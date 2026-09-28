# Pass 11101 — a localized light Higgs singles out the top: 12 models survive every test so far

Producer: `analysis/w33_pass11101_localized_higgs_heavy_top.py`
Frozen results (all 104): `data/w33_pass11101_localized_higgs_104.json`
Certificate: `data/w33_pass11101_localized_higgs_heavy_top.json`
Regression: `tests/test_w33_pass11101_localized_higgs_heavy_top.py`

Pass 11098 found no single heavy top in any of the 104 SO(16)×SO(16) A8 models. It used a **generic** combination
of all Higgs doublets. In twisted sectors, though, every Higgs doublet sits at a definite fixed point, and the
physical light Higgs is one field. This pass repeats the analysis **Higgs by Higgs**:
* for each individual H_k, build the tree-level Yukawa matrix Y^(k) under the exact rules;
* take the untwisted couplings with their ε/δ structure;
* take the twisted couplings as O(1) when the three fields share the fixed point in every torus, and suppressed by
  ε_t = e^{−area_t} in each torus t where the fixed points differ;
* use generic, distinct ε_t (ε^(1, 1.7, 2.9)) and read the singular values' exponents at ε = 10⁻³ and 10⁻⁶.

## Results

| sector | pattern for **every** individual Higgs | models |
|---|---|---|
| up, untwisted | (0, 0, –): top = charm, whatever the Higgs | 73 |
| up, twisted | **(0, w_t, w_t)**: one O(1) top; charm and up both suppressed by the same torus | **31** |
| down | (0, w_t, w_t): single heavy bottom | 27 |
| down | no tree-level coupling | 76 |

The single-heavy-top models are exactly the twisted-up models of Pass 11098 (31 = 31).

**Why charm and up come out degenerate.** In each of these models the three generations sit at the three θ-fixed
points of one SU(3) torus. Those points are **pairwise equidistant**: all three minimal distances equal √(2/3) in
lattice units, checked here. The Z3 twist fixes the shape of the lattice, so no modulus can deform the triangle.
This two-valued structure of Z3 Yukawas is known (Casas, Muñoz, Ibáñez and others); only its application to these
models is new. Splitting m_c from m_u needs higher-order couplings, through the scalar VEVs of Pass 11097, or
mixing.

## The survivors

Combining Passes 11095, 11097 and 11101 leaves 12 models: **2, 10, 13, 14, 15, 35, 53, 57, 69, 77, 78, 102**. Each is:
* non-supersymmetric, tachyon-free and anomaly-free;
* exactly three chiral generations;
* one where every fractionally charged fermion can be made massive with the hidden group unbroken under all rules,
  at order **T* = 3** in every one of the 12;
* one where a localized Higgs gives a single heavy top **and** a single heavy bottom at tree level.

Still open for these 12:
* m_c = m_u at tree level;
* the positive one-loop vacuum energy and dilaton tadpole (Pass 11099);
* whether the required VEVs form a vacuum;
* the lepton sector, which has mixed patterns.
