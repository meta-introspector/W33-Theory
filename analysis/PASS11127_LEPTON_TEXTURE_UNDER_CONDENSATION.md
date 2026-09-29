# Pass 11127 — the unlocked entries leave the light-Higgs lepton hierarchy unchanged (m_μ = m_e survives), correcting the reading of Pass 11124

Producer: `analysis/w33_pass11127_lepton_texture_under_condensation.py`
Scan: `analysis/w33_pass11127_scan_textures.py` (Pass 11102 `order_matrix`: hidden-gauge check and exact cubic rules)
Frozen: `data/w33_pass11127_textures_six.json`
Regression: `tests/test_w33_pass11127_lepton_texture_under_condensation.py`

For every H_d of the six, the 8×3 lepton and 3×7 down order matrices were computed with and without the condensate T.
Findings:
* **The light doublets already have lepton Yukawas at tree level.** For H6–H8, the doublets with a tree-level down
  coupling, the pattern is one same-fixed-point cubic entry and two distinct-fixed-point entries. The exponents are
  [0, 1, 1] for both down and lepton, with and without T.
* **m_μ = m_e survives.** The cubic lepton entries of a light doublet are symmetric under the point reflection about the
  heavy generation's fixed point (18/18 doublet cases), the Pass 11102 mechanism behind m_c = m_u. The condensate is a
  Δ(54) singlet and enters only at order ≥ 2, so it cannot split them at leading order.
* **Exponents change only for H3–H5.** Down goes [3,3,3] → [2,2,2] in four models, and the lepton exponents drop
  likewise. No light doublet is affected.
* **The 11124 counts stand.** With the hidden check the new entries are lepton 27/27/81/108/81/81 and down
  0/0/81/81/81/81, identical to Pass 11124, which had omitted the check.

**Correction to Pass 11124.** Its reading, "a new source for exactly the Yukawas this class lacks (down/lepton)", was an
over-read. Those Yukawas exist at tree level for the light doublets. What the condensate unlocks are forbidden
*entries* of matrices that already have allowed ones. Pass 11130 shows where they do matter: for the conjugate of the top
Higgs. 11124's selection-rule statements stand: one insertion, q_T outside the singlet span, up and μ untouched.
