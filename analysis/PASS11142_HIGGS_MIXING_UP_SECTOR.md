# Pass 11142 — no mixture of Higgs doublets splits charm from up

Producer: `analysis/w33_pass11142_higgs_mixing_up_sector.py`
Scans: `analysis/w33_pass11142_scan_all_doublets.py`, `analysis/w33_pass11142_scan_mixing_search.py`
Frozen: `data/w33_pass11142_all_doublet_textures.json`, `data/w33_pass11142_mixing_search.json`
Regression: `tests/test_w33_pass11142_higgs_mixing_up_sector.py`

Search: the light Higgs h = Hu_k + ε^a·Hu_j.
* Hu_k is a tree-level top doublet; Hu_j is any other doublet; a = 1, 2, 3.
* The condensate is included.
* The geometric factor of cubic entries is counted as small (g = 1) or O(1) (g = 0, as at T* = ρ, where it is ½).
* Models: 36621, 24165, 40521 and 2233 (model 57).

Result:
* **0 of 576 mixtures give an up hierarchy.** The exponents are always [0, 1, 1] or [0, 0, 0].
* **Why:** admixtures enter at order ≥ a ≥ 1, on top of a charm/up pair already fixed at tree level by the point
  reflection (Pass 11102). They shift the pair together.

Together with Passes 11102 (every vacuum alignment), 11114 (every direction in the tree triplet), 11127 (the condensate)
and 11133 (two VEV scales), this closes the up-sector hierarchy question at the level of selection rules. m_c = m_u at
leading order in every construction tried. Subleading O(ε^a) splittings exist but are not hierarchies.
