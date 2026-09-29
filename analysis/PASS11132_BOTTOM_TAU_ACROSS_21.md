# Pass 11132 — single light Higgs: the condensate lowers the bottom to ε² in exactly the 14 models with unlocked down entries

Producer: `analysis/w33_pass11132_bottom_tau_across_21.py` (shared with Pass 11133)
Scan: `analysis/w33_pass11132_scan_conjugate_textures.py`
Frozen: `data/w33_pass11132_conjugate_textures_21.json`
Regression: `tests/test_w33_pass11132_bottom_tau_across_21.py`

Pass 11130's single-Higgs analysis, extended to all 21 neutral-exit survivors:
* **Conjugation map.** Tree-level up doublet → its conjugate Hd, found by negating the U(1) charges and space-group
  classes. It is a bijection in 21/21.
* **Without the condensate** the bottom is at ε³ in all 21.
* **With it** the bottom drops to ε² in exactly the 14 models with unlocked down entries (Pass 11131). The two lists
  match exactly.
* **Tau at the same order.** In each of the 14, at least two of the three candidate light Higgs (all three in 24165)
  also put the tau at ε², so m_b/m_τ = O(1) as observed. With the other choice the tau sits at ε.
* **Scale.** ε_S·ε_T ≈ m_b/m_t ≈ 0.0135 (SM at 2·10¹⁶ GeV).
