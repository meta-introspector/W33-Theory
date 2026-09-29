# Pass 11130 — with one light Higgs, the condensate brings the bottom and tau to the same order, ε², in four of the six

Producer: `analysis/w33_pass11130_single_higgs_bottom_tau.py`
Scan: `analysis/w33_pass11130_scan_higgs_conjugation.py`
Frozen: `data/w33_pass11130_higgs_conjugation_six.json` (+ the Pass 11127 textures)
Regression: `tests/test_w33_pass11130_single_higgs_bottom_tau.py`

**The conjugation map.** In all six the doublets with a tree-level up coupling are Hu6, Hu7 and Hu8, one per fixed point
of the family torus. Their complex conjugates are exactly Hd3, Hd5 and Hd4 (U(1) charges and space-group classes
negated). They are *not* the tree-level down doublets Hd6–8, which are conjugates of Hu0–2.

With a single light Higgs h = Hu_k (non-supersymmetric, so no holomorphy), the bottom and tau masses come from Q d^c h* and
L e^c h*, i.e. from the textures of conj(Hu_k):

| models | down (no T → T) | lepton, conj of Hu6 / Hu7 / Hu8 (no T → T) |
|---|---|---|
| 36621, 46043 | [3,3,3] → [3,3,3] | [1,3,3] / [2,2,3] / [2,2,3], unchanged |
| 24165 | [3,3,3] → [2,2,2] | [2,3,3] → [2,2,2] for all three |
| 40521, 5904, 17224 | [3,3,3] → [2,2,2] | [1,3,3] → [1,2,2] (Hu6); [2,3,3] → [2,2,2] (Hu7, Hu8) |

**Comparison with data.** Reference: the SM at 2·10¹⁶ GeV (Xing–Zhang–Zhou 2008), m_b/m_t ≈ 0.0135 and
m_τ/m_t ≈ 0.0228.
* **With the condensate:** n_b = n_τ = 2 needs ε ≈ 0.116 for the bottom and 0.151 for the tau. That is the same order, so
  m_b/m_τ = O(1), as observed (0.59). It requires ⟨T⟩⟨S⟩/M² ≈ 0.014–0.023.
* **Without it:** the bottom sits at ε³ (ε ≈ 0.24), while the tau sits at ε² or ε.
* **Where it fails:** in each sector all three singular values share one exponent, so m_s ~ m_d ~ m_b at this order. The
  light generations are not reproduced: the Pass 11102 degeneracy recurs in the conjugate sector.

Scope:
* Exponents only, with random O(1) coefficients and one common ε for S and T insertions.
* ⟨T⟩ is not computed (Pass 11123: beyond perturbation theory).
* The 2·10¹⁶ GeV values are an order-of-magnitude reference.
