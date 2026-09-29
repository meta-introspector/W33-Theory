# Pass 11133 — no ratio of the two VEV scales splits the light generations: the conjugate sector is anarchic

Producer: `analysis/w33_pass11132_bottom_tau_across_21.py`
Regression: `tests/test_w33_pass11133_two_scale_anarchy.py`

The condensate T and the singlets S are independent scales. With r = log ε_T / log ε_S = 0.5, 1, 2, 3, each entry of the
conjugate-sector matrices was minimised as n_S + r·n_T.

* **One option set for every entry.** The lowered entries use either (n_S, n_T) = (1, 1) or three singlets and no
  condensate, and every entry of the 3×3 down block has the same options. The three singular values therefore always
  share one exponent, min(1 + r, 3).
* **No splitting anywhere.** There are 0 cases with three distinct exponents, across 21 models × 3 Higgs × 2 sectors ×
  4 values of r.

Reading:
* The conjugate sector is anarchic. It sets the third-generation scale but not the hierarchy among generations.
* The light-doublet sector gives [0, 1, 1] with m_c = m_u and m_μ = m_e (Passes 11102, 11127).
* So neither sector produces m_d ≪ m_s ≪ m_b.
* A hierarchy would need VEVs that differ *among the singlets themselves*, or couplings the exponent counting cannot see.
