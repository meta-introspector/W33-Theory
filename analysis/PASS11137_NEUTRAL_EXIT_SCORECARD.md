# Pass 11137 — scorecard of the neutral-exit SO(16)×SO(16) W(3,3) vacua (Passes 11117–11136)

Producer: `analysis/w33_pass11137_neutral_exit_scorecard.py` (reads the certificates; no new computation)
Regression: `tests/test_w33_pass11137_neutral_exit_scorecard.py`

**Works**
* **The first instability is SM-neutral on the whole B circle, exactly.** The domains are hyperbolic disks at the Γ0(3)
  points; the neutral top lies at least 0.197 above the charged one (11136). The potential drives the moduli there
  steeply (11122), and B is nearly flat (11129).
* **Third-generation masses.**
  * A heavy top comes from one localised Higgs (11101).
  * With that single Higgs and the condensate, the bottom and tau sit at the same order ε² in 14/21 models (11132), so
    m_b/m_τ = O(1).
  * This needs ε_S·ε_T ≈ 0.0135.
* **The condensate never touches the up sector or μ** (11131).

**Fails or open**
* **Light generations.**
  * m_c = m_u and m_μ = m_e are protected by point reflection (11102, 11127).
  * The conjugate sector is anarchic at every ratio of VEV scales (11133).
* **The size of ⟨T⟩.**
  * It is not fixed by D-flatness: there is no FI mechanism without supersymmetry (11135).
  * The supersymmetric-analogue scale of 0.15–0.20 M_P is a coincidence, not a prediction.
  * Nor is it fixed by a tree-level quartic (11123).
* **The endpoint.** An exact winding law runs toward the Δ = −½ state of a T-dual string: one complex species, with no
  tachyon-free island (11128, 11134).
* **Vacuum energy.** It is positive at one loop (11106).
