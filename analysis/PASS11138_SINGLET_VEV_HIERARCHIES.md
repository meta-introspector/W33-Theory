# Pass 11138 — no pattern of singlet VEVs splits the light generations

Producer: `analysis/w33_pass11138_singlet_vev_hierarchies.py`
Scan: `analysis/w33_pass11138_scan_singlet_vevs.py`
Frozen: `data/w33_pass11138_singlet_vev_draws.json`
Regression: `tests/test_w33_pass11138_singlet_vev_hierarchies.py`

Pass 11133 gave the singlets one common scale and the condensate another. Here every singlet gets its **own** scale, a
random weight in [1, 3], in 12 draws. Setup: three models with the bottom at ε² (24165, 40521, model 57), each tree-level
top Higgs, and the single-Higgs (conjugate) down and lepton textures, condensate included.

* **Down quarks stay degenerate** in all 108 draws. The largest exponent spread is 0.12, i.e. about a factor 1.3 at
  ε = 0.1; a hierarchy would need a spread of order 1.
* **The two lighter charged leptons stay degenerate** (spread ≤ 0.10). The tau separates only for one Higgs choice (as
  in 11130).
* **So the conjugate-sector degeneracy is structural, not an artefact of equal VEVs.** The same singlet monomials
  compensate the generation label (the fixed-point class of the SM fields) for all three generations.

Note on the request: "fixed by their own D-terms" does not apply here. Pass 11135 showed there is no D-term potential
without supersymmetry, so the VEVs were treated as free and sampled.
