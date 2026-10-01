# Pass 11248: in a sequestered two-flavon model, the sign of one coupling selects TM₁

Producer: `analysis/w33_pass11248_sequestered_tm1.py`
Certificate: `data/w33_pass11248_sequestered_tm1.json`
Regression: `tests/test_w33_pass11248_sequestered_tm1.py`

**Model.** An S₄-invariant toy model in the triplet 3′ of the line stabiliser (Pass 11243). The charged-lepton
potential's minimum is the point direction d_X; the neutrino potential's minimum is the chord d_X − d_Q. The neutrino
potential needs degree-8 terms, since no renormalisable one-triplet potential has the chord as its minimum (Pass
11230). The two sectors couple through λ|φ|²|χ|² + ε(φ·χ)², and the potential is minimised globally from all 96
orbit pairs.

**Results.**
* **ε < 0 selects TM₁; ε > 0 selects the excluded θ₁₃ = 0 orientation.** This holds at every ε tested, from ±0.005 to
  ±0.08. It confirms Pass 11243's tie-breaker, (φ·χ)² = 16 against 0.
* **Size of the misalignment.** The cross term tilts both flavons away from their symmetric directions at first
  order: about 84°·|ε| for φ and 33°·|ε| for χ at these couplings (b = 0.3, k = 2, h = 3, λ = 0.2). For ε = −0.005
  that is 0.42° and 0.17°.
* **What this means for the mixing.** The TM₁ column is predicted up to corrections linear in ε. This model does not
  fix the size of ε. The values of θ₁₃ and δ beyond TM₁ require a full Yukawa model: OPEN.

**Reading.** The vacuum-alignment problem for TM₁ in W33 now has a minimal concrete solution: sequestered sectors, a
non-renormalisable chord potential, and one negative cross-coupling. These are inputs the geometry does not derive.
