# Pass 11028 — phase-scanned interferometer sees the hidden signed clock

Producer: `analysis/w33_pass11028_signed_clock_interferometric_discriminator.py`
Certificate: `data/w33_pass11028_signed_clock_interferometric_discriminator.json`
Regression: `tests/test_w33_pass11028_signed_clock_interferometric_discriminator.py`

The central support element −I in GL₂(3) is invisible on the four projective
clock labels, so the coarse projective model treats it as the identity.

The exact 24-mode signed carrier does not. Let v_d be the uniform six-mode
state on clock direction d, and let z be the exact signed central operation.
The four signed normalized overlaps are

(+1/3, −1/3, +1/3, −1/3).

The signs depend on gauge/direction convention, but the phase-scanned
visibility does not:

|⟨v_d, z v_d⟩| / ⟨v_d, v_d⟩ = 1/3

for every one of the four clock directions.
An equal-arm identity-versus-z interferometer therefore predicts

signed carrier: Pmax = 2/3, Pmin = 1/3,

projective identity model: Pmax = 1, Pmin = 0.

Operationally z is a signed permutation of the 24 noncentral modes:
12 disjoint swaps, with 6 swap-pairs carrying a negative sign. A coherent
photonic network can represent the ideal operation by mode routing plus
calibrated π phase shifts.

The statistic is permutation-blind: S₄ relabelings may permute the four
fibres and exchange the signs, but the scanned visibility remains 1/3.

Pass 10952 supplies an independent numerical cross-check: the qutrit Clifford
representation of the same abstract central element has |tr U(−I)|/3 = 1/3.
That numerical equality is not an identification of representations.

This is an ideal coherent-mode discriminator, not a hardware error budget.
Loss imbalance, phase noise, state impurity, and detector visibility must be
calibrated in any experiment. The test distinguishes the exact signed carrier
from its coarse projective shadow; it is not a direct experimental test of the
full TOE.
