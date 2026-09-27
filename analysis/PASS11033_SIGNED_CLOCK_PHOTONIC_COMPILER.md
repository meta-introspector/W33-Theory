# Pass 11033 — exact 24-mode photonic compiler for the signed clock center

Producer: `analysis/w33_pass11033_signed_clock_photonic_compiler.py`
Certificate: `data/w33_pass11033_signed_clock_photonic_compiler.json`
Regression: `tests/test_w33_pass11033_signed_clock_photonic_compiler.py`

The central signed operation from Pass 11028 is an exact fixed-point-free involution on the 24 noncentral modes.

It compiles to:

- 12 disjoint mode swaps;
- 6 negative swap blocks;
- equivalently, an unsigned swap layer followed by 12 single-mode π phases.

With free physical placement, pair-adjacent modes give one simultaneous layer of 12 swaps plus one phase layer.

In the repository’s canonical linear ordering the unsigned permutation is exactly two reversals: a 16-mode reversal and an 8-mode reversal. The inversion number is 148, so a nearest-neighbor line needs exactly 148 adjacent swaps. Both reversals can run in parallel with minimum swap depth 15.

Each clock fibre has six interfering contributions whose signed sum is ±2, giving ideal visibility 1/3 versus 1 for the projective-identity model.

For bounded termwise phase error |δᵢ|≤ε,

V_signed ≤ 1/3 + 2 sin(ε/2),

while the projective cone bound gives V_projective ≥ cos ε.

The two envelopes remain disjoint for ε < 0.5338 rad ≈ 30.59°. This is a substantial ideal phase margin.

Boundary: that margin assumes equal amplitudes. Loss, source impurity, crosstalk, and detector imbalance still require a calibrated hardware likelihood model.
