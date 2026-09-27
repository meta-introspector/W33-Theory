# Pass 11047 — the cubic phase weld hits the exact dressed 54D target

Producer: `analysis/w33_pass11047_phase_weld_exact_54_right_inverse.py`
Certificate: `data/w33_pass11047_phase_weld_exact_54_right_inverse.json`

The diagonal E6 cubic weld was already known to surject onto an abstract 54D quotient. Pass 11046 identifies that quotient concretely as the two non-scalar latent sectors

`S2(27) + L(27)`.

Pass 11047 builds the full `S1 + S2 + L` K-Fourier matrix and expresses the cubic Jacobian in those coordinates.

For both diagonal weld orientations `c+p` and `c-p`:

- the quotient map onto `S2+L` has rank 54;
- the `S2` projection has rank 27;
- the `L` projection has rank 27;
- the same deterministic set of 54 input columns gives a nonzero 54×54 minor at both split primes 103 and 109;
- an explicit modular right inverse is verified.

Therefore the representation target and the cubic tangent mechanism are now aligned, not merely dimension-matched.

The remaining gap is physical synthesis: exponentiating this tangent compiler into a finite-time coherent unitary while preserving the FI orientation.
