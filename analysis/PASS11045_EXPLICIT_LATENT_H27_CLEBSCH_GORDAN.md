# Pass 11045 — explicit commutant generators and Clebsch–Gordan compiler

Producer: `analysis/w33_pass11045_explicit_latent_h27_clebsch_gordan.py`
Certificate: `data/w33_pass11045_explicit_latent_h27_clebsch_gordan.json`

The latent `A9` action is made explicit in the multiplicity commutant:

`X_A = diag(1,1,1,X,X)`,
`Z_A = diag(1,omega,omega^2,Z,Z^2)`.

It commutes with the old internal-qutrit execution because the two actions occupy different tensor factors.

A closed-form 27×27 Clebsch–Gordan transform decomposes `A9 tensor V` into the regular Fourier sectors:

- scalar-character × `V` gives the three compatible `V` copies;
- `V × V` gives three `Vbar` copies by a difference-coordinate permutation;
- `Vbar × V` gives all nine one-dimensional characters by a difference coordinate plus qutrit Fourier transform.

The transform intertwines `X,Z,C` and has nonzero determinant at split Eisenstein primes 103 and 109.
