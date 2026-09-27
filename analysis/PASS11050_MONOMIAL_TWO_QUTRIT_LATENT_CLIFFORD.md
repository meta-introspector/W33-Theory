# Pass 11050 — the latent A9 dressing is two-qutrit Clifford, not generic U(9)

Producer: `analysis/w33_pass11050_monomial_two_qutrit_latent_clifford.py`
Certificate: `data/w33_pass11050_monomial_two_qutrit_latent_clifford.json`

For the induced direction `L=<Z>`, label the nine cosets by `|b,c>` with `b,c in F3`.

The latent H27 generators are

`X_A: |b,c> -> |b+1,c>`

`C_A: |b,c> -> |b,c+1>`

`Z_A: |b,c> -> omega^theta |b,c+b>`.

Thus `X_A` and `C_A` are qutrit translations and `Z_A` is a qutrit SUM shear plus a global `mu3` phase.

The shear normalizes the two-qutrit Weyl group:

`X_b -> X_b X_c`, `X_c -> X_c`, `Z_b -> Z_b`, `Z_c -> Z_b^-1 Z_c`.

Therefore the missing symmetry-changing latent action can be synthesized inside the two-qutrit Clifford group. A generic 9-mode unitary is not required.
