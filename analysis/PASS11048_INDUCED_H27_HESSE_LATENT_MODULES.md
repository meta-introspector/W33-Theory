# Pass 11048 — the minimal latent module is induced from a noncentral C3

Producer: `analysis/w33_pass11048_induced_h27_hesse_latent_modules.py`
Certificate: `data/w33_pass11048_induced_h27_hesse_latent_modules.json`

The previously discovered minimal latent module

`A9 = 3 chi + V_omega + V_omega2`

has a canonical representation-theoretic meaning.

Choose one of the four noncentral projective directions `L < H27`, with `L ~= C3`, and one character `theta` of `L`. Frobenius reciprocity gives

`A9(d,theta) = Ind_L^H27(theta)`.

The three one-dimensional H27 characters that occur are exactly one affine line in the dual `F3^2`; both 3D Schrodinger irreps occur once because each restricts to the regular representation of `L`.

Thus the 12 minimal latent choices are exactly

`4 projective directions x 3 theta sectors = 12 Hesse affine lines`.

For any fixed direction, the three theta sectors direct-sum to `Reg(H27)`.
