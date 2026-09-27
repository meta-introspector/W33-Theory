# Pass 11049 — the twelve induced sectors are spectral projectors with Hesse incidence

Producer: `analysis/w33_pass11049_right_regular_hesse_projectors.py`
Certificate: `data/w33_pass11049_right_regular_hesse_projectors.json`

In the regular H27 carrier, let `R_g` be right translation by a generator of one noncentral C3 direction. The three spectral projectors

`P_(d,theta) = (I + omega^-theta R_g + omega^-2theta R_g^2)/3`

have rank 9 and commute with the left regular H27 action.

Each of the four directions gives an orthogonal `9+9+9` decomposition of `C[H27]`. Projectors from distinct directions intersect in exactly one dimension. That common ray is precisely the unique one-dimensional H27 character at the intersection of the corresponding dual affine lines.

The cross-projector principal-angle spectrum is rigid:

`PQP | Im(P): 1^1, (1/3)^6, 0^2`.

So this is not an isoclinic fusion construction: two nonparallel compiler sectors share one exact ray, have six qutrit-unbiased directions, and two orthogonal directions.
