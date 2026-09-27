# Pass 11044 — latent dimension nine is minimal, and it explains 27 + 54

Producer: `analysis/w33_pass11044_minimal_latent_dimension_and_54_budget.py`
Certificate: `data/w33_pass11044_minimal_latent_dimension_and_54_budget.json`

Write a general latent module as `N` one-dimensional characters plus `p V` plus `q Vbar`.

Requiring `V tensor A = Reg(H27)` forces exactly

`N=3, p=1, q=1`.

So the latent dimension is uniquely minimal at

`3 + 3 + 3 = 9`.

The old carrier is `9V`; the regular carrier contains only `3V`. Their maximal common H27 submodule therefore has dimension 9. Tensoring by the external regular `C3` gives 27 compatible dimensions out of 81, hence exactly 54 dimensions must change representation type.

That reproduces the old equivariant rank ceiling and the minimal compiler budget from a single representation-ring theorem.

Choosing the three one-dimensional characters as an affine line in dual `F3^2` gives 12 natural Hesse-compatible gauges; the theorem itself does not require that geometric choice.
