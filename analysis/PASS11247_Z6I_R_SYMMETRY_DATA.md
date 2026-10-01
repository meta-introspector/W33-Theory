# Pass 11247: the Z6-I R-symmetries are not in the data (blocked, stated exactly)

Producer: `analysis/w33_pass11247_z6i_r_symmetry_data.py`
Certificate: `data/w33_pass11247_z6i_r_symmetry_data.json`

**Why it matters.** Pass 11241 found that in Z6-I, the W33 vacuum's class, no vacuum symmetry encoded in the frozen
charges forbids μ while allowing the top Yukawa. The search was exhaustive over every FI ray and every field choice.
The only loophole left is an R-symmetry.

**What the data contain.** Of the 87 Z6-I models in Pass 10967's discrete-charge file, **none** carries R-charges;
only the point-group Z₆ is recorded. All 128 Z6-II models carry Z₆^R × Z₃^R × Z₂^R.

**What is needed.** For every field: the shifted H-momenta q_sh^i and the oscillator numbers N^i, N̄^i in each
twisted sector. Then R^i = q_sh^i − N^i + N̄^i, modulo the plane-rotation orders. For the Z6-I twist
v = (1/6, 1/6, −1/3) these orders are (6, 6, 3). The values come from the orbifolder's raw spectrum files, which are
not in this repository.

**Status.** BLOCKED on data. The Z6-I μ statement stays scoped to "U(1)s, their discrete remnants and the point group".
