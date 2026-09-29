# Pass 11159 — an explicit state achieving 4/√15, with the value derived

Producer: `analysis/w33_pass11159_explicit_polygamy_state.py`
Regression: `tests/test_w33_pass11159_explicit_polygamy_state.py`

In the eigenbases of its reduced states, the numerical optimum is W-like. That suggests the family

    |ψ(x)⟩ = √x |000⟩ + √((1−x)/4) Σ_{a=1,2} |a⟩(|0a⟩ + |a0⟩).

* **Closed form (exact, sympy):** the partial-transpose spectrum gives

      N_AB(x) = N_AC(x) = √((1 + 15x)(1 − x)) / 4.

* **Maximum:** 14 − 30x = 0 gives x = 7/15, where (1 + 15x)(1 − x) = 64/15. So **N_AB = N_AC = 2/√15** and
  N_AB + N_AC = **4/√15**, now derived rather than fitted.
* **The explicit optimum:** |ψ*⟩ = √(7/15)|000⟩ + √(2/15) Σ_{a=1,2} |a⟩(|0a⟩ + |a0⟩). It reproduces every invariant of
  the numerical optimum: spectra (7,4,4)/15 and (11,2,2)/15, N_BC = 1/15.

**Status:** 4/√15 is exactly the maximum of this family and coincides with the numerical global optimum from 46 random
restarts. A proof that no state outside the family does better is still open.
