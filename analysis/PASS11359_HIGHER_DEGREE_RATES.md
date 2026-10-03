# Pass 11359 — higher-degree contraction rates are exact rationals: 3/8, 11/24, 1/2, growing with degree

Producer: `analysis/w33_pass11359_higher_degree_rates.py`
Certificate: `data/w33_pass11359_higher_degree_rates.json`
Regression: `tests/test_w33_pass11355_11360.py`

**Setup.** This continues Pass 11354, now matrix-free:
* π(C) acts leg by leg on tensors;
* the Clifford-invariant subspace is the range of the Clifford average applied to random vectors, with its dimension
  known exactly from characters (the computed rank matches every time);
* π(T) is diagonal;
* ρ is the largest modulus below 1 of B† π(T) B.

| (p, q) | dim | Clifford invariants | SU(3) invariants | ρ (second modulus) |
|---|---|---|---|---|
| (3,3) | 729 | 7 | 6 | **3/8** |
| (4,4) | 6561 | 40 | 23 | **3/8** |
| (6,3) | 19 683 | 105 | 47 | **3/8** |
| (9,0) | 19 683 | 91 | 42 | **1/2** (then 3/8) |
| (5,5) | 59 049 | 301 | 103 | **11/24** (then 3/8) |

**Reading.**
* The contraction rates of the Clifford-averaged cubic gate are exact simple rationals.
* New, slower modes appear as the degree grows: 3/8 → 11/24 → 1/2.
* That is consistent with the reversible fraction's empirical 0.72 per gate (Pass 11312) being governed by high-degree
  representations, as Pass 11354 argued.
* Whether the supremum over all degrees equals the observed rate, and what it is in closed form, are open. Degree
  (6,6) needs 3¹² = 531 441-dimensional tensors with 2542 invariants, beyond this matrix-free approach on one machine.
