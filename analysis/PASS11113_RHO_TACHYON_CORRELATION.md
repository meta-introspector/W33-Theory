# Pass 11113 — "tachyon-free at ρ ⇒ untwisted quarks" at scale: suggestive, not significant

Producer: `analysis/w33_pass11113_rho_tachyon_correlation.py`
Data: `data/w33_pass11113_quark_sectors_387.json` (up-quark sector of all 387 new models, Pass 11098 engine, full spectra
from both engines)
Certificate: `data/w33_pass11113_rho_tachyon_correlation.json`
Regression: `tests/test_w33_pass11113_rho_tachyon_correlation.py`

The table below covers all 491 models, the 104 of Pass 11095 plus the 387 of the Pass 11108 rescan:

|  | twisted up-quarks | untwisted up-quarks |
|---|---|---|
| tachyon-free with the Wilson-line tori at ρ | **0** | 5 |
| tachyonic at ρ | 152 | 334 |

Under independence the expected count of ρ-safe twisted models is 1.55. The Poisson probability of seeing 0 is 0.21,
and Fisher's exact one-sided p is **0.16**. **The correlation is not significant.** Pass 11108 drew its reading
("tachyon-freedom and a realistic quark sector exclude each other") from 0/31, a convenience sample, and that reading
is withdrawn.

What stands:
* no model found so far is both ρ-safe and has twisted quarks;
* the ρ-safe rate is low everywhere, about 1%;
* the lowest ρ level is universal: Δ ∈ {−1/18, −1/6} in every tachyonic model.

A theorem either way would need the ρ tachyon condition in closed form in terms of the Wilson lines. That is not
derived here.
