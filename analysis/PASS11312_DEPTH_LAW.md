# Pass 11312 — one-qutrit T-violation versus the number of cubic gates: exact through four, then geometric

Producer: `analysis/w33_pass11312_depth_law.py`
Certificate: `data/w33_pass11312_depth_law.json`
Regression: `tests/test_w33_pass11308_11312.py`

**Reduction.**
* The nine diagonal Cliffords commute with T.
* So a uniform word C_k T ⋯ C_1 T has the same distribution as R_k T ⋯ R_2 T·C_1 T, with each R_i one of 24 coset
  representatives of Cl/Diag and C_1 ranging over all 216.
* That leaves 216·24^{k−1} equally weighted words.
* The reduction reproduces the exhaustive values of Pass 11266 for k = 1, 2, 3.

| cubic gates k | probability of violating substrate time reversal | reversible fraction |
|---|---|---|
| 1 | **1/12** (exact) | 0.917 |
| 2 | **13/32** (exact) | 0.594 |
| 3 | **421/768** (exact) | 0.452 |
| 4 | **1425/2048** = 0.69580 (exact, 2 985 984 words) | 0.304 |
| 5 | 0.7764 ± 0.0016 (66 000 samples) | 0.224 |
| 6 | 0.8425 ± 0.0014 | 0.157 |
| 8 | 0.9212 ± 0.0010 | 0.079 |
| 12 | 0.9788 ± 0.0006 | 0.021 |

**Reading.**
* From k = 4 the reversible fraction decays geometrically, by a factor of about **0.72 per cubic gate**.
* **Corrected (Codex's audit; Pass 11332).** An earlier version said that every long enough magic circuit breaks the
  arrow. That is false: T⁹ = I, so reversible circuits exist at every depth 9m. What is true, and proved in Pass 11332,
  is the probabilistic statement: under this random walk P(reversible) → 0 (Kawada–Itô equidistribution plus the
  portmanteau theorem on the closed, Haar-null reversible set). The 0.72-per-gate rate is an empirical fit, not a
  theorem.
* The exact rate, and a closed form for the sequence 1/12, 13/32, 421/768, 1425/2048, are OPEN. The denominators are
  12, 2⁵, 2⁸·3 and 2¹¹.
