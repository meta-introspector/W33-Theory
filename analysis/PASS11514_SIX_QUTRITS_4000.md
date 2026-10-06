# Pass 11514 — six qutrits, one magic gate, pushed toward 4000 classes

Producer: `analysis/w33_pass11514_six_qutrits_4000.py` (arguments: classes per kind, workers; `--summarise-only`)
Certificate: `data/w33_pass11514_six_qutrits_4000.json` (per-class checkpoint `data/w33_pass11514_checkpoint.json`)
Regression: `tests/test_w33_pass11511_11515.py`

**Method.** This continues the streamed census of Passes 11489 and 11501: same seeds, same interleaving, every completed
class used. The first 2500 classes are 11501's.

## Result

4000 classes: 2000 collinear and 2000 non-collinear. Three were beyond the decider and are excluded.

| cell | bad | good | bad fraction |
|---|---|---|---|
| collinear, different lines | 183 | 1147 | 0.138 |
| same line, other | 95 | 573 | 0.142 |
| collinear (both) | 278 | 1720 | 0.139 |
| non-collinear | 254 | 1745 | 0.127 |

* **P₆ = 0.131 ± 0.006.**

| n | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| bad-class fraction | 1/8 | 1/8 | 437/3276 = 0.1334 (exact) | 0.128 ± 0.005 | 0.138 ± 0.004 | **0.131 ± 0.006** |

**Reading.**
* Still inconclusive. P₆ is 1.1σ above 1/8 and 1.0σ below P₅.
* The sequence from n = 3 on (0.133, 0.128, 0.138, 0.131) is consistent with a constant near 0.132.
* It is also consistent with fluctuations around 1/8 plus a small positive correction.
* With the n = 4 to 6 error bars, no trend can be claimed.
* The magic-axis and good-rule cells (Passes 11498, 11499, 11503) have share 6/(3^{2n} − 1) and do not move this number.
