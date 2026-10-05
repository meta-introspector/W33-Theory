# Pass 11501 — six qutrits, one magic gate, extended toward the n = 5 precision

Producer: `analysis/w33_pass11501_six_qutrits_extended.py` (arguments: classes per kind, workers; `--summarise-only`)
Certificate: `data/w33_pass11501_six_qutrits_extended.json` (per-class checkpoint `data/w33_pass11501_checkpoint.json`)
Regression: `tests/test_w33_pass11498_11505.py`

**Method.**
* This continues Pass 11489's streamed and checkpointed census, with the same seeds and the same interleaving of
  collinear and non-collinear classes. The first 1000 classes are Pass 11489's.
* Pass 11489's committed certificate is left unchanged.
* Every completed class is used.

## Result

2500 classes: 1250 collinear and 1250 non-collinear. One non-collinear class was beyond the decider and is excluded.

| cell | bad | good | bad fraction |
|---|---|---|---|
| collinear, different lines | 116 | 715 | 0.140 |
| same line, other | 66 | 353 | 0.158 |
| collinear (both) | 182 | 1068 | 0.146 |
| non-collinear | 149 | 1100 | 0.119 |

* **P₆ = 0.128 ± 0.007.**
* The first 1000 classes (Pass 11489) gave 0.123 ± 0.011.

| n | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| bad-class fraction | 1/8 | 1/8 | 437/3276 (exact) | 0.128 ± 0.005 | 0.138 ± 0.004 | **0.128 ± 0.007** |

**Reading.**
* **Still inconclusive.** P₆ agrees with 1/8 (z = +0.44) and with P₅ = 0.138 ± 0.004 (z = −1.26).
* The error bar is now 1.7 times the n = 5 one.
* No n-dependence beyond n = 3's exact 437/3276 is established.
* **Why the run stopped here.** It was planned for 4000 classes. It was stopped at 2500 so the round could be
  published, and every completed class is used. The checkpoint resumes it.
