# Pass 11489 — six qutrits, one magic gate: the bad-class fraction at n = 6

Producer: `analysis/w33_pass11489_six_qutrits.py` (arguments: classes per kind, workers; rerunning resumes from the
checkpoint)
Certificate: `data/w33_pass11489_six_qutrits.json` (checkpoint `data/w33_pass11489_six_qutrits_checkpoint.json`)
Regression: `tests/test_w33_pass11486_11492.py`

**Method.**
* The estimator is the one from Passes 11437 and 11460: exact cell shares times sampled within-cell fractions.
* The decider is Pass 11421's sparse decider at n = 6: 729-dimensional states and 531 441 frames per class.
* Collinear and non-collinear classes are interleaved, so the sample stays balanced whenever the run is stopped.
* The run streams and checkpoints every 100 classes.

## Result

1000 classes: 500 collinear and 500 non-collinear, all decided.

| cell | bad | good | bad fraction |
|---|---|---|---|
| collinear, different lines | 42 | 290 | 0.127 |
| same line, other | 29 | 139 | 0.173 |
| collinear (both) | 71 | 429 | 0.142 |
| non-collinear | 57 | 443 | 0.114 |

* Cell shares are exact: collinear 22143/66430, non-collinear 177147/265720.
* **P₆ = 0.123 ± 0.011.**

## The sequence

| n | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| bad-class fraction | 1/8 | 1/8 | 437/3276 = 0.1334 (exact) | 0.128 ± 0.005 | 0.138 ± 0.004 | **0.123 ± 0.011** |

**Reading.**
* **Inconclusive.** P₆ agrees with 1/8 (z = −0.15) and with P₅ (z = −1.3).
* The sample is about a quarter of the n = 5 one, so its error bar is 2.7 times larger. Nothing about the n-dependence is
  claimed.
* The n = 5 excess over 1/8 (3.3σ, Pass 11460) is neither confirmed nor refuted.
* **The checkpoint is the asset.** Rerunning the producer with a larger first argument continues from the 1000 classes,
  with no recomputation. About 4000 classes would reach the n = 5 precision.

## Process note

* A first, unstreamed attempt ran Pool.map over 5000 classes on 6 workers.
* It produced no output for 2.5 hours, and it starved every other computation of commit memory. It was stopped and
  nothing from it is used.
* The rerun streams, checkpoints, interleaves the two kinds and caps the workers.
* Stopping it early is a sample-size choice, not a selection. Every completed class is used.
