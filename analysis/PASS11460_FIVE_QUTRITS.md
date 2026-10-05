# Pass 11460 — five qutrits: the one-gate bad-class fraction is 0.138 ± 0.004, not 1/8

Producer: `analysis/w33_pass11460_five_qutrits.py`
Certificate: `data/w33_pass11460_five_qutrits.json`
Regression: `tests/test_w33_pass11456_11460.py`

**Method.**
* Pass 11437's estimator: P_n = Σ_cells (exact share) × (sampled bad fraction of the cell).
* Here at n = 5, with Pass 11421's sparse decider: 243-dimensional states and 59 049 frames per class.
* 4000 classes per cell type. The n = 3 control in Pass 11437 reproduced the exact 437/3276.
* The shares are exact: collinear 2460/7381, non-collinear 19683/29524.

| cell | n = 3 (exact) | n = 4 (Pass 11437) | **n = 5** |
|---|---|---|---|
| same line, other | 0.225 | 0.140 | **194/1327 = 0.146** |
| collinear, different lines | 0.111 | 0.125 | **380/2672 = 0.142** |
| all collinear cells | 0.1505 | 0.130 | **0.1435** |
| non-collinear | 0.1235 | 0.127 | **0.1356** |
| **P_n** | **0.13339** | **0.128 ± 0.005** | **0.1382 ± 0.0041** |

Three classes were undecided.

**The sequence:** P₁ = P₂ = 1/8, P₃ = 437/3276 = 0.13339, P₄ ≈ 0.125–0.128, **P₅ = 0.1382 ± 0.0041**.

**Reading.**
* **The conjecture is refuted.** P₅ lies **3.3σ above 1/8**, so the conjecture "P_n → 1/8" recorded in Pass 11437
  does not survive. (A correction is added to that note.)
* **No settling at 1/8.** P₅ is consistent with the exact n = 3 value (1.2σ) and 1.6σ above the n = 4 cell estimate.
  The bad-class fraction does not settle at 1/8; it hovers between about 0.12 and 0.14, with no clear trend in n.
* **The cells equalise.** At n = 5 every cell sits near 0.14. The cell structure that made n = 3 special (77/342 in
  the same-line cell) has washed out.
* **Open.** The large-n behaviour, i.e. whether P_n converges and to what, is open. Six-qutrit sampling would need a
  729-dimensional decider (531 441 frames per class).

**Fixed-axis law at n = 5.** 40 targeted fixed-axis classes were checked with the decider capped at 3⁹:
* **24 satisfy F′ with exact equality** (violating frames = {a_{x₁} ≠ 0});
* 16 are beyond the cap;
* none fails.

The law holds on every decided five-qutrit class tested.
