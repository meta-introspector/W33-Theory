# Pass 11459 — the reversible fraction decays at a constant ≈ 0.72 per gate to depth 30, and the ρ₍₈,₈₎ coincidence is back (with a correction of Pass 11422's sampled tail)

Producer: `analysis/w33_pass11459_deep_decay.py` (main run; `--nomerge`; `--sens`)
Certificate: `data/w33_pass11459_deep_decay.json`
Regression: `tests/test_w33_pass11456_11460.py`

## A threshold artefact, found and corrected

* **The cut matters.** The deciders of Passes 11422 and 11436 call a word reversible when its best Clifford overlap
  max_C|tr((CUC†)†Uᵀ)| exceeds 2.999.
  * Exactly reversible words reach **3 to rounding**, with a deficit below 10⁻¹⁰.
  * Clifford+T words are dense, so deep violating words come **arbitrarily close**.

| depth | count above 2.999 | count above 3 − 10⁻⁶ = 3 − 10⁻⁸ = 3 − 10⁻¹⁰ | false at 2.999 |
|---|---|---|---|
| 12 | 13 129 | 13 042 | 0.7% |
| 24 | 412 | 281 | **32%** |
| 30 | 185 | 39 | **79%** |

* No word has a deficit between 10⁻¹⁰ and 10⁻⁶, so the tight cut **3 − 10⁻⁷** separates cleanly.
* **What this invalidates.** A first run of this pass with the loose cut appeared to show the decay slowing: local
  rates 0.7285, 0.771 and 0.875. That was this artefact. It also invalidates Pass 11422's "upward drift" reading;
  corrections are added to Passes 11422 and 11436.
* **What stays valid.** The exact enumerations (k ≤ 7) stay valid: their largest violating overlap is 2.99861.

## Corrected measurement (10 million samples per depth, K-reduced walk of Pass 11436, cut 3 − 10⁻⁷)

| k | 12 | 18 | 24 | 30 |
|---|---|---|---|---|
| reversible fraction | 0.021095 | 0.0029618 | 4.2327×10⁻⁴ | 6.291×10⁻⁵ |
| local per-gate rate | — | **0.7209 ± 0.0007** | **0.7231 ± 0.0020** | **0.728 ± 0.005** |

* **No slowdown.** The rate is essentially constant at about 0.722 per gate from depth 12 to 30. The depth-12 value
  agrees with Pass 11422's independent sample (0.021059 ± 0.00007).
* **The coincidence returns.** ρ₍₈,₈₎ = 0.72276 (Pass 11372; the largest root of 2592x³ − 1332x² − 585x + 140)
  **lies within one standard error of the last two windows**. The first window is 2.4σ below and approaching it.
* **No mechanism.** No Fourier rate theorem governs a Haar-null event, and the largest irrep rate is 0.781 at (18,0),
  not 0.7228. So this is a **sharpened, unexplained coincidence**, not an identity.

## Mechanism test: excluding the identity coset (no T-merging)

* With the identity coset removed from the 8 per gate, adjacent T gates can never merge.
* Local rates: **0.6987 ± 0.0012, 0.7019 ± 0.0025, 0.6963 ± 0.0050** (depths 12→16→20→24). Also constant, and a
  little faster than the full walk.
* T-merging contributes a small slowing (0.70 → 0.72), not the dominant effect.

**Lesson.** A threshold that separates cleanly at shallow depth can fail at depth when the objects are dense. Here
the near-miss mass grew until it outweighed the signal. Every sampled count must be re-checked against the tightest
cut that still separates, at the deepest depth used.

**Name clash, not prior art.** "Deficit" here is the overlap shortfall 3 − max_C|tr((CUC†)†Uᵀ)|. The deficits in
`PASS11160_TEMPORAL_CGLMP_HIGH_D.md` (a temporal CGLMP gap), `w33_clock_magic_renewal.py` and
`w33_magic_resource_accounting.py` (a Kochen–Specker deficit) are different quantities.
