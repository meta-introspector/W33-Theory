# Pass 11147 — one qudit across time beats two qudits in space for every d ≥ 3, and the gap grows

Producer: `analysis/w33_pass11147_temporal_cglmp_by_dimension.py`
Scan: `analysis/w33_pass11147_scan_cglmp_by_dimension.py`
Frozen: `data/w33_pass11147_cglmp_optima_d2_5.json` (optimal parameters included)
Regression: `tests/test_w33_pass11147_temporal_cglmp_by_dimension.py`

Setup: the CGLMP inequality I_d (two settings, d outcomes, local bound 2, algebraic maximum 4).
* **Temporal:** a pure start and projective Lüders measurements at two times.
* **Spatial:** a pure bipartite state and projective measurements.

| d | spatial optimum | temporal optimum | gap |
|---|---|---|---|
| 2 | 2.8284271 (2√2) | 2.8284271 | 0 |
| 3 | 2.9148542 (1 + √(11/3)) | 3.1628065 | 0.2480 |
| 4 | 2.9726983 | 3.2830617 | 0.3104 |
| 5 | 3.0157105 | 3.3773031 | 0.3616 |

**Validation.**
* The spatial optima reproduce the known values. For d = 3, an integer-relation search recovers **1 + √(11/3)** exactly.
* d = 2 reproduces Fritz's theorem that temporal and spatial coincide for qubits.

**Findings.**
* **Qubits are the only coincidence.** For every d ≥ 3 a single system measured twice exceeds the best pair, and the gap
  grows monotonically with d. This matches Budroni–Emary's finding that temporal quantum values grow with dimension.
* **No closed form.** No low-height algebraic form was found for the temporal optima; only spurious high-coefficient
  fits appear at 13-digit precision. They are reported numerically.
* **Why time wins.** A measurement can hand information forward in time, and more levels carry more of it. Space
  forbids this.
