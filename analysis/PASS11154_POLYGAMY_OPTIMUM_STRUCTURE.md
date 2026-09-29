# Pass 11154 — the 4/√15 optimum: every invariant is a fifteenth

Producer: `analysis/w33_pass11154_polygamy_optimum_structure.py`
Scan: `analysis/w33_pass11154_scan_polygamy_optimum.py`
Frozen (optimal state + invariants): `data/w33_pass11154_polygamy_optimum_state.json`
Regression: `tests/test_w33_pass11154_polygamy_optimum_structure.py`

The best a qutrit can split its entanglement between two partners in space is N_AB + N_AC = 4/√15 (Pass 11148). The
state achieving it has, to about 10⁻⁸:

| invariant | value |
|---|---|
| spectrum of A | (7, 4, 4)/15 |
| spectrum of B and of C | (11, 2, 2)/15 (B ↔ C symmetric) |
| N_AB = N_AC | **2/√15** |
| N_BC | 1/15 |
| spectrum of ρ_AB^{T_B} | (1 − √15)/15 ×2, −2/15, 2/15 ×3, (1 + √15)/15 ×2, 7/15 |

The three negative partial-transpose eigenvalues sum to 2/√15.

**Contrast with time:** the same qutrit, persisting under the identity, has negativity **1** with its past and 1 with its
future (Pass 11148).

**Status:**
* The invariants are identified exactly.
* That 4/√15 is the global maximum is supported by 46 random restarts (11148 and 11154) but **not proved**.
* An analytic proof, for example from the rational spectra, is open.
