# Pass 11149 — where the magic of spooky action sits, and how much it needs: no threshold

Producer: `analysis/w33_pass11149_magic_for_spookiness.py`
Scan (minimum-mana curve): `analysis/w33_pass11149_scan_min_mana.py`
Frozen: `data/w33_pass11149_min_mana_curve.json`
Regression: `tests/test_w33_pass11149_magic_for_spookiness.py`

Magic is measured by mana, M = log Σ|W| of the discrete Wigner function on F₃^{2n}; stabilizer states have M = 0. The test
is the qutrit CGLMP inequality (local bound 2).

| resources | best CGLMP value | magic |
|---|---|---|
| stabilizer states + Pauli measurements (W(3,3)) | 2 | none (Pass 11144) |
| **any state + Pauli measurements only** | **2√2·cos(⅓·arccos(1/(2√2))) = 2.6016791** | state mana 0.392 |
| stabilizer \|Ω⟩ + magic (CGLMP) measurements | 2.8729 | measurement vectors, mana up to 0.437 |
| optimal state (Acín et al.) + CGLMP measurements | 2.9149 | state mana 0.173 |

* **A closed form.** The optimal Pauli-measurement Bell operator's spectrum splits into the roots of x³ − 6x − 2 and of
  x³ − 3x + 1 (2cos(π/9), 2cos(2π/9), 2cos(4π/9)). The maximum is the largest root of **x³ − 6x − 2 = 0**.
* **Magic can sit in either place.** The state alone suffices, and so do the measurements alone.
* **There is no threshold.** The minimum mana needed to reach 2 + ε rises linearly from zero, at about 0.63–0.66 per unit
  of violation:

| violation I − 2 | minimum mana |
|---|---|
| 0.019 | 0.0128 |
| 0.099 | 0.064 |
| 0.299 | 0.183 |
| 0.598 | 0.378 |

This is the state-side counterpart of Pass 11144's quadratic onset on the measurement side. W(3,3) sits exactly on the
boundary of classical correlations, and any magic at all crosses it.
