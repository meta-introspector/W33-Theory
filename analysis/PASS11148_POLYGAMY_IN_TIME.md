# Pass 11148 — entanglement across time is polygamous: a qutrit is maximally entangled with its past and its future at once

Producer: `analysis/w33_pass11148_polygamy_in_time.py`
Scan: `analysis/w33_pass11148_scan_polygamy.py`
Frozen: `data/w33_pass11148_polygamy.json`
Regression: `tests/test_w33_pass11148_polygamy_in_time.py`

| | pair (t₁,t₂) + pair (t₂,t₃) | all three pairs |
|---|---|---|
| **time**: one qutrit, identity evolution, maximally mixed start | **2** | **3** |
| **space**: best three-qutrit state (numerical optimum) | **1.0327956 = 4/√15** | 1.1394 |
| references | antisymmetric state: ⅓ per pair; GHZ: 0 per pair | |

* **Time.** Every pair of times has the same two-time pseudo-density operator, SWAP/3: the temporal Bell line of
  Pass 11143, with negativity 1. The Bell lines Δ₁₂, Δ₂₃ and Δ₁₃ hold simultaneously.
* **Space.** Negativity is convex, so pure states give the spatial optimum. Its value 1.0327955589886 equals 4/√15 to
  10⁻¹⁵: an integer-relation search gives (x − 1)(15x² − 16) = 0. This is an identification, not a proof.
* **Reading.** In space a qutrit's maximal entanglement cannot be shared (monogamy). Across time it can, with its entire
  history, because it is persistence, not a bond divided between partners. This is also why the "spooky" spatial
  correlations are scarce while temporal ones are free. Temporal non-monogamy of pseudo-density operators is known in
  general; the qutrit numbers and the W(3,3) reading are specific here.
