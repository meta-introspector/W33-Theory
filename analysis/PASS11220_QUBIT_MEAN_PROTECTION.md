# Pass 11220: how many qubits a random Clifford tick protects

Producer: `analysis/w33_pass11220_qubit_mean_protection.py`
GAP classes of Sp(12,2): `analysis/gap/w33_pass11220_class_reps_q2_n6.g` → `data/w33_pass11220_gap_class_reps_q2_n6.txt`
Certificate: `data/w33_pass11220_qubit_mean_protection.json`
Regression: `tests/test_w33_pass11220_qubit_mean_protection.py`

**Background.**
* Pass 11217 proved the qubit arrow law A = n − c and gave the closed form c = m₁(1)/2 + χ₂·m₂(1) + m₁(x²+x+1). That
  form is exact on all 320 classes with n ≤ 5.
* Master's Pass 11215 found the qutrit mean E_n[c] for every n, with limit 0.72553535…
* The qubit analogue of master's cycle-index method needs Hesselink's index. This pass instead sums exactly over GAP's
  conjugacy classes, now including Sp(12,2): 477 classes, |Sp(12,2)| = 208 114 637 736 580 743 168 000.

**Results (exact).**

| n | classes | E_n[c] | decimal |
|---|---:|---|---|
| 2 | 11 | 29/40 | 0.725 |
| 3 | 30 | 5981/8960 | 0.667522 |
| 4 | 81 | 12945139/19496960 | 0.663957 |
| 5 | 198 | 6791923649683/10212039720960 | 0.665090 |
| 6 | 477 | 180847576047477961/271885345530839040 | **0.665161** |

* The successive differences are −0.0575, −0.00357, +0.00113 and +0.0000714. The last ratio is 0.063 ≈ 1/16, so
  geometric extrapolation gives **lim E_n[c] ≈ 0.665166**. This is an estimate, not a proof.
* A random qubit tick therefore protects about 0.665 qubits. A random qutrit tick protects 0.7255 qutrits.
* At n = 6:
  * the maximal arrow (c = 0) has probability 0.4927;
  * an arrow-free tick (c = 6) has probability 7.7·10⁻⁶;
  * c = 5 never occurs, as A ≠ 1 requires.
* **Cross-check.** The n = 6 closed form is compared with brute-force c over all 1 397 760 nondegenerate planes of F₂¹²
  (`--brute`). The outcome is recorded in the certificate when the run completes.
* **Open.** An exact qubit limit needs Fulman's cycle index refined by Hesselink's index (Fulman–Neumann–Praeger treat
  even characteristic).
