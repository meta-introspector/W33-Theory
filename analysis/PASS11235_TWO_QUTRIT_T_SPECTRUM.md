# Pass 11235 — random two-qutrit Clifford+cubic dynamics: T-violation becomes generic, and 2π/9 is the mildest violation only at low depth

Producer: `analysis/w33_pass11235_two_qutrit_t_spectrum.py`
Certificate: `data/w33_pass11235_two_qutrit_t_spectrum.json`
Regression: `tests/test_w33_pass11234_11235.py`

**Method.**
* Random words W = C_d T_{q_d} ⋯ C₁ T_{q₁} C₀, with uniformly random two-qutrit Cliffords C_i and the cubic phase T on a
  random qutrit, for d = 1, …, 4 cubic gates; 60 words per depth.
* For each word, the exact best time-reversal fidelity is computed over all 51 840 × 81 anti-unitary two-qutrit
  Cliffords. It is evaluated as max_{r,a} |tr(B_r Y_a)|/9, with B_r = C_r U* C_r† and Y_a = P_a† U P_a. That is about
  100× faster than Pass 11227's evaluator.
* Controls: (T⊗T)·SUM gives 0.8440296287, and SUM gives 1.

**Result.**

| cubic gates | T-violating | mildest violation (max F_T) | strongest (min F_T) |
|---|---|---|---|
| 1 | 5/60 | 0.8440296 | 0.8440296 |
| 2 | 39/60 | 0.8440296 | 0.7123860 |
| 3 | 56/60 | **0.9392626** | 0.5948590 |
| 4 | 59/60 | 0.9392626 | 0.5669408 |

* **T-violation becomes generic.** On two qutrits the violating fraction is 8%, 65%, 93% and 98% for 1–4 cubic gates.
  That is faster than on one qutrit (Pass 11213).
* **The fidelity spectrum is discrete at each depth**, and widens with depth.
* **The 2π/9 value (1 + 2cos 2π/9)/3 is the mildest violation only with at most two cubic gates.**
  * Every violator with 1 or 2 cubic gates has F_T ≤ 0.8440296.
  * With 3 or more, milder violations occur: F_T = 0.9392626 in 3 words.
* So "the quantum of T-violation" is precisely the common fidelity of every minimal violator (Passes 11227 and 11228),
  and the bound on violations with at most two cubic gates. It is **not** a lower bound on T-violation in general.

**Correction this records.** The Pass 11227–11228 paper paragraph said every violating tick *found* had fidelity
(1 + 2cos 2π/9)/3. That was already contradicted by deeper one-qutrit words in Pass 11228's own data (minimum 0.7258).
The paper sentence is corrected in this commit to the scoped statement above.
