# Pass 11360 — the Jarlskog spectrum: every J₆ value lies in Q(cos 2π/9); one magic gate gives exactly 9/8

Producer: `analysis/w33_pass11360_jarlskog_arithmetic.py`
Certificate: `data/w33_pass11360_jarlskog_arithmetic.json`
Regression: `tests/test_w33_pass11355_11360.py`

**Setup.**
* J₆ (Pass 11355) is recomputed at 60 digits. The 216 Cliffords are rebuilt exactly (Pass 11237) and T is exact.
* It is evaluated on every violating one-qutrit word with one and two cubic gates (coset-reduced, Pass 11312).
* Each distinct value is identified as a + b·c + e·c² with c = 2cos(2π/9), by PSLQ. A value is accepted only if the
  relation holds to 10⁻⁵⁰.

**Results: the 9 distinct values.**

| gates | J₆ | exact form |
|---|---|---|
| 1 | 1.125 | **9/8** |
| 2 | 1.125 | 9/8 |
| 2 | 1.875 | **15/8** |
| 2 | 0.36832867… | 119/216 − (19/108)c + (1/27)c² |
| 2 | 0.49429471… | 227/216 − (1/27)c − (23/108)c² |
| 2 | 0.51734347… | 311/216 − (4/27)c − (8/27)c² |
| 2 | 0.67169740… | 119/216 − (4/27)c + (4/27)c² |
| 2 | 1.01237662… | 59/216 + (23/108)c + (19/108)c² |
| 2 | 1.35262580… | 119/216 + (8/27)c + (4/27)c² |

**Reading.**
* The substrate Jarlskog invariant takes values in the **same cubic field** as the time-reversal fidelity levels of
  Pass 11253.
* All denominators divide 216, the order of the projective Clifford group.
* A single magic gate breaks time reversal with the **rational** strength J₆ = 9/8.
