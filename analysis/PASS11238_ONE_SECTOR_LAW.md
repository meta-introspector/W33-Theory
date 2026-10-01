# Pass 11238 — the "both sectors" law is false as stated; what is true, exactly, on all short words

Producer: `analysis/w33_pass11238_one_sector_law.py`
Certificate: `data/w33_pass11238_one_sector_law.json`
Regression: `tests/test_w33_pass11236_11239.py`

**What was claimed (Pass 11227).** One cubic phase with an entangler, or two cubic phases on the same qutrit, is always
inverted by a substrate time reversal, so T-violation needs the cubic phase on both qutrits. That came from ten circuits
tested.

**Exhaustive test.**
* Every distinct word (as a matrix up to phase) containing at least one cubic gate, of length ≤ 5 over the target-only
  and control-only gate sets, and of length ≤ 4 over the both-sector set.
* Each word gets the exact best time-reversal fidelity over all 51 840 × 81 anti-unitary two-qutrit Cliffords.

| cubic phases on | gate set | distinct words | T-violating | fidelity levels | first violation |
|---|---|---:|---:|---|---|
| control only | {T⊗I, SUM, SUM†} | 13 | **0** | — | never |
| target only | {I⊗T, SUM, SUM†} | 66 | **12** | 0.712386014 | **3 cubic gates**, e.g. (I⊗T)·SUM·(I⊗T²) |
| both | {T⊗I, I⊗T, SUM, SUM†} | 70 | 20 | 0.8440296 (16), 0.712386 (4) | 2 cubic gates, (T⊗T)·SUM |

**What is true.**
* **Cubic phases only on the control qutrit never break T.** A diagonal gate on the control commutes with SUM, so every
  such word reduces to (T^k ⊗ I)·SUM^m. All 13 classes are reversible.
* **Cubic phases only on the target can break T**, but not before three cubic gates. The minimal examples, such as
  (I⊗T)·SUM·(I⊗T²) and (I⊗T²)·SUM·(I⊗T), all sit at the stronger violation level F = 0.712386.
* **Phases on both qutrits break T already with two gates**, at the 2π/9 level F = (1 + 2cos 2π/9)/3.

**Correction.** "T-violation needs the phase in both sectors" (Pass 11227 and its paper paragraph) is false as stated.
The correct statement is: **the cheapest T-violation needs both sectors (two cubic gates); a single sector needs three,
and the control sector alone never violates.** Pass 11227 tested (I⊗T)·SUM·(I⊗T), which is reversible, but not
(I⊗T)·SUM·(I⊗T²). The paper and ledger are corrected in this commit.
