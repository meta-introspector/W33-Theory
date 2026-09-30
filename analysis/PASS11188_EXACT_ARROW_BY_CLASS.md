# Pass 11188 — the intrinsic arrow of time, exactly, on every Clifford class

Producer: `analysis/w33_pass11188_exact_arrow_by_class.py` (GAP: `analysis/gap/w33_pass11188_class_reps.g`)
Certificate: `data/w33_pass11188_exact_arrow_by_class.json` (GAP output `data/w33_pass11188_gap_class_reps.txt`)
Regression: `tests/test_w33_pass11188_exact_arrow_by_class.py`

**What was open.** Pass 11183 defined the intrinsic arrow A(S): the minimum, over all stabilizer tensor
factorisations, of the information a tick exports (Σ_q rank S_F[others, q], in trits). The result was exact for two
qutrits (A ∈ {0, 2}) but only a 240-element sample for three qutrits (A ∈ {0, 2, 3}).

**Method.**
* A is a class function: S_{gF}(gSg⁻¹) = S_F(S), and −S changes no rank.
* GAP prints one matrix per conjugacy class of PSp(4,3) (20 classes) and PSp(6,3) (74 classes). The matrices are
  written in an adapted symplectic basis (ω(uᵢ, vᵢ) = 1).
* Each representative is checked three ways:
  * symplectic for our form;
  * of GAP's projective order;
  * local in exactly GAP's number of fixed factorisations (the Pass 11180 cross-check, all 94 classes).
* Each representative is then scanned over all 45 or 110 565 factorisations.

**Result (exact, with class sizes).**

| | A = 0 | A = 2 | A = 3 |
|---|---|---|---|
| two qutrits, PSp(4,3) | 133/360 | 227/360 | — |
| three qutrits, PSp(6,3) | 55241/884520 | 581/1080 | 8836/22113 |

The three-qutrit counts are 286 369 344, 2 466 749 376 and 1 832 232 960 elements of 4 585 351 680.

* **A is never 1 and never exceeds the number of qutrits (A ≤ n for n = 2, 3).** This is now a theorem over all
  classes, not a sample.
* **Every element whose order is divisible by 7 or 13 has the maximal arrow 3.** These are the Singer-type ticks
  (3³ ± 1 = 26, 28), which no split can make leak-free.
* A = 3 also holds for the order-3 class that permutes the three qutrits cyclically. In characteristic 3 the 3-cycle
  is unipotent: it is local in 27 splits, yet in none of them does it export less than three trits. SWAP-type ticks
  have A = 0: SWAP is (I) ⊕ (−I) in the symmetric/antisymmetric split.
* A = 2 holds for every class of order 5, 10, 15, 20 or 30 (the prime 5 of 3² + 1), and for some classes of order
  8, 9, 12, 18, 24 or 36.

**Scope.** The arrow is a property of Clifford ticks on the stabilizer factorisations of two and three qutrits. The
minimisation is over stabilizer (symplectic) splits only.

**Prior art.** The idea of choosing subsystems by minimising what the dynamics does across them is Zanardi, Dallas,
Andreadakis & Lloyd, *Operational quantum mereology and minimal scrambling*, Quantum **8**, 1406 (2024),
arXiv:2212.14340, where scrambling is minimised over all tensor product structures. The Clifford/finite-field version
and the exact class distribution are not in that paper.
