# Pass 11234 — independent brute-force confirmation of the qubit arrow law and protection numbers

Producer: `analysis/w33_pass11234_qubit_bruteforce_check.py`
Certificate: `data/w33_pass11234_qubit_bruteforce_check.json`
Regression: `tests/test_w33_pass11234_11235.py`

**What is checked.**
* The parallel session's Pass 11217 proved the qubit arrow law A(S) = n − c(S).
* Its Pass 11220 computed the mean number of protected qubits from GAP's conjugacy classes and the Jordan closed form.
* Its own brute-force check at n = 6 had to be stopped.
* Following the repo rule that a cross-check must be independent, this pass uses **neither the closed form nor the class
  list**:
  * every element of Sp(4,2) and Sp(6,2) is enumerated, by breadth-first search from transvections;
  * c(S) is computed directly, as the largest mutually orthogonal family of S-invariant nondegenerate planes of F₂^{2n};
  * A(S) is computed directly, as the minimum over **all** qubit factorisations of Σ_q rank S_F[others, q]. There are
    10 factorisations for two qubits and 1120 for three.

**Result.**

| n | elements | distribution of c | E_n[c] | Pass 11220 | A = n − c |
|---|---:|---|---|---|---|
| 2 | 720 | c=0: 459, c=2: 261 | **29/40** | 29/40 ✓ | all 720 ✓ |
| 3 | 1 451 520 | c=0: 686 160, c=1: 663 579, c=3: 101 781 | **5981/8960** | 5981/8960 ✓ | 1500 random elements, all 1120 splits each ✓ |

* Both protection means agree exactly with Pass 11220.
* c = n − 1 never occurs, as the law requires.
* The qubit arrow law holds, with A taken directly as a minimum over all splits, on every element of Sp(4,2) and on
  every sampled element of Sp(6,2).

**Scope.** Exhaustive in c for n ≤ 3. The direct A is exhaustive for n = 2 and sampled (1500 of 1 451 520) for n = 3.
