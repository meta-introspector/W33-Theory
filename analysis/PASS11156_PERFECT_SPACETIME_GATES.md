# Pass 11156 — perfect space–time gates: a complementarity law, and the tetracode as the one perfect gate

Producer: `analysis/w33_pass11156_perfect_spacetime_gates.py` (exact over F₃, all 51 840 elements of Sp(4,3))
Regression: `tests/test_w33_pass11156_perfect_spacetime_gates.py`

**Setup.** A two-qutrit Clifford gate is S = [[S_AA, S_AB], [S_BA, S_BB]] ∈ Sp(4,3). Its Choi state lives on A_in, B_in,
A_out, B_out, and its three 2|2 cuts measure three kinds of entanglement (in units of log 3):

| cut | meaning | entropy |
|---|---|---|
| {A_in B_in \| A_out B_out} | **time** | 2 always (unitarity) |
| {A_in A_out \| B_in B_out} | **space** (the gate's entangling power) | rank S_BA |
| {A_in B_out \| B_in A_out} | **space–time diagonal** (A's past with B's future) | rank S_AA |

These formulas were checked against brute-force Choi entropies of random Clifford unitaries.

**Results.**
* **Complementarity law: det S_AA + det S_BA = 1 (mod 3)** for every gate. This is just the symplectic condition on the
  first block column. So the space entanglement and the diagonal entanglement can never both be deficient. The census:
  local gates (0,2): 576; swap-type (2,0): 576; (1,2): 18 432; (2,1): 18 432; perfect (2,2): 13 824.
* **Perfect gates** (Choi state = AME(4,3), maximally entangled across every cut) are exactly those with
  det S_AA = det S_BA = −1. There are **13 824 = 24³**, i.e. 4/15 of Sp(4,3). All four of their blocks are invertible,
  so their Choi Lagrangian is an MDS code of length 4.
* **One gate up to local operations.** All 13 824 form a single orbit under local Clifford operations on both sides,
  with a local stabilizer of order 24. This is the Clifford-level instance of the known uniqueness of the four-qutrit AME
  state up to local unitaries (Tan, arXiv:2601.19677).
* **The tetracode.** For time-reversing (anti-unitary) Clifford operations the law reads det S_AA + det S_BA = −1, and
  again 13 824 are perfect. The **tetracode gate** K⊗I, with K = [[1,1],[1,−1]] the tetracode's redundancy matrix, is one
  of them. Composed with local time reversal it is the unitary perfect gate.

**Reading: entanglement that is spatial, temporal and diagonal at once.** A perfect gate hides everything in pairs: every
two-qutrit marginal of its Choi state is maximally mixed. A's past is uncorrelated with A's future and with B's future
separately, and completely correlated with them jointly. The tripartite information I₃(A_in : A_out : B_out) = −2, the
minimum possible: maximal scrambling. Qubits have no such gate (AME(4,2) does not exist). For qutrits it exists, and in
W(3,3) it is the tetracode.

Prior art: AME(4,3) from the tetracode / orthogonal array OA(9,4,3,2) (Goyeneche–Życzkowski; Helwig et al.); uniqueness
up to local unitaries (Tan 2026). New here and specific to W(3,3):
* the complementarity law;
* the 24³ census across Sp(4,3);
* the unitary/anti-unitary split with the tetracode as the time-reversing representative;
* the space / time / diagonal reading, which ties back to the repo's "every point of W(3,3) carries a tetracode".

Related in the corpus:
* Passes 10946, 10954 and 10955 use the tetracode as a clock code, and the determinant character of GL(2,3)
  (det = −1 swaps S±) as a CP/Pin grading.
* Here the determinant appears differently, as the space/diagonal complementarity law det S_AA + det S_BA = 1 on
  two-qutrit gates.
* No earlier repo file treats AME(4,3) or perfect two-qutrit gates. The rediscovery hook's "clifford+tetracode" and
  "qutrit+tetracode" matches were read and are clock-code results, not this one.

