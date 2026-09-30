# Pass 11188 — the intrinsic arrow, exactly: at most one trit per qutrit, and maximal exactly when no subsystem is left in place

Producer: `analysis/w33_pass11188_exact_arrow_by_class.py`
GAP: `analysis/gap/w33_pass11188_class_reps.g` (n = 2, 3), `analysis/gap/w33_pass11188_class_reps_n4.g` (n = 4)
Frozen: `data/w33_pass11188_gap_class_reps*.txt`, `data/w33_pass11188_exact_arrow_by_class.json`
Regression: `tests/test_w33_pass11188_exact_arrow_by_class.py`

**Setup.** Pass 11183 defined the intrinsic arrow A(S) of a Clifford tick as the minimum, over all subsystem splits F, of
the information E(S, F) = Σ_q export_q exported per tick. It computed A exactly for two qutrits and on a 240-element
sample for three. A is a class function on PSp(2n, 3): the splits form one Sp-orbit, and the central sign changes no
rank. GAP prints one matrix per conjugacy class in our adapted basis. Each matrix is checked: symplectic, with GAP's
projective order and GAP's count of fixed splits (Pass 11180).

**The plane formula.** export_q = 2 − d(P_q), where d(P) = dim(P ∩ SP). Hence

    A(S) = 2n − max { Σ d(P) : P₁, …, P_k mutually orthogonal nondegenerate planes }.

Any orthogonal family extends to a split, and planes with d = 0 add nothing. This turns a minimum over 110 565 splits
(n = 3), or about 8.3·10⁹ splits (n = 4), into a branch and bound over the few thousand planes that meet their image. It
agrees with direct enumeration on all 20 + 74 classes for n = 2, 3.

**Consequences for every n.**
* d(P) = 2 iff SP = P.
* E ≤ n − 1 forces Σ d ≥ n + 1, so some plane is invariant; therefore A(S) ≤ n − 1 ⇒ S fixes a nondegenerate plane.
* E = 1 is impossible: all planes but one invariant forces the last one to be invariant too. So **A ≠ 1** for every n.
* Conversely, if A ≤ n − 1 for all (n−1)-qutrit ticks, a fixed plane P gives A(S) ≤ A(S|P^⊥) ≤ n − 1.

**Results (exact, all conjugacy classes).**

| n | classes | values of A | fraction of PSp(2n, 3) with A = 0 / 2 / 3 / 4 |
|---|---|---|---|
| 2 | 20 | {0, 2} | 133/360, 227/360 (= Pass 11183) |
| 3 | 74 | {0, 2, 3} | 55241/884520, 581/1080, 8836/22113 |
| 4 | 278 | {0, 2, 3, 4} | 43815911/8703676800, 1769071/10497600, 219917/597051, 43815421/95644800 |

* **A ≤ n for n = 2, 3, 4.** No tick forces the loss of more than one trit per qutrit. At most half of the 2n trits a
  split could export are unavoidable.
* **A = n ⇔ S fixes no nondegenerate plane**, checked class by class for n = 2, 3, 4. The maximal arrow belongs exactly
  to the ticks that leave no qutrit, in any description, in place. For those ticks there is always a *half-moving split*:
  every subsystem meets its own image in exactly one line.
* Orders with A = n: {4, 5, 6, 9} for n = 2; {3, 6, 7, 9, 13, 14} for n = 3, which includes every tick whose order is
  divisible by 7 or 13; for n = 4 the orders include 41 (the irreducible ticks, 41 | 3⁴ + 1).

**Conjecture (all n).** A(S) ≤ n, with equality iff S has no invariant nondegenerate plane. By the argument above,
this is equivalent to: every tick with no invariant nondegenerate plane has a half-moving split.

**Scope.** This is exact for n ≤ 4. The GAP class lists are complete: class sizes sum to |PSp(2n, 3)|. The plane counts
are 90, 7371 and 597 780. The branch and bound is exact, not a heuristic.
