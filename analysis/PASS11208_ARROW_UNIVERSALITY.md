# Pass 11208 — the arrow law is universal (qubits, qutrits, 5- and 7-dits), and a large system loses one trit per qutrit less ≈ 0.7255

Producer: `analysis/w33_pass11208_arrow_universality.py`
GAP: `analysis/gap/w33_pass11208_class_reps_q.g`, `analysis/gap/w33_pass11208_class_reps_q_big.g`,
`analysis/gap/w33_pass11207_class_reps_n5.g` (PSp(10,3): 940 classes, 115 min)
Frozen: `data/w33_pass11208_arrow_universality.json` and the three GAP class lists in `data/`
Regression: `tests/test_w33_pass11208_arrow_universality.py`

**Numbering.** First committed on the branch as Pass 11200; renumbered after master's reservation of 11199–11206 (see
Pass 11207). Master's reserved item "a second law from counting (arrow-free fraction vs n)" overlaps §2. The exact
arrow-free fractions are P(c = n): 133/360, 55241/884520 and 43815911/8703676800 for n = 2, 3, 4, and 0.0002 for n = 5.

**Background.** Pass 11207 proved the law A(S) = n − c(S) for every qutrit Clifford tick. Here A is the intrinsic arrow
and c the number of protected qutrits. The proof uses only odd characteristic.

## 1. Universality: exact, on every conjugacy class

A is computed exactly by branch and bound over planes, independently of the theorem, and compared with n − c:

| group | local dimension | classes | law A = n − c | fraction with maximal arrow A = n |
|---|---|---|---|---|
| Sp(4,2), Sp(6,2), Sp(8,2), Sp(10,2) | qubits, n ≤ 5 | 11, 30, 81, 198 | holds on all | 51/80, 953/2016, 3366429/6963200, 2319469541/4693032960 |
| PSp(4,5) | 5-dits | 34 | holds on all | 2317/3900 |
| PSp(4,7) | 7-dits | 52 | holds on all | 33487/58800 |
| PSp(10,3) | five qutrits | 940 | holds on all (see below) | 0.45029 |

* **Qubits are not covered by the proof.** The Cayley transform and the diagonalisation of symmetric forms fail in
  characteristic 2. The law nevertheless holds on every class through five qubits.
* The odd-characteristic mechanism is not what makes it work there. Over F₂, the blocks V(4) and W(2) have transverse
  bi-Lagrangians, but b(x, y) = ω(x, Sy) is **alternating** on every one of them, so the radical lemma cannot diagonalise
  it. For V(6), V(8) and W(4), non-alternating ones exist.
* The {f, f*} components do carry over to characteristic 2. If every admissible Φ·r(A) were alternating, A would be
  self-adjoint for a symplectic form, which forces paired elementary divisors; that is impossible on a cyclic
  decomposition.
* **Conjecture:** the law holds for qubits for every n.

**Five qutrits, all 940 classes.** For each class:
* take the nondegenerate parts of the ±1-eigenspaces (every plane there is invariant);
* add a maximum family of invariant planes of type span(x, Sx);
* find a half-moving split of the rest, allowing the last plane to pass through an eigenvector, as the W(k) mechanism
  requires.

This gives a split with E = 5 − c, which meets the lower bound, so A = 5 − c exactly. In every class c equals the Jordan
closed form. The exact distribution over PSp(10,3) is:

| A | 5 | 4 | 3 | 2 | 0 |
|---|---|---|---|---|---|
| fraction | 0.45029 | 0.39905 | 0.12589 | 0.02456 | 0.00020 |

**Burnside.** Sp(2n, q) is transitive on nondegenerate planes. So the **average number of planes a tick leaves invariant
is exactly 1**, for every n and q. This is confirmed from the class data in every case above and for qutrits n ≤ 4.

## 2. The law for large systems (qutrits)

The expected number of protected qutrits is nearly independent of n. Exactly:

| n | 2 | 3 | 4 | 5 |
|---|---|---|---|---|
| E[c] | 133/180 = 0.73889 | 106927/147420 = 0.725322 | 0.725518 | 0.725535 |

Random ticks up to n = 32, using the polynomial-time closed form of Pass 11207, give P(c = 0) ≈ 0.45, P(c = 1) ≈ 0.39,
P(c = 2) ≈ 0.13, P(c = 3) ≈ 0.02 and E[c] ≈ 0.72 (values frozen in the JSON). So

    E[A] = n − E[c] = n − 0.7255… ,

and a typical large Clifford dynamics loses exactly one trit per qutrit per tick, less a bounded number of protected
qutrits (fewer than one on average). About 45% of ticks have the maximal arrow A = n.

**Reading.** Taken with the "half the record" reading of Pass 11207 (each unprotected qutrit exports exactly one of its
two trits in the best description), the unavoidable dissipation of a generic Clifford dynamics is k_BT ln 3 per qutrit
per tick, up to O(1). The law is not a special property of qutrits. The W(3,3) substrate inherits it as a property of
Clifford dynamics.

**Scope.**
* Exact class-by-class for the groups listed.
* Large-n statistics are Monte Carlo over products of random transvections, with standard errors of about 0.01.
* The limit E[c] → 0.7255… is observed, not proved. A closed form would follow from cycle-index methods (Fulman)
  applied to the Jordan formula for c.
* PSp(6,5) was started and dropped: its 406 875-plane enumeration is too slow in pure Python, and odd q is covered by
  the proof.
