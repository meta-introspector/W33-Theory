# Pass 11511 — the good rules for every n, reduced to the unipotent part of M at eigenvalue 1

Producer: `analysis/w33_pass11511_good_rules_all_n.py`
Certificate: `data/w33_pass11511_good_rules_all_n.json`
Regression: `tests/test_w33_pass11511_11515.py`

**Setting.**
* The good rule says: if Mz₁ = −z₁ or M²z₁ = −z₁, every frame is reversible.
* At k = 0, a reversal is an anti-symplectic A = QJ with AMA⁻¹ = M⁻¹ (equivalently A = MAM) and Az₁ = −z₁.
* Pass 11499 showed that a frame a is reversible through such an A whenever ω((A⁻¹ − I)x, a) = 0 for all
  x ∈ K = ker(M − I).
* At n = 2 every good class was covered. At n = 3 about 4% of the good mass stayed open.

## The reduction

**Lemma 1 (decoupling).** Let V₁ = ker(M − I)^{2n} be the generalised 1-eigenspace of M, and V′ the sum of its other
generalised eigenspaces.
* V₁ and V′ are M-invariant and ω-orthogonal.
* K lies in V₁. z₁ lies in V′, since its eigenvalue is −1 or a root of x² + 1.
* A reverser maps the generalised λ-eigenspace to the λ⁻¹-eigenspace, so A = A₁ ⊕ A′.
* Conversely, for **any** reverser A₁ of M₁ = M|_{V₁} and the V′-block A′ of one k = 0 solution, A₁ ⊕ A′ is again a
  k = 0 solution: it is anti-symplectic, satisfies A = MAM and Az₁ = −z₁.
* The criterion involves only A₁ and the V₁-component of a.

So the good rule holds for M as soon as two conditions hold:
* **(E)** (S) has a k = 0 solution;
* **(U)** for every a₁ ∈ V₁, some reverser A₁ of M₁ satisfies ω((A₁⁻¹ − I)x, a₁) = 0 for all x ∈ ker(M₁ − I).

**Lemma 2 (the semisimple case).** If M₁ = I, (U) holds in any dimension.
* Take the anti-symplectic involution that is −1 on a Lagrangian containing a₁ and +1 on a complementary Lagrangian.
* Its Im(A₁⁻¹ − I) is that Lagrangian, which is isotropic and contains a₁, so it is orthogonal to a₁.

**Exhaustion.** (U) holds for **every** unipotent element of Sp(2, 3) and of Sp(4, 3). There are 3² and 3⁸ = 6561 of them
(Steinberg's count), each checked against all 51 840 anti-symplectic maps. All pass (`exhaustion_U`). The minimum number of reversers is 6 on Sp(2,3) and 18 on Sp(4,3).

**Lemma 3 ((E) at block level).**
* Refine the split as V = V₁ ⊕ V_z ⊕ V_rest, where V_z is the generalised block of z₁:
  * ker(M + I)^{2n} if Mz₁ = −z₁;
  * ker(M² + I)^{2n} if M²z₁ = −z₁.
* All three pieces are M-invariant, nondegenerate and ω-orthogonal, and a reverser preserves each of them.
* By Wonenburger (characteristic ≠ 2), M is reversed by an anti-symplectic involution on V₁ and on V_rest.
* So (E) holds as soon as **some reverser of M|V_z negates z₁**. For M²z₁ = −z₁ such a reverser then fixes Mz₁
  automatically, since AM = M⁻¹A.
* **Exhaustion.** This block statement holds for **every** element of Sp(2, 3) and Sp(4, 3) of the relevant type, and
  every admissible z:

| block | elements | vectors z | failures |
|---|---|---|---|
| Sp(2,3), M + I nilpotent, Mz = −z | 9 | 24 | 0 |
| Sp(2,3), M² + I nilpotent, M²z = −z | 6 | 48 | 0 |
| Sp(4,3), M + I nilpotent, Mz = −z | 6 561 | 19 440 | **0** |
| Sp(4,3), M² + I nilpotent, M²z = −z | 4 860 | 77 760 | **0** |

> **Theorem.** For every n: if M z₁ = −z₁ or M²z₁ = −z₁, and
> * the unipotent part at 1 has dim V₁ ≤ 4, or M is semisimple at eigenvalue 1, and
> * the block of z₁ has dimension ≤ 4,
>
> then every frame is reversible.

No per-class search remains: (U) is covered by Lemma 2 and the unipotent exhaustion, (E) by Lemma 3.

## Coverage

**n = 2 (all 2592 good classes).**
* The block theorem applies to **every one**.
* The decider confirms all frames reversible on every class it settles: 2591. The exception is M = −I, which lies beyond
  its cap.

**n = 3 (Pass 11373's orbits in the good cells).**

| route | Mz₁ = −z₁ (mass) | M²z₁ = −z₁ (mass) |
|---|---|---|
| block theorem (dim V₁ ≤ 4 or semisimple, z₁-block ≤ 4) | 93 orbits, 87.3% | 82 orbits, 85.4% |
| z₁-block = whole space (dim 6), (E) shown per class | 33 orbits, 12.7% | 9 orbits, 14.6% |
| **uncovered** | **0** | **0** |

* **The per-class (E) cases.** Here V₁ = 0, so (U) is vacuous, and a k = 0 solution was exhibited:
  * Q = I for 5 orbits;
  * the single-map hypothesis for 19;
  * random sampling of the solution space for 18.
* So **every good orbit at n = 3 is proved**.
* Pass 11499 independently confirmed the decided n = 3 good orbits against the decider.

**What remains for a proof at every n.**
* **(E) for z₁-blocks of dimension ≥ 6.** The block statement of Lemma 3 needs to hold in larger blocks.
* **(U) for unipotent parts of dimension ≥ 6.** Sp(6, 3) has 3¹⁸ unipotent elements, too many for exhaustion. It needs
  the structure theory of reversers (Wonenburger: every symplectic matrix is a product of two anti-symplectic
  involutions).

**Prior art.**
* The rules are Pass 11373's, and the frame-level criterion is Pass 11499's.
* Reversibility of symplectic matrices by anti-symplectic involutions is classical (Wonenburger 1966).
* The decoupling lemma and the exhaustions are new here.
