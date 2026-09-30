# Pass 11181 — which qutrit dynamics entangle in every subsystem split: local dynamics only knows the primes 2 and 3

Producer: `analysis/w33_pass11181_intrinsic_entanglement.py`
GAP: `analysis/gap/w33_pass11180_mereology_classes.g` (all classes of PSp(4,3) and PSp(6,3) acting on the
factorisations)
Regression: `tests/test_w33_pass11181_intrinsic_entanglement.py`

**Exact fractions** of ticks that are local in 0, exactly 1, or at least 2 subsystem splits:

| register | splits | classes | local in 0 (intrinsically entangling) | exactly 1 | ≥ 2 |
|---|---|---|---|---|---|
| two qutrits | 45 | 20 | **19/45** (42.2%) | **19/48** | 131/720 |
| three qutrits | 110 565 | 74 | **7922/12285** (64.5%) | **7/48** | 41143/196560 |

Burnside check: the average number of fixed splits is exactly 1 in both cases, as it must be for a transitive action.
As the register grows, intrinsic entanglement becomes the typical case and dynamically selected subsystems become rare.

**The prime rule.** A tick is local in F iff it lies in F's stabiliser (SL(2,3)ⁿ ⋊ Sₙ)/⟨−1⟩. That stabiliser is a
{2,3}-group, of order |PSp|/#splits: 25920/45 = 2⁶3² for two qutrits and 4585351680/110565 = 2⁹3⁴ for three. Its
element orders are exactly the orders of the classes with fixed points. Hence:
* **every tick whose period has a prime factor ≥ 5 is intrinsically entangling**: orders 5, 7, 10, 13, 14, 15, 20, 30;
* for **two qutrits** the stabiliser's Sylow-3 subgroup is C₃×C₃, of exponent 3, so every order-9 tick is intrinsically
  entangling too. The fixed-point-free orders are exactly {5, 9};
* for **three qutrits** the Sylow-3 subgroup is C₃ wr C₃, of exponent 9. So order 9 splits: 6 of the 8 order-9 classes
  are fixed-point-free, and 2 fix 3 splits each. The fixed-point-free orders are {5, 7, 9, 10, 13, 14, 15, 18, 20, 30, 36}.
* Every order that does admit a local split is a {2,3}-number.

**Reading.** On the W(3,3) substrate, whether a dynamics can be described without generating entanglement, in *some*
choice of subsystems, is an arithmetic property of its period. Local qutrit dynamics only knows the primes 2 and 3.
Any tick whose period involves a larger prime entangles, whichever way the register is carved.
