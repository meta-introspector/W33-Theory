# Pass 11311 — the perfect-tick magic edge law tested on the whole Clifford group: it fails except on the F₉-linear class

Producer: `analysis/w33_pass11311_perfect_edge_law.py`
Certificate: `data/w33_pass11311_perfect_edge_law.json`
Regression: `tests/test_w33_pass11308_11312.py`

**Claim tested.** Pass 11268 found, for two fixed perfect circulant representatives, that:
* magic on exactly one leg never breaks substrate time reversal;
* magic on 2n − 1 legs always does.

Codex's audit scoped this to those representatives. Here it is tested over every two-qutrit Clifford tick, split by
the compiler classes of Pass 11193 (perfectness verified on both nontrivial 2|2 leg cuts).

**Setup.**
* Random ticks V = P_a·V_M, with the Pauli frame uniform and the symplectic class uniform within its class.
* Each tick is dressed with cubic phases on k of its 4 legs; the other legs carry Z^j Paulis.
* Every dressing is decided exactly. There are 6000 samples per cell.

| class (symplectic classes) | k = 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| local incl. swap (1152) | 0 | 0.081 | 0.215 | 0.435 | 0.314 |
| perfect, F₉-linear (64) | 0 | 0.086 | 0.431 | **1.000 (6000/6000)** | 0.630 |
| perfect, other (13760) | 0 | 0.068 | 0.533 | 0.906 | 0.779 |
| other entangling (36864) | 0 | 0.099 | 0.508 | 0.847 | 0.733 |

**What survives and what falls.**
* **"One magic leg never" is false** for perfect ticks in general: 6.8% of them violate, and 8.6% on the F₉-linear
  class. It held only for the two representatives of Pass 11268.
* **"2n − 1 magic legs always" holds on the F₉-linear perfect class**, the class of Pass 11170's circulants (6000/6000,
  Pauli frames included), and **fails** for other perfect ticks (90.6%).
* A k = 3 tick violates at least 85% of the time in every class.

So the clean statement is about F₉-linearity, not about perfectness. Perfect but non-F₉-linear ticks behave like
generic entangling ones.

**Open.** A proof of the F₉-linear k = 2n − 1 law, and its n = 3 analogue (Pass 11268 verified 139 968/139 968 for one
F₉-linear three-qutrit representative).
