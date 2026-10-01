# Pass 11268 — canonical perfect magic ticks: one magic leg never breaks the arrow; all but one always do

Producer: `analysis/w33_pass11268_perfect_magic_ticks.py`
Certificate: `data/w33_pass11268_perfect_magic_ticks.json`
Regression: `tests/test_w33_pass11265_11269.py`

**Setup.**
* Perfect (maximally scrambling) qutrit ticks exist for 1, 2, 3 and 5 qutrits (Pass 11170).
* Perfectness is invariant under local unitaries on every leg. So dressing a perfect Clifford tick V with cubic phases,
  U(a, b) = (⊗ T^{a_i}) V (⊗ T^{b_i}) with a_i, b_i ∈ Z₉, gives perfect **magic** ticks. A leg carries magic iff its
  exponent is not 0 mod 3, since T³ = Z is Clifford.
* The perfect Clifford ticks are built from Pass 11170's perfect F₉-circulants: circulant (4,5) for n = 2 and (1,1,4)
  for n = 3. Each is embedded as a symplectic matrix, realised by the Weil twirl, and verified perfect on every n|n cut
  (also rechecked after dressing).
* Every dressing is decided exactly with Pass 11252's criterion.

## Two qutrits: all 9⁴ = 6561 dressings

| legs carrying magic | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| violating / total | 0/81 | **0/648** | 216/1944 = **1/9** | **2592/2592** | 648/1296 = **1/2** |

## Three qutrits

| legs carrying magic | 1 | 5 |
|---|---|---|
| violating / total, exhaustive | **0/8748** | **139968/139968** |

A random sample of 1500 dressings fills in the middle: 2 legs 69/120, 3 legs 311/326, 4 legs 435/499, 6 legs 86/132.

## Pattern (computer-verified for the fixed n = 2 and n = 3 circulant representatives)

* **Magic on exactly one leg of either tested perfect tick never breaks substrate time reversal.**
* **Magic on all but one leg (2n − 1) always does for either tested tick.**
* With magic on every leg the outcome is mixed: exactly 1/2 for n = 2.

**Comparison.**
* Generic Clifford+T ticks break the arrow with probability 8–10% for a single magic gate (Pass 11266).
* A perfect scrambler with one magic leg never does.
* With 2n − 1 magic legs it always does.

For these representatives, maximal scrambling makes the arrow an all-or-nothing function at the two edge weights.

**Open.** A proof for all perfect ticks (n = 5 is the next case), and the meaning of the exact 1/9 and 1/2 at even
numbers of magic legs.
