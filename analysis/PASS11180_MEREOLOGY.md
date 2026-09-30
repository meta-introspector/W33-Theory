# Pass 11180 — quantum mereology on W(3,3): the dynamics decides the subsystems, and 42% of two-qutrit ticks entangle in every one

Producer: `analysis/w33_pass11180_mereology.py`
GAP: `analysis/gap/w33_pass11180_mereology_classes.g`, frozen as `data/w33_pass11180_gap_classes.txt`
Regression: `tests/test_w33_pass11180_mereology.py`

**Idea.** A theory of everything has no given subsystems. Which tensor factorisation counts as "physical" must come
from the dynamics: this is quantum mereology (Zanardi 2001; Zanardi–Lidar–Lloyd 2004; Carroll–Singh, PRA 103, 022213
(2021)). On the W(3,3) substrate the candidate factorisations are finite. For two qutrits they are the 45 tritangent
planes (Pass 11177), so the question can be answered exactly.

**Method.** For a tick S and a factorisation F with adapted symplectic basis B_F, S_F = B_F⁻¹ S B_F is "S as seen by
F's subsystems":
* S is **local** in F iff S_F is block-monomial;
* S is **perfect** in F iff every block of S_F is invertible;
* otherwise S is **partial** in F.

This reproduces the octet computation of Pass 11177 exactly.

**Census** (all 51 840 elements of Sp(4,3)), as (#local, #perfect, #partial) over the 45 splits:

| profile | projective order | ticks |
|---|---|---|
| (45, 0, 0) | 1 | 2 |
| (13, 0, 32) | 2 | 90 |
| (5, 24, 16) | 2 | 540 |
| (9, 0, 36) | 3 | 160 |
| (3, 0, 42) | 3 | 960 |
| **(6, 27, 12)** | 3 | 480 |
| (1, 12, 32) | 4, 12 | 16 200 |
| (1, 0, 44), (1, 18, 26), (2, 15, 28), (4, 9, 32) | 6 | 1440, 2880, 4320, 2880 |
| **(0, 15, 30)** | **5** | 10 368 |
| **(0, 9, 36)** | **9** | 11 520 |

**Readings.**
* **Intrinsically entangling:** local in no split, exactly when the projective order is 5 or 9. That is 19/45 of all
  ticks (GAP, exact). No choice of subsystems makes them non-entangling.
* **Unique emergent subsystems:** local in exactly one split. That is 19/48 of all ticks (orders 4, 6, 12, including
  the E₆ Coxeter order 12). These ticks single out their own two-qutrit structure.
* **Perfection is relative:** 480 ticks are product gates for 6 subsystem splits and maximal scramblers for 27.
* A side check: perfect gates are spread evenly over orders. The order-12 fraction is 4/15, the global average, so
  they have no special tie to the E₆ Coxeter class.

**Prior art.** That entanglement generation depends on the tensor product structure is known (Zanardi, Lidar, Lloyd;
Carroll–Singh; "Tensor Product Structure Geometry under Unitary Channels", Quantum 2025). New here: the exact finite
census on the W(3,3) substrate, with tritangent planes as the candidate subsystems.
