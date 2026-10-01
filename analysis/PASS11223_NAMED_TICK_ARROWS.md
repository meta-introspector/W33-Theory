# Pass 11223: the arrow of the substrate's named ticks

Producer: `analysis/w33_pass11223_named_tick_arrows.py`
Certificate: `data/w33_pass11223_named_tick_arrows.json`
Regression: `tests/test_w33_pass11223_named_tick_arrows.py`

**Background.**
* The arrow A(S) = n − c(S) (Passes 11207 and 11217) does not depend on a factorisation.
* A = 0 exactly when S is a product of single-qutrit gates in **some** tensor factorisation.
* Pass 11193 classified two-qutrit Cliffords relative to the **standard** factorisation: local (1152), perfect (13 824),
  partial (36 864). This pass crosses the two classifications.

## Two qutrits: all 51 840 elements of Sp(4,3)

| relative to the standard split | A = 0 | A = 2 | P(arrow-free) |
|---|---:|---:|---|
| local | 816 | 336 | 17/24 |
| perfect | 4 896 | 8 928 | **17/48** |
| partial | 13 440 | 23 424 | 35/96 |
| all | 19 152 | 32 688 | 133/360 |

The total arrow-free fraction is 133/360, Pass 11204's value. All 576 non-swapping local gates are arrow-free. Of the
576 swapping local gates, 240 are arrow-free and 336 have A = 2.
* **A third of the maximally entangling gates have no intrinsic arrow.** They are perfect for the standard split and
  products in another one.
* **Exchanges are arrows.** Some local gates have A = 2. h\* is one: up to local gates it swaps the two qutrits, so no
  plane of the standard split is invariant. The bare SWAP has A = 0, because it is local in the symmetric/antisymmetric
  split.

**Named two-qutrit ticks.**
* Every one of the 80 transvections (the W33 point symmetries) has A = 0.
* SUM has A = 0.
* **The fixed perfect gate p of Pass 11193 has A = 0**, and so does p².
* The local word h\* (order 12) has A = 2, and so does p h\* p.

## Three qutrits: Pass 11182's paper ticks

| tick | A | c | order |
|---|---:|---:|---:|
| clock K, K², dual kick V | 0 | 3 | 3 |
| perfect interacting tick V K V | **2** | 1 | 9 |
| F₉ gate | **3** | 0 | 3 |

* V K V is local in no split (Pass 11182). Still, it keeps one qutrit dynamically protected (c = 1).
* The F₉ gate is "local" in 27 splits, mapping planes onto planes, yet it has no invariant plane. So it cycles the three
  qutrits there, and carries the maximal arrow.

## Reading

Whether a gate is maximally entangling depends on both the gate and the chosen factorisation. The arrow is what remains
once the factorisation is chosen optimally. The paper's perfect two-qutrit gate carries no intrinsic arrow. Its
three-qutrit perfect tick carries two trits per tick. The largest arrow belongs to a gate that only relabels the qutrits.
