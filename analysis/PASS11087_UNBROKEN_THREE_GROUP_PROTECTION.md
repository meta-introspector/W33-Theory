# Pass 11087 — what protects the massless fractional charges: an interlocking 3-group, not one symmetry

Producer: `analysis/w33_pass11087_unbroken_three_group_protection.py`
Certificate: `data/w33_pass11087_unbroken_three_group_protection.json`
Regression: `tests/test_w33_pass11087_unbroken_three_group_protection.py`
Builds on: Pass 11024. Solvers: `w33_exact_monomial_orders.py` (lattice test) and `unbroken_group`.

In every Z3×Z3 parity vacuum with the hidden SU(4) condensed, Pass 11024 found at least 120 fractionally
charged multiplets massless at all orders. Most are protected by an exact unbroken discrete 3-group of
order 3⁶ or 3⁴. This pass takes that group apart. The hypothesis to test: the protection is simply the
orbifold's space-group selection rules surviving the vacuum.

**Method.** The charge coordinates fall into blocks:
* **U1**: the continuous U(1)s; the vacuum breaks them to discrete remnants.
* **PG**: the two point-group sector labels.
* **SG**: the three space-group translation charges, i.e. the fixed-point labels.
* **R**: the three Z₃^R plane rotations.
* **centre**: the N-ality of the hidden SU(4) and SU(3)′.

For every vacuum the pass computes two things. First, the unbroken group G = torsion(L_all / L_vac) and
its order when one block is dropped. Second, for every lattice-forbidden mass term in the
hidden-singlet charge-1/3 class, which blocks are *necessary*: dropping a necessary block makes the
term allowed.

## Result

**The group**

| vacua | unbroken U(1)s | G | order |
|---|---|---|---|
| 55 hidden-condensate (Pass 11024 D) | 1 (U(1)_Y) in 54; 5 in the c1 near miss | Z₃²×Z₉² (36), Z₃⁴×Z₉ (12), Z₃²×Z₉ (6), Z₃³ (1) | 3⁶ (48), 3⁴ (6), 3³ (1) |
| 28 singlet (Pass 10979) | **2**: U(1)_Y and an extra U(1)′ | Z₃⁴ (16), Z₃²×Z₉ (12) | 3⁴ |

Dropping one block (typical 3⁶ vacuum):

| block dropped | effect on G |
|---|---|
| PG | unchanged |
| centre | unchanged |
| SG | 729 → 81 |
| R | 729 → 81 |
| U1 | 729 → 9–81 (the Z₉ factors are U(1) remnants) |

The c1 near-miss vacuum (order 3³, five unbroken U(1)s) is the exception: dropping SG leaves G unchanged,
dropping R trivialises it, and dropping U1 gives order 18, the factor 2 coming from the hidden SU(2)
centres. For the 54 SU(4)-condensed vacua, G is assembled from three sources: remnants of the broken U(1)s, the space-group translation Z₃s,
and the Z₃^R rotations. The point-group labels and the hidden centres add nothing independent.

**The protection** (81,876 forbidden charge-1/3 mass terms in the 55 condensate vacua)

| which blocks are necessary | terms | share |
|---|---|---|
| none: redundant, independent symmetries each forbid it | 38,190 | 47% |
| U(1) remnants alone | 21,206 | 26% |
| U(1) + R | 9,636 | 12% |
| space group alone | 8,256 | 10% |
| U(1) + space group | 3,794 | 5% |
| U(1) + space group + R | 742 | 1% |
| R alone | 52 | 0.1% |

In the 28 singlet vacua the pattern is similar. Of 40,820 forbidden terms, 57% are redundant, 21% need
the space group and 22% need the U(1) remnants, alone or with R or the space group.

## Reading

* **The hypothesis is refuted.** The protection is *not* "the space-group selection rules surviving".
  The space group is necessary for only 16% of the forbidden terms; the broken-U(1) remnants for 43%.
  Almost half are forbidden redundantly. No single symmetry of the orbifold is the switch that, if
  broken, would release the fractional charges. Pass 11090 quantifies what a smooth resolution could
  release.
* **A new defect of the singlet vacua.** All 28 hidden-unbroken vacua of Pass 10979 leave an extra
  unbroken U(1)′ besides hypercharge, i.e. a massless gauge boson beyond the Standard Model. **In all 28 it
  couples to Standard-Model matter**, since some quark or lepton carries nonzero U(1)′ charge. A
  gauge-strength long-range force on ordinary matter is excluded. That makes a second reason, independent of
  the fractional-charge index, to discard every hidden-unbroken Z3×Z3 vacuum. The U(1)′ is not the anomalous
  U(1): the FI VEVs break that one.

Scope: lattice-level (exact-symmetry) protection of the hidden-singlet charge-1/3 class. Holomorphy-only
obstructions are counted in Pass 11024.
