# Pass 11184 — the entanglement data classify how two three-qutrit subsystem splits relate, up to the direction of time

Producer: `analysis/w33_pass11184_merged_orbitals.py`
Regression: `tests/test_w33_pass11184_merged_orbitals.py`

**Independent recomputation.** The stabiliser of the standard split (local SL(2,3) and qutrit relabellings) acts on
the 110 565 splits. Union-find over its generators gives **20 orbits with exactly GAP's sizes** (Pass 11178). This is a
second source for the rank-20 geometry.

**A finer invariant.** For a split F, the block of F's plane Qⱼ projected onto the standard plane Pᵢ has, when it has
rank one, an image point in Pᵢ and a kernel point in Qⱼ. The invariant records three things:
* for each row, the number of distinct image points;
* for each column, the number of distinct kernel points;
* an incidence count linking kernels to images.

It is constant on orbits and separates more than the signature did:
* the pair {6912, 6912} splits;
* the five-orbital class {256, 256, 2304, 6912, 6912} splits into {256, 256}, {2304} and {6912, 6912}.

**Direction of time.** For each orbital, the reversed relation (F, F₀), which is the inverse gate, was located:
* **14 relations are symmetric**;
* **6 form three pairs of mutually reversed relations**, of sizes 256, 2304 and 6912;
* the two ties the fine invariant leaves are **exactly two of these reversal pairs**. The third pair (2304) is already
  split by the signature, because its two orbitals have transposed signatures.

**Result.** Block ranks, orientations and image/kernel incidences classify the relation between two three-qutrit splits
completely, **except for its direction in time**. For two pairs of relations, a gate and its inverse look identical to
every one of these entanglement measures, and only the relation's geometry records which way it runs.
