# Pass 11177 — a two-qutrit gate is perfect exactly when it moves the tensor factorisation to a collinear point of GQ(4,2): a tritangent plane sharing a line

Producer: `analysis/w33_pass11177_perfect_gates_tritangent.py`
GAP: `analysis/gap/w33_pass11177_factorisation_orbitals.g`, frozen as `data/w33_pass11177_gap_orbitals.txt`
Regression: `tests/test_w33_pass11177_perfect_gates_tritangent.py`

**Known ingredients** (Holotrade track, in the memory dictionary):
* A tensor factorisation C⁹ = C³⊗C³ is a pair of orthogonal nondegenerate planes {L, L⊥} of F₃⁴: an "octet" of 8
  points of W(3,3). There are 45.
* A complete factorisation frame is 5 factorisations whose octets partition the 40 points. There are 27: the lines of
  the cubic surface, GQ(2,4). The 45 factorisations are its tritangent planes.

**New (exact, all 51 840 elements of Sp(4,3)).**
* The octet-disjointness graph is **SRG(45,12,3,3)**, the collinearity graph of GQ(4,2). Every disjoint pair lies in
  exactly one of the 27 frames, and each factorisation lies in 3 frames.
* **The Choi entanglement of a gate is the relation it induces:**

| rank (S_AA, S_BA) | I₃ (Pass 11161) | S(F₀) relative to F₀ | gates |
|---|---|---|---|
| (2, 2) — perfect | −2 | **collinear** (octets disjoint: they share a frame) | 13 824 |
| (1, 2) or (2, 1) | −1 | non-collinear (octets meet) | 36 864 |
| (2, 0) or (0, 2) — local or swap | 0 | equal | 1 152 |

* The stabiliser H of F₀ (order 1152) is transitive on the 12 collinear factorisations, with point stabiliser 96. So
  **the perfect gates form one double coset HgH.**
* GAP: Sp(4,3) acts on the 45 with **rank 3, subdegrees 1, 12, 32**. This is the rank-3 action of U₄(2) on the 45
  tritangent planes; its index-45 subgroup is unique up to conjugacy.

**Reading.**
* How a two-qutrit register splits into "A" and "B" is a tritangent plane.
* A perfectly scrambling tick moves the split to a neighbouring tritangent plane, one sharing a line with it.
* A merely entangling tick moves it to a non-neighbour.
* The E₆ cubic has one monomial per tritangent plane, so a perfect tick joins two monomials that share a variable.
* This joins the space-time entanglement census (Passes 11156, 11161) to the W(E₆) geometry that runs through the
  repository.

**Prior art.** The factorisations, frames and GQ(2,4) structure are the Holotrade dictionary (commits 5419c27,
e92a047), and the U₄(2) rank-3 action is classical. We found no statement linking perfect (AME) gates to GQ(4,2)
collinearity.

**Related corpus files (read).**
* `2026-09-01_regulus_e8_completion_bridge.md` and `BT1792_BT1794_index_h27_execution.md` state the classical
  45/27 GQ(4,2) tritangent geometry; neither involves gates.
* `2026-07-15_pass84_e6_w33_explicit_iso.md` pairs the 40 points of W(3,3) with "40 tritangent planes". A smooth cubic
  surface has **45** tritangent planes: they are the 45 factorisations here. The 40 points correspond to a different
  W(E₆)-set. That count looks like a slip, flagged here rather than edited.
